"""Reproducible HW7 CTLE design with matched-source voltage normalization.

Run ``python homework/hw07/src/hw7.py`` from any directory.  The input data
remain unchanged; plots and machine-readable results are written separately.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import irfft
from scipy.signal import oaconvolve
from scipy.special import i0

BIT_RATE = 40e9
RISE_TIME = 5e-12
SOURCE_VPP = 2.0
MATCHING_DIVIDER = 0.5


@dataclass(frozen=True)
class TouchstoneData:
    frequency_hz: np.ndarray
    s: np.ndarray
    reference_ohm: float


def read_touchstone(path: str | Path) -> TouchstoneData:
    """Read a two-port Touchstone 1.x file without resampling any data.

    HZ/KHZ/MHZ/GHZ and RI/MA/DB are accepted. Continuation lines and comments
    are allowed. Unsupported network types, Touchstone 2.0 sections, malformed
    records, nonfinite values, and nonincreasing frequencies are rejected.
    Matrix ordering in a two-port file is S11, S21, S12, S22.
    """
    tokens: list[float] = []
    option = None
    for number, raw in enumerate(Path(path).read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            raise ValueError("Only Touchstone 1.x two-port files are supported")
        if line.startswith("#"):
            if option is not None or tokens:
                raise ValueError("Exactly one option line before data is required")
            option = line[1:].upper().split()
            if len(option) != 5 or option[1] != "S" or option[3] != "R":
                raise ValueError("Expected '# <unit> S <RI|MA|DB> R <impedance>'")
            continue
        if option is None:
            raise ValueError("Missing Touchstone option line")
        try:
            tokens.extend(float(x.replace("D", "E").replace("d", "e")) for x in line.split())
        except ValueError as exc:
            raise ValueError(f"Invalid numeric data on line {number}") from exc
    if option is None or not tokens or len(tokens) % 9:
        raise ValueError("Incomplete two-port data: each record needs nine numbers")
    scales = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}
    if option[0] not in scales or option[2] not in {"RI", "MA", "DB"}:
        raise ValueError("Unsupported frequency unit or complex representation")
    reference = float(option[4])
    records = np.asarray(tokens).reshape(-1, 9)
    if not np.all(np.isfinite(records)) or not np.isfinite(reference) or reference <= 0:
        raise ValueError("Data and positive reference impedance must be finite")
    f = records[:, 0] * scales[option[0]]
    if np.any(f < 0) or np.any(np.diff(f) <= 0):
        raise ValueError("Frequencies must be nonnegative and strictly increasing")
    a, b = records[:, 1::2], records[:, 2::2]
    if option[2] == "RI":
        values = a + 1j * b
    else:
        if option[2] == "MA" and np.any(a < 0):
            raise ValueError("MA magnitudes cannot be negative")
        magnitude = 10.0 ** (a / 20) if option[2] == "DB" else a
        values = magnitude * np.exp(1j * np.deg2rad(b))
    s = np.empty((len(f), 2, 2), dtype=complex)
    s[:, 0, 0], s[:, 1, 0], s[:, 0, 1], s[:, 1, 1] = values.T
    return TouchstoneData(f, s, reference)


def ctle_response(f_hz, fp1_ghz=20.6, fp2_ghz=22.8, adc_db=-4.91, aac_db=11.48):
    """USB-form CTLE: Aac*p2*(s+Adc*p1/Aac)/((s+p1)*(s+p2))."""
    if fp1_ghz <= 0 or fp2_ghz <= 0:
        raise ValueError("Pole frequencies must be positive")
    f = np.asarray(f_hz, dtype=float) / 1e9
    adc, aac = 10 ** (adc_db / 20), 10 ** (aac_db / 20)
    return aac * fp2_ghz * (1j * f + adc * fp1_ghz / aac) / ((1j * f + fp1_ghz) * (1j * f + fp2_ghz))


def ctle_peak(fp1_ghz, fp2_ghz, adc_db, aac_db):
    """Analytic global maximum on 0 <= f < infinity, including the DC case."""
    if fp1_ghz <= 0 or fp2_ghz <= 0:
        raise ValueError("Pole frequencies must be positive")
    adc, aac = 10 ** (adc_db / 20), 10 ** (aac_db / 20)
    z = adc * fp1_ghz / aac
    z2 = z * z
    discriminant = z2 * z2 - z2 * (fp1_ghz ** 2 + fp2_ghz ** 2) + (fp1_ghz * fp2_ghz) ** 2
    x = max(0.0, -z2 + np.sqrt(max(0.0, discriminant)))
    f = float(np.sqrt(x) * 1e9)
    peak_db = float(20 * np.log10(abs(ctle_response(f, fp1_ghz, fp2_ghz, adc_db, aac_db))))
    return {"peak_gain_db": peak_db, "peak_frequency_hz": f, "zero_frequency_hz": float(z * 1e9)}


def kaiser_weights(f_hz, alpha=1.0):
    f = np.asarray(f_hz, dtype=float)
    if f.ndim != 1 or len(f) < 2 or f[0] != 0 or f[-1] <= 0 or alpha < 0:
        raise ValueError("Need nonnegative alpha and a frequency array from DC")
    return i0(np.pi * alpha * np.sqrt(np.maximum(0.0, 1 - (f / f[-1]) ** 2))) / i0(np.pi * alpha)


def discrete_impulse(f_hz, response, samples_per_ui=50, alpha=1.0, bit_rate=BIT_RATE):
    """Return a real discrete convolution kernel, whose sum equals H(0).

    irfft supplies both conjugate frequency halves and its 1/N normalization.
    H bins above the measured fmax are zero. Neither an extra factor of two nor
    a separate DC half-weight is required. No source divider is applied here.
    """
    f, response = np.asarray(f_hz, float), np.asarray(response, complex)
    if len(f) != len(response) or len(f) < 2 or f[0] != 0:
        raise ValueError("Response must include DC and match frequency length")
    df = f[1] - f[0]
    if df <= 0 or not np.allclose(np.diff(f), df, rtol=1e-10, atol=1e-4):
        raise ValueError("Frequency samples must be uniform")
    if samples_per_ui < 2 or bit_rate <= 0:
        raise ValueError("Need at least two samples/UI and positive bit rate")
    dt = 1 / (bit_rate * samples_per_ui)
    n = int(round(1 / (df * dt)))
    if n % 2 or not np.isclose(n * df * dt, 1, rtol=1e-10):
        raise ValueError("Time step must form an even integer FFT period")
    if len(f) > n // 2 + 1:
        raise ValueError("Requested time step violates measured Nyquist frequency")
    spectrum = np.zeros(n // 2 + 1, complex)
    spectrum[:len(f)] = response * kaiser_weights(f, alpha)
    if abs(spectrum[0].imag) > 1e-10:
        raise ValueError("DC response must be real")
    if len(f) == len(spectrum) and abs(spectrum[-1].imag) > 1e-10:
        raise ValueError("Nyquist response must be real")
    return np.arange(n) * dt, irfft(spectrum, n=n)


def make_nrz(bits, samples_per_ui=50, rise_time_s=RISE_TIME, source_vpp=SOURCE_VPP, bit_rate=BIT_RATE):
    """Construct 0..source_vpp NRZ; transitions include their initial endpoint.

    A ramp occupies rise_time_s with samples k/Nrise for k=0..Nrise-1.
    The source is low before the first bit. No matching division occurs here.
    """
    bits = np.asarray(bits, dtype=float)
    if bits.ndim != 1 or not len(bits) or np.any((bits != 0) & (bits != 1)):
        raise ValueError("bits must be a nonempty one-dimensional binary array")
    dt = 1 / (bit_rate * samples_per_ui)
    n_rise = int(round(rise_time_s / dt))
    if not 1 <= n_rise <= samples_per_ui or not np.isclose(n_rise * dt, rise_time_s, rtol=1e-9):
        raise ValueError("Rise time must be an integer number of samples within one UI")
    values = np.repeat(bits, samples_per_ui).reshape(len(bits), samples_per_ui)
    previous = np.r_[0.0, bits[:-1]]
    values[:, :n_rise] = previous[:, None] + (bits - previous)[:, None] * np.arange(n_rise)[None, :] / n_rise
    return (source_vpp * values).ravel()


def simulate(bits, h, samples_per_ui=50, rise_time_s=RISE_TIME, source_vpp=SOURCE_VPP, bit_rate=BIT_RATE):
    """Linear overlap-add convolution; preserves every bit, including a partial block."""
    # One following zero bit completes a finite falling ramp after the last bit.
    x = make_nrz(np.r_[bits, 0], samples_per_ui, rise_time_s, source_vpp, bit_rate)
    return oaconvolve(x * MATCHING_DIVIDER, np.asarray(h), mode="full")


def pulse_response(h, samples_per_ui=50, rise_time_s=RISE_TIME, source_vpp=SOURCE_VPP, bit_rate=BIT_RATE):
    x = make_nrz([1, 0], samples_per_ui, rise_time_s, source_vpp, bit_rate)
    return oaconvolve(x * MATCHING_DIVIDER, np.asarray(h), mode="full")


def pda_phase_scan(pulse, samples_per_ui=50):
    """Conservative binary worst-case margin from all symbol-spaced SBR taps."""
    p = np.asarray(pulse)
    largest = int(np.argmax(p))
    candidates = np.arange(max(0, largest - samples_per_ui), min(len(p), largest + samples_per_ui + 1))
    absolute_sums = np.array([np.abs(p[phase::samples_per_ui]).sum() for phase in range(samples_per_ui)])
    margins = 2 * p[candidates] - absolute_sums[candidates % samples_per_ui]
    best = int(np.argmax(margins))
    return int(candidates[best]), candidates, margins


def analyze_eye(y, bits, pulse, samples_per_ui=50, skip_bits=800, bit_rate=BIT_RATE):
    cursor, scan_indices, scan_margins = pda_phase_scan(pulse, samples_per_ui)
    ids = np.arange(skip_bits, len(bits))
    indices = ids * samples_per_ui + cursor
    valid = indices < len(y)
    ids, indices = ids[valid], indices[valid]
    if not len(ids) or len(np.unique(np.asarray(bits)[ids])) != 2:
        raise ValueError("Insufficient steady-state data containing both symbol values")
    levels = np.asarray(y)[indices]
    labels = np.asarray(bits)[ids]
    one_min = float(levels[labels == 1].min())
    zero_max = float(levels[labels == 0].max())
    taps = np.asarray(pulse)[cursor % samples_per_ui::samples_per_ui].copy()
    main = cursor // samples_per_ui
    main_value = float(taps[main])
    taps[main] = 0.0
    high_min = float(main_value + taps[taps < 0].sum())
    low_max = float(taps[taps > 0].sum())
    return {
        "center_eye_height_v": one_min - zero_max,
        "pda_eye_height_v": high_min - low_max,
        "center_time_s": cursor / (bit_rate * samples_per_ui),
        "cursor_index": cursor,
        "phase_samples": cursor % samples_per_ui,
        "one_min_v": one_min,
        "zero_max_v": zero_max,
        "pda_high_min_v": high_min,
        "pda_low_max_v": low_max,
        "main_cursor_v": main_value,
        "isi_absolute_sum_v": float(np.abs(taps).sum()),
        "analyzed_bits": int(len(ids)),
        "skip_bits": int(skip_bits),
        "scan_indices": scan_indices,
        "scan_margins_v": scan_margins,
    }


def _json_eye(result):
    return {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}


def _db(values):
    return 20 * np.log10(np.maximum(np.abs(values), 1e-15))


def _save(fig, output, name):
    fig.savefig(output / name, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _eye_histogram(y, info, bit_count, samples_per_ui, skip_bits, voltage_edges):
    start = skip_bits * samples_per_ui + info["cursor_index"] - samples_per_ui // 2
    stop = bit_count * samples_per_ui + info["cursor_index"] - samples_per_ui // 2
    values = y[start:min(stop, len(y))]
    columns = np.arange(len(values)) % samples_per_ui
    rows = np.searchsorted(voltage_edges, values, side="right") - 1
    valid = (rows >= 0) & (rows < len(voltage_edges) - 1)
    hist = np.zeros((len(voltage_edges) - 1, samples_per_ui), dtype=np.int64)
    np.add.at(hist, (rows[valid], columns[valid]), 1)
    normalized = hist / np.maximum(hist.sum(axis=0, keepdims=True), 1)
    return np.hstack([normalized, normalized])


def _tail_info(full_h, dt, cutoff_ns):
    index = int(round(cutoff_ns * 1e-9 / dt))
    if index <= 0 or index >= len(full_h):
        raise ValueError("Tail cutoff must lie strictly inside the FFT period")
    energy = float(np.dot(full_h, full_h))
    tail = full_h[index:]
    return {"cutoff_ns": cutoff_ns, "samples": index,
            "discarded_energy_fraction": float(np.dot(tail, tail) / energy),
            "discarded_l1": float(np.abs(tail).sum()),
            "retained_dc_sum": float(full_h[:index].sum()),
            "full_dc_sum": float(full_h.sum())}


def run(args):
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = Path(args.input)
    network = read_touchstone(source)
    f, channel = network.frequency_hz, network.s[:, 1, 0]
    if not np.isclose(network.reference_ohm, 50):
        raise ValueError("This assignment requires a 50 ohm reference")
    settings = dict(fp1_ghz=args.fp1_ghz, fp2_ghz=args.fp2_ghz, adc_db=args.adc_db, aac_db=args.aac_db)
    peak = ctle_peak(**settings)
    if not (1 <= args.fp1_ghz <= 40 and 1 <= args.fp2_ghz <= 40 and -20 <= args.adc_db <= 6):
        raise ValueError("CTLE poles or DC gain violate assignment bounds")
    if peak["peak_gain_db"] > 6 + 1e-10:
        raise ValueError("The analytic full-frequency CTLE maximum exceeds 6 dB")
    ctle = ctle_response(f, **settings)
    time, h_before_full = discrete_impulse(f, channel, args.samples_per_ui, args.alpha)
    _, h_after_full = discrete_impulse(f, channel * ctle, args.samples_per_ui, args.alpha)
    dt = 1 / (BIT_RATE * args.samples_per_ui)
    cutoff = int(round(args.tail_ns * 1e-9 / dt))
    if not 0 < cutoff < len(time):
        raise ValueError("Tail cutoff must lie inside the FFT period")
    h_before, h_after = h_before_full[:cutoff], h_after_full[:cutoff]
    skip_bits = int(np.ceil(args.tail_ns * 1e-9 * BIT_RATE))
    if args.bits < skip_bits + 100:
        raise ValueError(f"Need at least {skip_bits + 100} bits for steady-state analysis")
    bits = np.random.default_rng(args.seed).integers(0, 2, size=args.bits, dtype=np.uint8)
    y_before = simulate(bits, h_before, args.samples_per_ui)
    y_after = simulate(bits, h_after, args.samples_per_ui)
    p_before = pulse_response(h_before, args.samples_per_ui)
    p_after = pulse_response(h_after, args.samples_per_ui)
    before = analyze_eye(y_before, bits, p_before, args.samples_per_ui, skip_bits)
    after = analyze_eye(y_after, bits, p_after, args.samples_per_ui, skip_bits)

    tail_convergence = []
    for duration in sorted(set([args.tail_ns, 40.0, 80.0])):
        if duration * 1e-9 >= time[-1]:
            continue
        entry = {"before": _tail_info(h_before_full, dt, duration), "after": _tail_info(h_after_full, dt, duration)}
        count = int(round(duration * 1e-9 / dt))
        for key, full in [("before", h_before_full), ("after", h_after_full)]:
            pulse = pulse_response(full[:count], args.samples_per_ui)
            center, _, margins = pda_phase_scan(pulse, args.samples_per_ui)
            entry[key].update(pda_eye_height_v=float(margins.max()), cursor_index=center)
        tail_convergence.append(entry)

    # Keep the complete measured-frequency FFT period as a conservative check:
    # this also includes very small periodic/precausal ringing near its end.
    full_period_pda = {}
    for key, full in [("before", h_before_full), ("after", h_after_full)]:
        pulse = pulse_response(full, args.samples_per_ui)
        center, _, margins = pda_phase_scan(pulse, args.samples_per_ui)
        full_period_pda[key] = {"pda_eye_height_v": float(margins.max()), "cursor_index": center,
                                "center_time_s": center * dt}
    fine_sps = 2 * args.samples_per_ui
    fine_time, fine_h = discrete_impulse(f, channel * ctle, fine_sps, args.alpha)
    fine_cutoff = int(round(args.tail_ns * 1e-9 * BIT_RATE * fine_sps))
    fine_pulse = pulse_response(fine_h[:fine_cutoff], fine_sps)
    fine_center, _, fine_margins = pda_phase_scan(fine_pulse, fine_sps)
    fine_full_pulse = pulse_response(fine_h, fine_sps)
    fine_full_center, _, fine_full_margins = pda_phase_scan(fine_full_pulse, fine_sps)
    time_step_convergence = {
        "coarse_dt_s": dt, "fine_dt_s": 1 / (BIT_RATE * fine_sps),
        "coarse_samples_per_ui": args.samples_per_ui, "fine_samples_per_ui": fine_sps,
        "coarse_cutoff_pda_v": after["pda_eye_height_v"],
        "fine_cutoff_pda_v": float(fine_margins.max()),
        "cutoff_pda_difference_v": float(fine_margins.max() - after["pda_eye_height_v"]),
        "coarse_full_period_pda_v": full_period_pda["after"]["pda_eye_height_v"],
        "fine_full_period_pda_v": float(fine_full_margins.max()),
        "full_period_pda_difference_v": float(fine_full_margins.max() - full_period_pda["after"]["pda_eye_height_v"]),
        "fine_cutoff_center_time_s": fine_center / (BIT_RATE * fine_sps),
        "fine_full_period_center_time_s": fine_full_center / (BIT_RATE * fine_sps),
    }

    metrics = {
        "source": {"filename": source.name, "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                   "points": len(f), "df_hz": float(f[1] - f[0]), "fmax_hz": float(f[-1]),
                   "reference_ohm": network.reference_ohm},
        "settings": {"bit_rate_bps": BIT_RATE, "rise_time_s": RISE_TIME, "source_vpp": SOURCE_VPP,
                     "matching_divider": MATCHING_DIVIDER, "samples_per_ui": args.samples_per_ui,
                     "dt_s": dt, "bits": args.bits, "seed": args.seed, "alpha": args.alpha,
                     "random_bit_generator": "numpy.default_rng / PCG64, uint8 output",
                     "tail_ns": args.tail_ns, "skip_bits": skip_bits, "fft_period_s": float(len(time) * dt),
                     "window_fmax_hz": float(f[-1]), "convolution": "scipy.signal.oaconvolve, full linear convolution"},
        "ctle": {**settings, **peak},
        "before": _json_eye(before), "after": _json_eye(after),
        "tail_convergence": tail_convergence,
        "full_period_pda": full_period_pda,
        "time_step_convergence": time_step_convergence,
        "requirements": {"actual_peak_at_most_6db": bool(peak["peak_gain_db"] <= 6 + 1e-10),
                         "random_center_eh_at_least_300mv": bool(after["center_eye_height_v"] >= .300),
                         "pda_center_eh_at_least_300mv": bool(after["pda_eye_height_v"] >= .300),
                         "full_period_pda_at_least_300mv": bool(full_period_pda["after"]["pda_eye_height_v"] >= .300)},
        "report_checks": [
            {"label": "Analytic gain bound", "detail": f"Global CTLE maximum {peak['peak_gain_db']:.9f} dB, including frequencies beyond the measured band."},
            {"label": "Impulse normalization", "detail": f"Full channel kernel sum {h_before_full.sum():.12g}; measured channel DC gain {channel[0].real:.12g}. The source divider is applied separately once."},
            {"label": "All bits retained", "detail": f"All {args.bits} generated bits are convolved; the first {skip_bits} bits are excluded only from measurement."},
            {"label": "Finite numerical results", "detail": f"Both complete output waveforms are finite: {bool(np.all(np.isfinite(y_before)) and np.all(np.isfinite(y_after)))}."},
            {"label": "Full-period conservative margin", "detail": f"Full {len(time) * dt * 1e9:g} ns kernel PDA EH {full_period_pda['after']['pda_eye_height_v'] * 1e3:.6f} mV; no impulse renormalization."},
            {"label": "Time step convergence", "detail": f"Halving dt changes cutoff PDA EH by {time_step_convergence['cutoff_pda_difference_v'] * 1e3:.6f} mV and full-period PDA EH by {time_step_convergence['full_period_pda_difference_v'] * 1e3:.6f} mV."},
        ],
        "conventions": {"eye_height": "min(samples conditioned on transmitted 1) - max(samples conditioned on transmitted 0)",
                        "sampling_center": "maximum all-tap PDA eye height around the main cursor; independently aligned for each path",
                        "pda": "main cursor minus sum of absolute values of every other symbol-spaced SBR tap",
                        "random_bits": "all generated bits simulated; initial tail-length bits excluded from measurements",
                        "voltage": "unipolar source 0..2 V; a single 0.5 matched-source divider before convolution",
                        "window": "Kaiser positive-frequency taper using the measured 80 GHz endpoint; bins above measurement are zero",
                        "tail": "FFT-periodic impulse truncated without renormalization; 40/80 ns PDA convergence also reported"},
    }
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .25})
    fig, ax = plt.subplots(figsize=(9, 4.7))
    ax.plot(f / 1e9, _db(channel), label="Channel S21")
    ax.plot(f / 1e9, _db(ctle), label="CTLE")
    ax.plot(f / 1e9, _db(channel * ctle * MATCHING_DIVIDER), label="Overall voltage transfer: S21 x CTLE / 2")
    ax.axhline(6, color="black", ls="--", lw=.8, label="CTLE peak limit")
    ax.scatter([peak["peak_frequency_hz"] / 1e9], [peak["peak_gain_db"]], s=22, color="C1")
    ax.set(xlabel="Frequency (GHz)", ylabel="Magnitude (dB)", xlim=(0, f[-1] / 1e9), ylim=(-50, 10))
    ax.legend(ncol=2)
    _save(fig, output, "frequency_response.png")

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    for h, label in [(h_before_full, "Channel"), (h_after_full, "Channel x CTLE")]:
        axes[0].plot(time * 1e9, h, label=label, lw=.9)
        axes[1].plot(time * 1e9, _db(h / np.max(np.abs(h))), label=label, lw=.9)
    axes[0].set(ylabel="Discrete impulse coefficient")
    axes[0].legend()
    axes[1].set(xlabel="Time (ns)", ylabel="Relative amplitude (dB)", ylim=(-120, 5), xlim=(0, min(40, time[-1] * 1e9)))
    for ax in axes:
        ax.axvline(args.tail_ns, ls="--", color="black", lw=.8)
    _save(fig, output, "impulse_response.png")

    vmin = min(y_before[skip_bits * args.samples_per_ui:].min(), y_after[skip_bits * args.samples_per_ui:].min())
    vmax = max(y_before[skip_bits * args.samples_per_ui:].max(), y_after[skip_bits * args.samples_per_ui:].max())
    edges = np.linspace(vmin - .03, vmax + .03, 601)
    hist_before = _eye_histogram(y_before, before, args.bits, args.samples_per_ui, skip_bits, edges)
    hist_after = _eye_histogram(y_after, after, args.bits, args.samples_per_ui, skip_bits, edges)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True, constrained_layout=True)
    for ax, hist, info, title in zip(axes, [hist_before, hist_after], [before, after], ["Post-channel", "Post-CTLE"]):
        im = ax.imshow(np.log10(np.maximum(hist, 1e-6)), origin="lower", aspect="auto", extent=[0, 2, edges[0], edges[-1]], vmin=-6, vmax=-1, cmap="viridis")
        ax.set(title=f"{title}: center EH = {info['center_eye_height_v'] * 1e3:.1f} mV", xlabel="Aligned time (UI)", xlim=(0, 2))
        ax.axvline(.5, color="white", ls="--", lw=.7)
        ax.axvline(1.5, color="white", ls="--", lw=.7)
    axes[0].set_ylabel("Voltage (V)")
    fig.colorbar(im, ax=axes, label="log10 probability per voltage bin")
    _save(fig, output, "eye_diagrams.png")

    fig, ax = plt.subplots(figsize=(9, 4.3))
    for info, label in [(before, "Post-channel"), (after, "Post-CTLE")]:
        phase = (info["scan_indices"] - info["cursor_index"]) / args.samples_per_ui
        ax.plot(phase, 1e3 * info["scan_margins_v"], label=label)
    ax.axhline(300, ls="--", color="black", lw=.8, label="300 mV requirement")
    ax.axvline(0, ls=":", color="black", lw=.8)
    ax.set(xlabel="Sampling offset from selected center (UI)", ylabel="All-tap PDA eye height (mV)", xlim=(-.5, .5))
    ax.legend()
    _save(fig, output, "pda_eye_height.png")

    waveform_data = {}
    for start_ns in [100, 200, 300]:
        start = int(round(start_ns * 1e-9 / dt))
        count = int(round(1e-9 / dt))
        if start + count >= min(len(y_before), len(y_after)):
            continue
        indices = np.arange(start, start + count + 1)
        fig, ax = plt.subplots(figsize=(10, 3.6))
        ax.plot(indices * dt * 1e9, y_before[indices], label="Post-channel", lw=1)
        ax.plot(indices * dt * 1e9, y_after[indices], label="Post-CTLE", lw=1)
        ax.set(xlabel="Time (ns)", ylabel="Voltage (V)", xlim=(start_ns, start_ns + 1))
        ax.legend()
        _save(fig, output, f"waveform_{start_ns}ns.png")
        waveform_data[f"waveform_{start_ns}ns_before"] = y_before[indices]
        waveform_data[f"waveform_{start_ns}ns_after"] = y_after[indices]
    np.savez_compressed(output / "analysis_data.npz", frequency_hz=f, channel=channel, ctle=ctle,
                        impulse_time_s=time[:cutoff], h_before=h_before, h_after=h_after,
                        pulse_before=p_before, pulse_after=p_after, bits=bits,
                        histogram_voltage_edges=edges, histogram_before=hist_before, histogram_after=hist_after,
                        **waveform_data)
    print(json.dumps({"ctle": metrics["ctle"], "before": metrics["before"], "after": metrics["after"], "requirements": metrics["requirements"]}, indent=2))
    return metrics


def parser():
    base = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=base / "instructions" / "HW7_stripline_50cm.s2p")
    p.add_argument("--output", type=Path, default=base / "results")
    p.add_argument("--bits", type=int, default=2 ** 16)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--samples-per-ui", type=int, default=50)
    p.add_argument("--alpha", type=float, default=1.0)
    p.add_argument("--tail-ns", type=float, default=20.0)
    p.add_argument("--fp1-ghz", type=float, default=20.6)
    p.add_argument("--fp2-ghz", type=float, default=22.8)
    p.add_argument("--adc-db", type=float, default=-4.91)
    p.add_argument("--aac-db", type=float, default=11.48)
    return p


if __name__ == "__main__":
    run(parser().parse_args())
