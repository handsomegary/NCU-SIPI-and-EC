"""Reproducible HW8 transmitter FFE design and time-domain simulation.

The source waveform is constructed as discrete FFE levels followed by finite
linear ramps. The channel is applied once using sectional overlap-add. The
analytic single-bit response is a separate reference, not a shortcut used to
synthesize the random-data output. All voltages refer to a 0..2 V source and
the matched-source transfer H = S21 / 2.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.special import i0


CASES = (
    ("case1_no_ffe", "No FFE", 0, 0),
    ("case2_pre1", "1 precursor", 1, 0),
    ("case3_post1", "1 postcursor", 0, 1),
    ("case4_pre1_post1", "1 precursor + 1 postcursor", 1, 1),
    ("case5_pre1_post2", "1 precursor + 2 postcursors", 1, 2),
)


def read_touchstone(path):
    """Read a two-port Touchstone 1.x file (RI, MA or DB), preserving S21 order.

    Returns frequencies in Hz, S with shape (frequency, output, input), and
    real reference impedance. Continuation lines and inline comments work.
    The transform routines additionally require uniform samples including DC.
    """
    path = Path(path)
    options = ["ghz", "s", "ma", "r", "50"]
    numbers = []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            supplied = line[1:].lower().split()
            if supplied:
                options = supplied
            continue
        if line.startswith("["):
            raise ValueError("This reader supports two-port Touchstone 1.x only")
        numbers.extend(float(v.replace("D", "E").replace("d", "e")) for v in line.split())
    if not numbers or len(numbers) % 9:
        raise ValueError("Each two-port Touchstone record must contain 9 numbers")
    unit = next((x for x in options if x in {"hz", "khz", "mhz", "ghz"}), "ghz")
    fmt = next((x for x in options if x in {"ri", "ma", "db"}), "ma")
    if "s" not in options:
        raise ValueError("Only S parameters are supported")
    z0 = float(options[options.index("r") + 1]) if "r" in options else 50.0
    rows = np.asarray(numbers, float).reshape(-1, 9)
    f = rows[:, 0] * {"hz": 1, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}[unit]
    if np.any(~np.isfinite(rows)) or np.any(np.diff(f) <= 0) or f[0] < 0:
        raise ValueError("Frequencies must be finite, nonnegative and increasing")
    pairs = rows[:, 1:].reshape(-1, 4, 2)
    if fmt == "ri":
        values = pairs[:, :, 0] + 1j * pairs[:, :, 1]
    else:
        magnitude = 10 ** (pairs[:, :, 0] / 20) if fmt == "db" else pairs[:, :, 0]
        values = magnitude * np.exp(1j * np.deg2rad(pairs[:, :, 1]))
    s = np.empty((len(f), 2, 2), complex)
    s[:, 0, 0], s[:, 1, 0] = values[:, 0], values[:, 1]
    s[:, 0, 1], s[:, 1, 1] = values[:, 2], values[:, 3]
    return f, s, z0


def kaiser_window(f, alpha=0.5):
    """One-sided lecture Kaiser window; endpoints are not halved for irfft."""
    f = np.asarray(f, float)
    if f.ndim != 1 or len(f) < 2 or f[-1] <= 0 or alpha < 0:
        raise ValueError("A positive band limit and nonnegative alpha are required")
    return i0(np.pi * alpha * np.sqrt(np.maximum(0.0, 1 - (f / f[-1]) ** 2))) / i0(np.pi * alpha)


def _transform_spectrum(f, H, dt, alpha):
    f, H = np.asarray(f, float), np.asarray(H, complex)
    if len(f) < 2 or f.shape != H.shape or dt <= 0:
        raise ValueError("Frequency/transfer vectors and dt are invalid")
    df = f[1] - f[0]
    if f[0] != 0 or df <= 0 or not np.allclose(np.diff(f), df, rtol=1e-9, atol=1e-5):
        raise ValueError("Uniform increasing frequency samples including DC are required")
    n = round(1 / (df * dt))
    if n < 2 or not np.isclose(n * df * dt, 1, rtol=1e-10, atol=1e-12):
        raise ValueError("dt and df must define an integral inverse-transform length")
    if f[-1] > 1 / (2 * dt) * (1 + 1e-12) or len(f) > n // 2 + 1:
        raise ValueError("Measured band exceeds the output Nyquist frequency")
    spectrum = np.zeros(n // 2 + 1, complex)
    spectrum[:len(f)] = H * kaiser_window(f, alpha)
    # Hermitian symmetry requires real self-conjugate DC/Nyquist endpoints.
    spectrum[0] = spectrum[0].real
    if n % 2 == 0:
        spectrum[-1] = spectrum[-1].real
    return n, spectrum


def frequency_impulse(f, H, dt=0.5e-12, alpha=0.5):
    """Return the discrete channel FIR over one 1/df period.

    No dt factor belongs in discrete convolution with this FIR. A continuous
    impulse response would be h/dt. DC and Nyquist are counted exactly once.
    Frequencies above the measured band are zero, with no extrapolation.
    """
    n, spectrum = _transform_spectrum(f, H, dt, alpha)
    return np.arange(n) * dt, irfft(spectrum, n=n)


def get_sbr(f, H, tb=25e-12, tr=5e-12, dt=0.5e-12, vpp=1.0, alpha=0.5):
    """Analytic causal trapezoidal-pulse response, sampled every dt seconds.

    Pulse: ramp 0..Tr, plateau Tr..Tb, ramp Tb..Tb+Tr. Its spectrum is
    Vpp*Tb*sinc(f*Tb)*sinc(f*Tr)*exp[-j*pi*f*(Tb+Tr)].
    """
    if not 0 <= tr <= tb or tb <= 0 or vpp <= 0:
        raise ValueError("Require 0 <= Tr <= Tb and positive Tb and Vpp")
    n, spectrum = _transform_spectrum(f, H, dt, alpha)
    fgrid = np.arange(len(spectrum)) / (n * dt)
    pulse = vpp * tb * np.sinc(fgrid * tb) * np.sinc(fgrid * tr) * np.exp(-1j * np.pi * fgrid * (tb + tr))
    spectrum *= pulse
    spectrum[0] = spectrum[0].real
    if n % 2 == 0:
        spectrum[-1] = spectrum[-1].real
    return np.arange(n) * dt, irfft(spectrum, n=n) / dt


def trapezoid_pulse(tb, tr, dt, vpp=1.0):
    """Sample the same causal trapezoid used by get_sbr, including its zero end."""
    if not 0 <= tr <= tb or tb <= 0 or dt <= 0:
        raise ValueError("Require 0 <= Tr <= Tb and positive Tb and dt")
    n_ui, n_rise = round(tb / dt), round(tr / dt)
    if not np.isclose(tb / dt, n_ui, rtol=1e-10, atol=1e-10) or not np.isclose(tr / dt, n_rise, rtol=1e-10, atol=1e-10):
        raise ValueError("Pulse durations must be integral multiples of dt")
    if n_rise == 0:
        return np.r_[np.full(n_ui, vpp), 0.0]
    t = np.arange(n_ui + n_rise + 1) * dt
    return vpp * np.clip(np.minimum(t / tr, (tb + tr - t) / tr), 0, 1)


def get_cursor(sbr, n_per_ui, npre, npost):
    """Sample cursors around the largest positive single-bit response sample."""
    if n_per_ui < 1 or npre < 0 or npost < 0:
        raise ValueError("Invalid cursor counts")
    sbr = np.asarray(sbr, float)
    indices = int(np.argmax(sbr)) + np.arange(-npre, npost + 1) * n_per_ui
    result = np.zeros(len(indices))
    valid = (indices >= 0) & (indices < len(sbr))
    result[valid] = sbr[indices[valid]]
    return result


def design_zf(a, npre):
    """Lecture's truncated-cursor ZF, followed by sum(abs(taps))=1.

    The finite model sets every channel cursor outside a's range to zero.
    Therefore a small matrix residual does not guarantee zero physical ISI.
    """
    a = np.asarray(a, float)
    if a.ndim != 1 or not len(a) or not 0 <= npre < len(a):
        raise ValueError("Invalid cursor array or main cursor index")
    n = len(a)
    ix = np.arange(n)[:, None] - np.arange(n)[None, :] + npre
    A = np.zeros((n, n))
    valid = (ix >= 0) & (ix < n)
    A[valid] = a[ix[valid]]
    target = np.zeros(n)
    target[npre] = 1.0
    raw = np.linalg.solve(A, target)
    norm = float(np.sum(np.abs(raw)))
    taps = raw / norm
    response = A @ taps
    scaled_target = target / norm
    return {"taps": taps, "matrix": A, "target": scaled_target,
            "response": response, "norm": norm,
            "residual": response - scaled_target}


def zf_taps(a, npre):
    return design_zf(a, npre)["taps"]


def ffe_levels(bits, taps, vpp=2.0):
    """Apply causal taps to source levels before adding finite transition time."""
    bits, taps = np.asarray(bits, float), np.asarray(taps, float)
    if bits.ndim != 1 or taps.ndim != 1 or not len(bits) or not len(taps):
        raise ValueError("Bits and taps must be nonempty one-dimensional arrays")
    norm = np.sum(np.abs(taps))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("FFE taps must have a finite positive absolute sum")
    return np.convolve(bits * vpp, taps / norm, mode="full")


def levels_to_waveform(levels, n_per_ui, n_rise, previous=0.0):
    """All jumps, regardless of their magnitude, take exactly n_rise samples."""
    levels = np.asarray(levels, float)
    if levels.ndim != 1 or not 0 <= n_rise <= n_per_ui or n_per_ui < 1:
        raise ValueError("Require 0 <= n_rise <= n_per_ui and positive n_per_ui")
    x = np.repeat(levels, n_per_ui).reshape(-1, n_per_ui)
    if n_rise == 0:
        return x.ravel()
    previous_levels = np.r_[previous, levels[:-1]]
    ramp = np.arange(n_rise) / n_rise
    x[:, :n_rise] = previous_levels[:, None] + (levels - previous_levels)[:, None] * ramp
    return x.ravel()


def sectional_convolve(x, h, section_samples=16384):
    """True linear sectional overlap-add with an independently transformed FIR."""
    x, h = np.asarray(x, float), np.asarray(h, float)
    if x.ndim != 1 or h.ndim != 1 or not len(x) or not len(h) or section_samples < 1:
        raise ValueError("Nonempty vectors and a positive section size are required")
    nfft = next_fast_len(int(section_samples) + len(h) - 1)
    H = rfft(h, n=nfft)
    result = np.zeros(len(x) + len(h) - 1)
    for start in range(0, len(x), int(section_samples)):
        block = x[start:start + int(section_samples)]
        valid = len(block) + len(h) - 1
        result[start:start + valid] += irfft(rfft(block, n=nfft) * H, n=nfft)[:valid]
    return result


def eye_histogram(y, n_per_ui, voltage_bins=400, voltage_range=None):
    """Return a one-UI conditional voltage histogram, increasing voltage upward."""
    y = np.asarray(y, float)
    if len(y) < n_per_ui or n_per_ui < 1 or voltage_bins < 2:
        raise ValueError("Not enough eye samples or invalid bin counts")
    if voltage_range is None:
        margin = max(float(np.ptp(y)) * 0.04, 0.002)
        lo, hi = float(np.min(y) - margin), float(np.max(y) + margin)
    else:
        lo, hi = map(float, voltage_range)
    if not hi > lo:
        raise ValueError("Voltage limits must increase")
    v = np.floor((y - lo) / (hi - lo) * voltage_bins).astype(int)
    v = np.clip(v, 0, voltage_bins - 1)
    phase = np.arange(len(y)) % n_per_ui
    counts = np.bincount(v * n_per_ui + phase, minlength=voltage_bins * n_per_ui).reshape(voltage_bins, n_per_ui)
    totals = counts.sum(axis=0)
    probability = np.divide(counts, totals, out=np.zeros_like(counts, float), where=totals != 0)
    return probability, [0.0, 1.0, lo, hi]


def _equalized_sbr(sbr, taps, n_per_ui):
    result = np.zeros(len(sbr) + (len(taps) - 1) * n_per_ui)
    for j, tap in enumerate(taps):
        result[j * n_per_ui:j * n_per_ui + len(sbr)] += tap * sbr
    return result


def _physical_cursors(response, index, n_per_ui):
    offsets = np.arange(-index // n_per_ui - 1, (len(response) - index) // n_per_ui + 1)
    indices = index + offsets * n_per_ui
    keep = (indices >= 0) & (indices < len(response))
    return offsets[keep], response[indices[keep]]


def _jsonable(value):
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 180})
    return plt


def run_simulation(data, output, bits=65536, seed=0):
    """Run all five cases and write anonymous, report-ready figures and metrics."""
    if bits < 64:
        raise ValueError("Use at least 64 payload bits")
    data, output = Path(data), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt = _plotting()
    import platform
    import scipy
    import matplotlib
    tb, tr, dt, vpp, alpha, t2 = 25e-12, 5e-12, 0.5e-12, 2.0, 0.5, 20e-9
    n_per_ui, n_rise = round(tb / dt), round(tr / dt)
    f, s, z0 = read_touchstone(data)
    H = s[:, 1, 0] / 2  # The sole source-divider factor in the simulation.
    t, h_period = frequency_impulse(f, H, dt, alpha)
    fir_count = min(round(t2 / dt) + 1, len(h_period))
    h = h_period[:fir_count].copy()
    tsbr, sbr = get_sbr(f, H, tb, tr, dt, vpp, alpha)
    kmain = int(np.argmax(sbr))
    # Compare the analytic continuous pulse with the directly sampled ramps.
    pulse = trapezoid_pulse(tb, tr, dt, vpp)
    pulse_filtered = sectional_convolve(pulse, h_period)
    comparison_length = min(len(sbr), len(pulse_filtered))
    pulse_error = float(np.max(np.abs(pulse_filtered[:comparison_length] - sbr[:comparison_length])))
    pulse_rms_error = float(np.sqrt(np.mean((pulse_filtered[:comparison_length] - sbr[:comparison_length]) ** 2)))
    # Exact sampled-wave model: its periodic DFT and the same channel must agree.
    n, Hgrid = _transform_spectrum(f, H, dt, alpha)
    sampled_sbr = irfft(rfft(pulse, n=n) * Hgrid, n=n)
    sampled_boundary_error = float(np.max(np.abs(sampled_sbr[:fir_count] - pulse_filtered[:fir_count])))
    folded_linear = pulse_filtered[:n].copy()
    folded_linear[:len(pulse_filtered) - n] += pulse_filtered[n:]
    sampled_reference_error = float(np.max(np.abs(sampled_sbr - folded_linear)))
    # Refine dt while retaining the measured band and compare sampled-pulse
    # errors over the same physical main-response window.
    fine_dt = dt / 2
    _, fine_h = frequency_impulse(f, H, fine_dt, alpha)
    _, fine_sbr = get_sbr(f, H, tb, tr, fine_dt, vpp, alpha)
    fine_pulse = trapezoid_pulse(tb, tr, fine_dt, vpp)
    fine_filtered = sectional_convolve(fine_pulse, fine_h)
    response_window = slice(max(0, kmain - 8 * n_per_ui), kmain + 24 * n_per_ui)
    coarse_window_error = float(np.max(np.abs(pulse_filtered[:n][response_window] - sbr[response_window])))
    fine_window = slice(response_window.start * 2, response_window.stop * 2)
    fine_window_error = float(np.max(np.abs(fine_filtered[:len(fine_sbr)][fine_window] - fine_sbr[fine_window])))
    ola_probe = np.random.default_rng(17).normal(size=257)
    ola_reference = np.convolve(ola_probe, h[:251])
    ola_error = float(np.max(np.abs(sectional_convolve(ola_probe, h[:251], 61) - ola_reference)))
    guard_bits = int(np.ceil((len(h) + kmain) / n_per_ui)) + 16
    rng = np.random.default_rng(seed)
    common_bits = rng.integers(0, 2, bits + 2 * guard_bits, dtype=np.int8)
    source_indices = guard_bits + np.arange(bits)
    source_payload = common_bits[source_indices]
    summary = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "parameters": {"bit_rate_bps": 1 / tb, "tb_s": tb, "tr_s": tr,
                       "dt_s": dt, "samples_per_ui": n_per_ui, "rise_samples": n_rise,
                       "source_vpp_v": vpp, "source_levels_v": [0, vpp],
                       "kaiser_alpha": alpha, "fir_truncation_s": t2,
                       "payload_bits": bits, "guard_bits_each_end": guard_bits, "seed": seed,
                       "section_samples": 16384, "input_file": data.name,
                       "source_sha256": hashlib.sha256(data.read_bytes()).hexdigest()},
        "channel": {"reference_impedance_ohm": z0, "frequency_samples": len(f),
                    "frequency_step_hz": f[1] - f[0], "maximum_frequency_hz": f[-1],
                    "dft_samples": len(h_period), "period_s": len(h_period) * dt,
                    "transfer_convention": "H = S21 / 2, applied exactly once",
                    "dc_gain": float(H[0].real), "finite_fir_dc_gain": float(h.sum()),
                    "omitted_fir_absolute_sum": float(np.sum(np.abs(h_period[fir_count:]))),
                    "dc_gain_relative_truncation_error": float(abs(h.sum() - H[0].real) / abs(H[0].real)),
                    "omitted_fir_energy_fraction": float(np.sum(h_period[fir_count:] ** 2) / np.sum(h_period ** 2)),
                    "main_sbr_index": kmain, "main_sbr_time_s": kmain * dt},
        "verification": {"ola_vs_direct_max_abs_v": ola_error,
                         "continuous_sbr_vs_sampled_pulse_max_abs_v": pulse_error,
                         "continuous_sbr_vs_sampled_pulse_rms_v": pulse_rms_error,
                         "continuous_sbr_vs_sampled_pulse_relative_peak": pulse_error / float(np.max(sbr)),
                         "sampled_pulse_fft_vs_folded_linear_filter_max_abs_v": sampled_reference_error,
                         "sampled_pulse_periodic_vs_linear_boundary_max_abs_v_before_t2": sampled_boundary_error,
                         "dt_convergence_main_window_coarse_max_abs_v": coarse_window_error,
                         "dt_convergence_main_window_fine_max_abs_v": fine_window_error,
                         "dt_convergence_fine_dt_s": fine_dt,
                         "pulse_comparison_note": "Continuous-pulse SBR and sample-first ramp filtering differ by finite time sampling. Folding the full linear-filter result at the DFT period verifies exact agreement with the discrete periodic reference; the unfurled linear comparison separately exposes periodic boundary contributions.",
                         "metric_note": "Empirical sampled vertical gaps are extrema over this finite noiseless bit sequence, not BER or noise-margin estimates."},
        "cases": [],
    }
    eye_waves, sbr_cache = [], []

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), constrained_layout=True)
    axes[0].plot(f / 1e9, 20 * np.log10(np.maximum(np.abs(H), 1e-15)), label="S21 / 2")
    axes[0].plot(f / 1e9, 20 * np.log10(np.maximum(np.abs(H * kaiser_window(f, alpha)), 1e-15)), label="with Kaiser")
    axes[0].set(xlabel="Frequency (GHz)", ylabel="Magnitude (dB)", title="50 cm stripline transfer")
    axes[0].legend()
    axes[1].plot(t[:fir_count] * 1e9, h, lw=1)
    axes[1].set(xlabel="Time (ns)", ylabel="Discrete impulse coefficient", title="Channel FIR, first 20 ns")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.savefig(output / "channel_response.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.8), constrained_layout=True)
    ax.plot(tsbr * 1e9, sbr, label="Analytic causal-pulse SBR")
    ax.plot(np.arange(len(pulse_filtered)) * dt * 1e9, pulse_filtered, "--", lw=1, label="Sampled pulse filtered by channel")
    cursor_indices = kmain + np.arange(-3, 9) * n_per_ui
    ax.plot(cursor_indices * dt * 1e9, sbr[cursor_indices], "o", ms=3, label="UI-spaced cursors")
    ax.set(xlim=((kmain * dt - 4 * tb) * 1e9, (kmain * dt + 12 * tb) * 1e9), xlabel="Time (ns)", ylabel="Response (V)", title="Single-bit response for a 0 to 2 V input pulse")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(output / "single_bit_response.png")
    plt.close(fig)

    for key, label, npre, npost in CASES:
        a = get_cursor(sbr, n_per_ui, npre, npost)
        zf = design_zf(a, npre)
        taps = zf["taps"]
        eq = _equalized_sbr(sbr, taps, n_per_ui)
        # Search only around the intended main cursor to keep symbol association.
        design_index = kmain + npre * n_per_ui
        left, right = max(0, design_index - n_per_ui // 2), min(len(eq), design_index + n_per_ui // 2 + 1)
        sample_index = left + int(np.argmax(eq[left:right]))
        offsets, physical = _physical_cursors(eq, sample_index, n_per_ui)
        residual = float(np.sum(np.abs(physical[offsets != 0])))
        selected_offsets = np.arange(-npre, npost + 1)
        selected_physical = eq[design_index + selected_offsets * n_per_ui]
        selected_target = np.zeros(len(a))
        selected_target[npre] = eq[design_index]
        levels = ffe_levels(common_bits, taps, vpp)
        x = levels_to_waveform(np.r_[levels, 0.0], n_per_ui, n_rise)
        y = sectional_convolve(x, h)
        sample_positions = source_indices * n_per_ui + sample_index
        decisions = y[sample_positions]
        high_min, low_max = float(decisions[source_payload == 1].min()), float(decisions[source_payload == 0].max())
        eye_start = int(sample_positions[0] - n_per_ui // 2)
        eye_wave = y[eye_start:eye_start + bits * n_per_ui]
        eye_waves.append((key, label, eye_wave.copy()))
        sbr_cache.append((label, eq, sample_index))
        tx_payload = x[guard_bits * n_per_ui:(guard_bits + bits) * n_per_ui]
        metric = {
            "name": key, "label": label, "npre": npre, "npost": npost,
            "taps": taps, "tap_l1_norm": np.sum(np.abs(taps)), "design_cursors_v": a,
            "zf_matrix_v": zf["matrix"], "zf_selected_response_v": zf["response"],
            "zf_matrix_residual_v": np.max(np.abs(zf["residual"])),
            "zf_condition_number": np.linalg.cond(zf["matrix"]),
            "zf_predicted_main_cursor_v": zf["response"][npre],
            "main_cursor_v": eq[sample_index], "design_phase_main_cursor_v": eq[design_index],
            "physical_selected_cursors_at_design_phase_v": selected_physical,
            "physical_selected_residual_v": np.max(np.abs(selected_physical - selected_target)),
            "physical_vs_zf_selected_max_abs_v": np.max(np.abs(selected_physical - zf["response"])),
            "full_residual_isi_l1_v": residual,
            "analytic_cursor_gap_bound_v": float(eq[sample_index] - residual),
            "analytic_cursor_bound_note": "Main minus the L1 sum of other UI-spaced cursors in the analytic finite-period response; not an infinite physical-channel guarantee.",
            "empirical_sampled_vertical_gap_v": high_min - low_max,
            "empirical_high_min_v": high_min, "empirical_low_max_v": low_max,
            "empirical_high_mean_v": float(np.mean(decisions[source_payload == 1])),
            "empirical_low_mean_v": float(np.mean(decisions[source_payload == 0])),
            "tx_min_v": np.min(tx_payload), "tx_max_v": np.max(tx_payload),
            "tx_swing_v": np.ptp(tx_payload), "sample_index": sample_index,
            "sample_latency_s": sample_index * dt, "main_tap_delay_ui": npre,
            "sample_phase_ps": sample_index * dt * 1e12,
            "phase_shift_from_design_ps": (sample_index - design_index) * dt * 1e12,
            "eye_file": f"eye_{key}.png",
        }
        summary["cases"].append(metric)
        if npre == npost == 1:
            metric["snapshots"] = []
            window_bits = min(20, bits // 8)
            for snapshot, fraction in enumerate((0.25, 0.70), 1):
                begin_bit = guard_bits + int(bits * fraction)
                for signal_name, wave, delay in (("tx", x, 0), ("rx", y, sample_index - n_per_ui // 2)):
                    begin = begin_bit * n_per_ui + delay
                    values = wave[begin:begin + window_bits * n_per_ui + 1]
                    time_ps = np.arange(len(values)) * dt * 1e12
                    fig, ax = plt.subplots(figsize=(9, 3.2), constrained_layout=True)
                    ax.plot(time_ps, values, lw=1)
                    ax.set(xlabel="Time from snapshot start (ps)", ylabel="Voltage (V)",
                           title=f"Case 4 {signal_name.upper()} snapshot {snapshot}: {window_bits} UI")
                    ax.grid(alpha=0.25)
                    fig.savefig(output / f"case4_{signal_name}_snapshot{snapshot}.png")
                    plt.close(fig)
                    csv_file = f"case4_{signal_name}_snapshot{snapshot}.csv"
                    np.savetxt(output / csv_file, np.column_stack((time_ps, values)), delimiter=",",
                               header="relative_time_ps,voltage_v", comments="", fmt="%.12g")
                    metric["snapshots"].append({"signal": signal_name, "snapshot": snapshot,
                                                 "start_sample": begin, "start_time_s": begin * dt,
                                                 "source_payload_start_bit": begin_bit - guard_bits,
                                                 "span_ui": window_bits, "csv_file": csv_file,
                                                 "figure_file": f"case4_{signal_name}_snapshot{snapshot}.png"})
        print(f"{key}: taps={np.array2string(taps, precision=6)}, empirical gap={high_min-low_max:.6g} V", flush=True)

    # Common bin edges as well as limits make probabilities comparable.
    raw_low = min(float(np.min(e[2])) for e in eye_waves)
    raw_high = max(float(np.max(e[2])) for e in eye_waves)
    margin = max((raw_high - raw_low) * 0.04, 0.002)
    vlo, vhi = raw_low - margin, raw_high + margin
    eye_cache = []
    for key, label, wave in eye_waves:
        probability, extent = eye_histogram(wave, n_per_ui, voltage_range=(vlo, vhi))
        eye_cache.append((key, label, probability, extent))
    del eye_waves
    summary["parameters"]["eye_voltage_bins"] = 400
    summary["parameters"]["eye_common_voltage_limits_v"] = [vlo, vhi]
    summary["parameters"]["eye_voltage_bin_width_v"] = (vhi - vlo) / 400
    for key, label, probability, extent in eye_cache:
        fig, ax = plt.subplots(figsize=(6.2, 4), constrained_layout=True)
        im = ax.imshow(np.log10(np.maximum(np.hstack((probability, probability)), 1e-6)),
                       origin="lower", extent=[0, 2, extent[2], extent[3]], aspect="auto",
                       cmap="viridis", vmin=-6, vmax=0, interpolation="nearest")
        ax.set(xlabel="Time (UI)", ylabel="Voltage (V)", title=label, ylim=(vlo, vhi))
        ax.axvline(0.5, color="white", alpha=0.35, lw=0.8)
        ax.axvline(1.5, color="white", alpha=0.35, lw=0.8)
        fig.colorbar(im, ax=ax, label="log10 probability per voltage bin")
        fig.savefig(output / f"eye_{key}.png")
        plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for ax, (_, label, probability, extent), metric in zip(axes.flat, eye_cache, summary["cases"]):
        im = ax.imshow(np.log10(np.maximum(np.hstack((probability, probability)), 1e-6)),
                       origin="lower", extent=[0, 2, extent[2], extent[3]], aspect="auto",
                       cmap="viridis", vmin=-6, vmax=0, interpolation="nearest")
        ax.set(title=f"{label}\nSampled gap {metric['empirical_sampled_vertical_gap_v']*1e3:.1f} mV",
               xlabel="Time (UI)", ylabel="Voltage (V)", ylim=(vlo, vhi))
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0, 0.9, "40 Gbps NRZ | 50 cm stripline\n\n0..2 V before FFE; H = S21/2\nTr = 5 ps for every transition\nSame random bits in all cases\nSample phase centered at 0.5 UI\n\nFinite noiseless sequence:\nsampled gap is not a BER estimate.", va="top", fontsize=10)
    fig.colorbar(im, ax=list(axes.flat[:5]), label="log10 probability per voltage bin", shrink=0.8)
    fig.savefig(output / "eyes_overview.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for metric in summary["cases"]:
        offsets = np.arange(-metric["npre"], metric["npost"] + 1)
        axes[0].plot(offsets, metric["taps"], "o-", label=metric["label"])
        freq = np.linspace(0, 40e9, 1001)
        response = sum(tap * np.exp(-2j * np.pi * freq * tb * j) for j, tap in enumerate(metric["taps"]))
        axes[1].plot(freq / 1e9, 20 * np.log10(np.maximum(np.abs(response), 1e-15)), label=metric["label"])
    axes[0].set(xlabel="Tap offset from main (UI)", ylabel="Normalized coefficient", title="ZF taps: sum of absolute values = 1", xticks=[-1, 0, 1, 2])
    axes[1].set(xlabel="Frequency (GHz)", ylabel="FFE magnitude (dB)", title="FFE response (reference only)")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    fig.savefig(output / "ffe_taps.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True)
    for label, eq, main_index in sbr_cache:
        lo, hi = max(0, main_index - 3 * n_per_ui), min(len(eq), main_index + 10 * n_per_ui)
        ax.plot((np.arange(lo, hi) - main_index) / n_per_ui, eq[lo:hi], label=label)
    ax.set(xlabel="Time relative to each main cursor (UI)", ylabel="Single-bit response (V)", title="Full physical response after each FFE")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(output / "equalized_sbr.png")
    plt.close(fig)

    summary = _jsonable(summary)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    columns = ["name", "npre", "npost", "taps", "main_cursor_v", "full_residual_isi_l1_v",
               "empirical_sampled_vertical_gap_v", "empirical_high_min_v", "empirical_low_max_v",
               "tx_min_v", "tx_max_v", "tx_swing_v", "sample_phase_ps", "phase_shift_from_design_ps",
               "zf_matrix_residual_v", "physical_selected_residual_v"]
    with (output / "metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for metric in summary["cases"]:
            writer.writerow({key: json.dumps(metric[key]) if key == "taps" else metric[key] for key in columns})
    return summary


def main(argv=None):
    here = Path(__file__).resolve().parent
    candidates = [here / "instructions" / "HW7_stripline_50cm.s2p", here.parent / "data" / "HW7_stripline_50cm.s2p"]
    default_data = next((p for p in candidates if p.exists()), candidates[0])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=default_data)
    parser.add_argument("--output", "--output-dir", type=Path, default=here / "results")
    parser.add_argument("--bits", type=int, default=65536)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)
    run_simulation(args.data, args.output, args.bits, args.seed)


if __name__ == "__main__":
    main()
