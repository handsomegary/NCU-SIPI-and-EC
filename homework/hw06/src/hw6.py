"""Reproducible peak distortion analysis for the supplied 40 Gbit/s channel.

Only NumPy and Matplotlib are required. The default input/output paths are
relative to this homework directory, so the analysis works after any clone.
The source is unipolar, 0 to 2 V, and the matched transfer is S21 / 2.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def read_touchstone(path):
    """Read a Touchstone 1.x two-port S-parameter file (MA, DB, or RI).

    Two-port order is S11, S21, S12, S22. Comments, continuation lines,
    Fortran D exponents, and HZ/KHZ/MHZ/GHZ units are supported. Unsupported
    Touchstone 2 keywords or non-S parameters are rejected explicitly.
    Returns frequencies in Hz, complex S21, and the real reference impedance.
    """
    units = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}
    unit, parameter, encoding, resistance = "GHZ", "S", "MA", 50.0
    tokens = []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        line = line.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            raise ValueError("Only Touchstone 1.x two-port files are supported")
        if line.startswith("#"):
            fields = line[1:].upper().split()
            if len(fields) < 3:
                raise ValueError("Invalid Touchstone option line")
            unit, parameter, encoding = fields[:3]
            if "R" in fields:
                resistance = float(fields[fields.index("R") + 1])
            continue
        tokens.extend(float(t.replace("D", "E").replace("d", "e"))
                      for t in line.split())
    if unit not in units or parameter != "S" or encoding not in ("MA", "DB", "RI"):
        raise ValueError("Unsupported Touchstone units, parameter, or encoding")
    if not tokens or len(tokens) % 9:
        raise ValueError("Expected nine numeric fields per two-port record")
    data = np.asarray(tokens, dtype=np.float64).reshape(-1, 9)
    f = data[:, 0] * units[unit]
    a, b = data[:, 3], data[:, 4]
    if encoding == "RI":
        s21 = a + 1j * b
    else:
        mag = a if encoding == "MA" else 10.0 ** (a / 20.0)
        s21 = mag * np.exp(1j * np.deg2rad(b))
    if not np.all(np.isfinite(data)) or np.any(np.diff(f) <= 0):
        raise ValueError("Frequency records must be finite and strictly increasing")
    return f, s21, resistance


def read_bits(path):
    """Read binary characters, accepting only whitespace besides 0 and 1."""
    text = "".join(Path(path).read_text(encoding="utf-8-sig").split())
    if not text or set(text) - {"0", "1"}:
        raise ValueError("Bit pattern must contain only 0, 1, and whitespace")
    return np.frombuffer(text.encode("ascii"), dtype=np.uint8) - ord("0")


def single_bit_response(f, s21, *, samples_per_ui=250, bit_period=25e-12,
                        rise_time=5e-12, vpp=2.0, kaiser_alpha=1.0):
    """Return the complete real IDFT response and its sampling metadata.

    The causal unit-bit trapezoid is the convolution of a [0,Tb] rectangle
    and a normalized [0,Tr] rectangle. Its transform is
    Vpp*Tb*sinc(f*Tb)*sinc(f*Tr)*exp(-j*pi*f*(Tb+Tr)).
    The one-sided physical spectrum is inverted with irfft and scaled by
    M*df; this counts positive/negative frequencies once and DC only once.
    Bins beyond the measured 80 GHz band are zero, not extrapolated.
    """
    f, s21 = np.asarray(f, float), np.asarray(s21, complex)
    if f.ndim != 1 or len(f) < 2 or s21.shape != f.shape:
        raise ValueError("Expected matching one-dimensional frequency and S21 arrays")
    if samples_per_ui < 2 or int(samples_per_ui) != samples_per_ui:
        raise ValueError("samples_per_ui must be an integer of at least two")
    if not (bit_period > 0 and 0 < rise_time <= bit_period and vpp > 0):
        raise ValueError("Require Tb>0, 0<Tr<=Tb, and Vpp>0")
    if kaiser_alpha < 0 or not np.isfinite(kaiser_alpha):
        raise ValueError("Kaiser alpha must be finite and nonnegative")
    df = float(f[1] - f[0])
    if f[0] != 0 or df <= 0 or not np.allclose(np.diff(f), df, rtol=1e-8, atol=1e-3):
        raise ValueError("IDFT requires a uniform frequency grid beginning at DC")
    dt = bit_period / int(samples_per_ui)
    m_float = 1.0 / (dt * df)
    m = round(m_float)
    if not np.isclose(m, m_float, rtol=1e-10):
        raise ValueError("The time step must divide the frequency-grid period")
    if f[-1] >= 0.5 / dt:
        raise ValueError("Measured frequencies must be strictly below Nyquist")
    window = np.i0(np.pi * kaiser_alpha * np.sqrt(np.clip(1 - (f / f[-1]) ** 2, 0, 1)))
    window /= np.i0(np.pi * kaiser_alpha)
    vin = (vpp * bit_period * np.sinc(f * bit_period) * np.sinc(f * rise_time)
           * np.exp(-1j * np.pi * f * (bit_period + rise_time)))
    spectrum = np.zeros(m // 2 + 1, dtype=np.complex128)
    spectrum[:len(f)] = vin * s21 * window / 2.0
    if abs(spectrum[0].imag) > 1e-15:
        raise ValueError("DC must be real for a real-valued channel response")
    full_sbr = np.fft.irfft(spectrum, n=m) * (m * df)
    meta = {"dt_s": dt, "df_hz": df, "idft_points": m,
            "idft_period_s": m * dt, "frequency_points": len(f),
            "maximum_frequency_hz": float(f[-1]),
            "s21_dc": float(s21[0].real),
            "matched_dc_gain": float(s21[0].real / 2.0)}
    return full_sbr, dt, meta


def pda(sbr, samples_per_ui, tolerance=0.0002):
    """Compute two UI of conditional worst-case boundaries around SBR peak.

    wc1/wc0 omit interferers smaller than tolerance. wc1_exact/wc0_exact
    include every interferer of the supplied truncated SBR. The latter are
    rigorous for that finite deterministic model. The omissions quantify
    the relaxation from applying a tolerance; they are not random noise.
    """
    p = np.asarray(sbr, dtype=np.float64)
    n = int(samples_per_ui)
    if p.ndim != 1 or p.size < 2 * n or n < 1 or n != samples_per_ui:
        raise ValueError("Invalid SBR or samples_per_ui")
    if tolerance < 0 or not np.isfinite(tolerance) or not np.all(np.isfinite(p)):
        raise ValueError("SBR and tolerance must be finite; tolerance must be nonnegative")
    kmain = int(np.argmax(p))
    kstart = kmain - n
    if kstart < 0 or kstart + 2 * n > len(p):
        raise ValueError("SBR must contain one UI on both sides of its peak")
    result = {"kmain": kmain, "kstart": kstart,
              "main": p[kstart:kstart + 2 * n].copy(),
              "tolerance_v": float(tolerance)}
    names = ("wc0", "wc1", "wc0_exact", "wc1_exact", "omitted_positive",
             "omitted_negative", "outer_min", "outer_max")
    result.update({key: np.empty(2 * n) for key in names})
    for q, ix in enumerate(range(kstart, kstart + 2 * n)):
        # Every in-range index with this phase is a symbol-spaced cursor.
        # Excluding the main index explicitly avoids the negative-slice and
        # omitted-index-zero hazards of a reverse slice.
        indices = np.arange(ix % n, len(p), n)
        taps = p[indices[indices != ix]]
        positive, negative = taps[taps > 0], taps[taps < 0]
        large_positive = positive[positive >= tolerance]
        large_negative = negative[negative <= -tolerance]
        result["wc0_exact"][q] = np.sum(positive)
        result["wc1_exact"][q] = p[ix] + np.sum(negative)
        result["wc0"][q] = np.sum(large_positive)
        result["wc1"][q] = p[ix] + np.sum(large_negative)
        result["omitted_positive"][q] = np.sum(positive) - np.sum(large_positive)
        result["omitted_negative"][q] = np.sum(large_negative) - np.sum(negative)
        result["outer_min"][q] = np.sum(negative) + min(p[ix], 0.0)
        result["outer_max"][q] = np.sum(positive) + max(p[ix], 0.0)
    return result


def sectional_convolution(bits, sbr, samples_per_ui, block_bits=1024):
    """FFT overlap-add of an NRZ bit impulse train with the channel SBR.

    The SBR already includes equal linear rise/fall shaping and Vpp, so each
    bit is a unit impulse. Adjacent shaped pulses add to a continuous NRZ
    signal. This is equivalent to constructing that NRZ waveform and then
    convolving with its linear channel, without shaping or scaling twice.
    Every block is zero-padded for *linear* convolution; the last partial
    block is retained. The output includes the complete finite-model tail.
    """
    bits, sbr = np.asarray(bits), np.asarray(sbr, dtype=np.float64)
    n, block_bits = int(samples_per_ui), int(block_bits)
    if bits.ndim != 1 or not len(bits) or not np.all((bits == 0) | (bits == 1)):
        raise ValueError("bits must be a nonempty one-dimensional binary array")
    if n < 1 or block_bits < 1 or sbr.ndim != 1 or not len(sbr):
        raise ValueError("Invalid convolution configuration")
    block_samples = min(block_bits, len(bits)) * n
    fft_size = 1 << (block_samples + len(sbr) - 2).bit_length()
    h = np.fft.rfft(sbr, n=fft_size)
    output = np.zeros(len(bits) * n + len(sbr) - 1, dtype=np.float64)
    for first in range(0, len(bits), block_bits):
        current = bits[first:first + block_bits]
        impulse_train = np.zeros(len(current) * n, dtype=np.float64)
        impulse_train[::n] = current
        block_output = np.fft.irfft(np.fft.rfft(impulse_train, n=fft_size) * h,
                                     n=fft_size)
        length = len(impulse_train) + len(sbr) - 1
        start = first * n
        output[start:start + length] += block_output[:length]
    return output


def analyze_eye(y, bits, pda_result, samples_per_ui, voltage_edges,
                max_traces=512, sbr_length=None):
    """Fold one UI, duplicate only when plotting, and compare conditional extrema.

    Trace j begins at y[kstart + j*N]. Its underlying reference bit is
    bits[j]. For second-UI conditional extrema the same bit labels are
    retained, rather than relabeling to the next bit. This exactly matches
    PDA's fixed main bit throughout its two-UI time window.
    """
    n, bits = int(samples_per_ui), np.asarray(bits, dtype=np.uint8)
    if sbr_length is None:
        raise ValueError("sbr_length is required to exclude finite-record startup and tail")
    guard = math.ceil(sbr_length / n) + 2
    first, stop = guard, len(bits) - guard
    count = stop - first
    if count < 2 or np.unique(bits[first:stop]).size < 2:
        raise ValueError("Need both bit levels after excluding the finite-record guards")
    start = int(pda_result["kstart"]) + first * n
    if start < 0 or start + (count + 1) * n > len(y):
        raise ValueError("Output waveform is too short for the aligned eye windows")
    matrix = np.asarray(y[start:start + count * n]).reshape(count, n)
    next_matrix = np.asarray(y[start + n:start + (count + 1) * n]).reshape(count, n)
    labels = bits[first:stop].astype(bool)
    observed1 = np.concatenate((matrix[labels].min(axis=0), next_matrix[labels].min(axis=0)))
    observed0 = np.concatenate((matrix[~labels].max(axis=0), next_matrix[~labels].max(axis=0)))
    hist = np.empty((len(voltage_edges) - 1, n), dtype=np.uint32)
    for phase in range(n):
        hist[:, phase] = np.histogram(matrix[:, phase], bins=voltage_edges)[0]
    total = hist.sum(axis=0, dtype=np.uint64)
    if not np.all(total == count):
        raise ValueError("Voltage histogram range excludes waveform samples")
    probability = hist.astype(np.float64) / total[None, :]
    trace_indices = np.linspace(0, count - 1, min(int(max_traces), count), dtype=int)
    violation1 = np.maximum(pda_result["wc1_exact"] - observed1, 0)
    violation0 = np.maximum(observed0 - pda_result["wc0_exact"], 0)
    best = int(np.argmax(pda_result["wc1_exact"] - pda_result["wc0_exact"]))
    summary = {"input_bits": len(bits), "input_ones": int(bits.sum()),
               "stable_bits": count, "stable_ones": int(labels.sum()),
               "guard_bits_per_end": guard, "histogram_samples_per_phase": count,
               "histogram_probability_max_column_error": float(np.max(abs(probability.sum(axis=0) - 1))),
               "stable_min_voltage_v": float(matrix.min()),
               "stable_max_voltage_v": float(matrix.max()),
               "maximum_exact_bound_violation_v": float(max(violation1.max(), violation0.max())),
               "observed_min_one_at_exact_optimum_v": float(observed1[best]),
               "observed_max_zero_at_exact_optimum_v": float(observed0[best]),
               "observed_gap_at_exact_optimum_v": float(observed1[best] - observed0[best]),
               "one_clearance_at_exact_optimum_v": float(observed1[best] - pda_result["wc1_exact"][best]),
               "zero_clearance_at_exact_optimum_v": float(pda_result["wc0_exact"][best] - observed0[best])}
    return {"histogram": hist, "probability": probability,
            "observed_min_one": observed1, "observed_max_zero": observed0,
            "traces": matrix[trace_indices].copy(), "summary": summary}


def _pda_summary(boundaries, dt, n):
    exact_gap = boundaries["wc1_exact"] - boundaries["wc0_exact"]
    filtered_gap = boundaries["wc1"] - boundaries["wc0"]
    exact_best, filtered_best = int(np.argmax(exact_gap)), int(np.argmax(filtered_gap))
    result = {}
    for prefix, gap, best, one, zero in (
            ("exact", exact_gap, exact_best, "wc1_exact", "wc0_exact"),
            ("filtered", filtered_gap, filtered_best, "wc1", "wc0")):
        ix = int(boundaries["kstart"]) + best
        result[prefix] = {"maximum_vertical_gap_v": float(gap[best]),
                          "optimum_boundary_index": best,
                          "optimum_relative_time_ui": best / n,
                          "optimum_sbr_time_ps": ix * dt * 1e12,
                          "optimum_global_phase_ps": (ix % n) * dt * 1e12,
                          "minimum_one_v": float(boundaries[one][best]),
                          "maximum_zero_v": float(boundaries[zero][best])}
    omitted = boundaries["omitted_positive"] + boundaries["omitted_negative"]
    result["tolerance_v"] = boundaries["tolerance_v"]
    result["maximum_sum_of_omitted_taps_v"] = float(omitted.max())
    result["exact_model_scope"] = "All taps of the finite, bandwidth-limited, windowed, truncated SBR; no jitter or noise"
    return result


def _plot_figures(output, full_sbr, sbr, dt, n, boundaries, eyes, voltage_edges, truncation_time):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 130, "savefig.dpi": 180})
    t = np.arange(len(full_sbr)) * dt * 1e9
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    ax[0].plot(t, full_sbr, color="#2269a1", linewidth=1)
    ax[0].set(xlim=(0, 3), xlabel="Time (ns)", ylabel="Voltage (V)", title="Single-bit response")
    decimate = max(1, len(full_sbr) // 30000)
    ax[1].plot(t[::decimate], abs(full_sbr[::decimate]), color="#2269a1", linewidth=.8)
    ax[1].axvline(truncation_time * 1e9, color="#b84c2e", linestyle="--", label="Truncation")
    ax[1].set(xlim=(0, 30), yscale="log", ylim=(1e-9, 2), xlabel="Time (ns)",
              ylabel="Absolute voltage (V)", title="Tail and truncation")
    ax[1].legend(frameon=False)
    for a in ax:
        a.grid(alpha=.2)
    fig.savefig(output / "sbr.png")
    plt.close(fig)

    x = np.arange(2 * n) / n
    fig, ax = plt.subplots(figsize=(7.5, 4), layout="constrained")
    for key, label, color in (("wc1_exact", "WC1: all retained taps", "#ae2d40"),
                               ("wc0_exact", "WC0: all retained taps", "#12629b")):
        ax.plot(x, boundaries[key], color=color, linewidth=1.7, label=label)
    ax.plot(x, boundaries["wc1"], color="#ae2d40", linestyle="--", linewidth=1,
            label=f"WC1: {boundaries['tolerance_v']*1e3:g} mV tap tolerance")
    ax.plot(x, boundaries["wc0"], color="#12629b", linestyle="--", linewidth=1,
            label=f"WC0: {boundaries['tolerance_v']*1e3:g} mV tap tolerance")
    ax.set(xlim=(0, 2), xlabel="Time (UI; SBR peak at 1 UI)", ylabel="Voltage (V)",
           title="PDA inner boundaries")
    ax.grid(alpha=.2)
    ax.legend(loc="upper center", fontsize=8, frameon=False)
    fig.savefig(output / "pda.png")
    plt.close(fig)

    def draw_eye(ax, eye, statistical):
        if statistical:
            # Include phase zero at the right edge and place pixel centers
            # on r/N, exactly where the waveform/PDA samples are plotted.
            probability = np.hstack((eye["probability"], eye["probability"], eye["probability"][:, :1]))
            image = ax.imshow(np.ma.masked_equal(probability, 0), origin="lower", aspect="auto",
                              extent=(-.5/n, 2+.5/n, voltage_edges[0], voltage_edges[-1]),
                              cmap="viridis", norm=LogNorm(vmin=1 / eye["summary"]["stable_bits"], vmax=.15),
                              interpolation="nearest")
        else:
            phase = np.arange(n) / n
            for offset in (0, 1):
                ax.plot(phase + offset, eye["traces"].T, color="#385e8d", alpha=.13, linewidth=.4)
            image = None
        ax.plot(x, boundaries["wc1_exact"], color="#ef5d75", linewidth=1.7, label="PDA WC1")
        ax.plot(x, boundaries["wc0_exact"], color="#f1a43c", linewidth=1.7, label="PDA WC0")
        ax.set(xlim=(0, 2), ylim=(voltage_edges[0], voltage_edges[-1]),
               xlabel="Time (UI; SBR peak at 1 UI)", ylabel="Voltage (V)")
        return image

    for name, eye in eyes.items():
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
        draw_eye(axes[0], eye, False)
        image = draw_eye(axes[1], eye, True)
        axes[0].set_title(f"{name.title()}: duplicated 1-UI traces")
        axes[1].set_title(f"{name.title()}: fixed-rate statistical eye")
        axes[0].legend(fontsize=8, frameon=False)
        fig.colorbar(image, ax=axes[1], label="Probability per voltage bin")
        fig.savefig(output / f"{name}_eye.png")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for ax, (name, eye) in zip(axes, eyes.items()):
        image = draw_eye(ax, eye, True)
        ax.set_title(f"{name.title()}: {eye['summary']['stable_bits']:,} stable bits")
    fig.colorbar(image, ax=axes, label="Probability per voltage bin", shrink=.9)
    fig.savefig(output / "eye_comparison.png")
    plt.close(fig)


def run_analysis(args):
    """Run the homework, write reproducible results, and return its summary."""
    channel_path, pattern_path = Path(args.channel), Path(args.pattern)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    f, s21, resistance = read_touchstone(channel_path)
    bits_pattern = read_bits(pattern_path)
    full_sbr, dt, spectrum_meta = single_bit_response(
        f, s21, samples_per_ui=args.samples_per_ui, bit_period=args.bit_period_ps * 1e-12,
        rise_time=args.rise_time_ps * 1e-12, vpp=args.vpp, kaiser_alpha=args.kaiser_alpha)
    n = args.samples_per_ui
    length = round(args.truncation_ns * 1e-9 / dt) + 1
    if not (2 * n < length < len(full_sbr)):
        raise ValueError("Truncation must fit inside the IDFT period and contain the SBR peak")
    sbr = full_sbr[:length].copy()
    boundaries = pda(sbr, n, args.tolerance)
    low = math.floor((boundaries["outer_min"].min() - .025) / .05) * .05
    high = math.ceil((boundaries["outer_max"].max() + .025) / .05) * .05
    edges = np.linspace(low, high, args.voltage_bins + 1)
    random_bits = np.random.default_rng(args.seed).integers(0, 2, size=len(bits_pattern), dtype=np.uint8)
    eyes = {}
    for name, bits in (("random", random_bits), ("pattern", bits_pattern)):
        print(f"Convolving {name}: {len(bits):,} bits", flush=True)
        waveform = sectional_convolution(bits, sbr, n, args.block_bits)
        eyes[name] = analyze_eye(waveform, bits, boundaries, n, edges,
                                 max_traces=args.max_traces, sbr_length=len(sbr))
        if eyes[name]["summary"]["maximum_exact_bound_violation_v"] > 1e-10:
            raise AssertionError(f"{name} waveform violates the finite-model PDA bounds")
        del waveform

    convergence = []
    previous = None
    for t_ns in (3.0, 5.0, 10.0, 20.0):
        cut = round(t_ns * 1e-9 / dt) + 1
        if cut >= len(full_sbr):
            continue
        trial = pda(full_sbr[:cut], n, args.tolerance)
        item = {"truncation_ns": t_ns, **_pda_summary(trial, dt, n)}
        if previous is not None:
            item["maximum_exact_boundary_change_from_previous_v"] = float(max(
                np.max(abs(trial["wc1_exact"] - previous["wc1_exact"])),
                np.max(abs(trial["wc0_exact"] - previous["wc0_exact"]))))
        convergence.append(item)
        previous = trial
    # Separate interior late tails from the periodic wraparound near 100 ns.
    late_end = min(len(full_sbr), round(0.9 * len(full_sbr)))
    beyond = full_sbr[length:late_end]
    summary = {
        "model": "40 Gbit/s unipolar linear NRZ; matched transfer S21/2; equal rise/fall ramps",
        "config": {"bit_period_ps": args.bit_period_ps, "rise_time_ps": args.rise_time_ps,
                   "vpp_v": args.vpp, "samples_per_ui": n, "kaiser_alpha": args.kaiser_alpha,
                   "truncation_ns": args.truncation_ns, "sbr_samples": len(sbr),
                   "tolerance_v": args.tolerance, "block_bits": args.block_bits,
                   "random_seed": args.seed, "voltage_bins": args.voltage_bins,
                   "max_displayed_traces": args.max_traces, "reference_impedance_ohm": resistance},
        "spectrum": spectrum_meta,
        "sbr": {"peak_voltage_v": float(sbr.max()), "peak_index": boundaries["kmain"],
                "peak_time_ps": boundaries["kmain"] * dt * 1e12,
                "pda_start_index": boundaries["kstart"],
                "pda_start_time_ps": boundaries["kstart"] * dt * 1e12,
                "maximum_interior_tail_after_truncation_v": float(abs(beyond).max()) if len(beyond) else 0,
                "tail_inspection_end_ns": late_end * dt * 1e9},
        "pda": _pda_summary(boundaries, dt, n),
        "truncation_study": convergence,
        "eye_alignment": {"fold_origin_sample": boundaries["kstart"],
                          "pda_relative_time_formula": "q/samples_per_ui; q=0..2*samples_per_ui-1",
                          "fold_formula": "y[kstart+j*samples_per_ui+r], r=0..samples_per_ui-1",
                          "display": "Construct one UI, then duplicate the same traces/histogram into the second UI"},
        "cases": {name: eye["summary"] for name, eye in eyes.items()},
        "source_sha256": {channel_path.name: hashlib.sha256(channel_path.read_bytes()).hexdigest(),
                          pattern_path.name: hashlib.sha256(pattern_path.read_bytes()).hexdigest()},
        "limitations": ["PDA is exact for the finite SBR used by convolution; it is not an exact bound on the untruncated physical channel.",
                        "Tap tolerance relaxes the boundary; solid overlays include all retained taps.",
                        "Finite random and supplied bit records need not attain the worst-case symbol combinations.",
                        "No jitter, random noise, receiver threshold model, or nonlinear device model is included."]}
    arrays = {"full_sbr": full_sbr, "sbr": sbr, "dt_s": np.array(dt),
              "frequency_hz": f, "s21": s21, "voltage_edges_v": edges,
              "random_bits": random_bits, "pattern_bits": bits_pattern,
              "pda_time_ui": np.arange(2 * n) / n}
    for key in ("main", "wc0", "wc1", "wc0_exact", "wc1_exact", "omitted_positive",
                "omitted_negative", "outer_min", "outer_max"):
        arrays[f"pda_{key}"] = boundaries[key]
    for name, eye in eyes.items():
        for key in ("histogram", "probability", "observed_min_one", "observed_max_zero", "traces"):
            arrays[f"{name}_{key}"] = eye[key]
    (output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    np.savez_compressed(output / "arrays.npz", **arrays)
    if not args.no_plots:
        _plot_figures(output, full_sbr, sbr, dt, n, boundaries, eyes, edges, args.truncation_ns * 1e-9)
    print(f"Exact finite-model maximum eye gap: {summary['pda']['exact']['maximum_vertical_gap_v']:.9f} V", flush=True)
    return summary


def argument_parser():
    base = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channel", type=Path, default=base / "instructions" / "HW6_channel.s2p")
    parser.add_argument("--pattern", type=Path, default=base / "instructions" / "HW6_wcbitpattern.txt")
    parser.add_argument("--output", type=Path, default=base / "results")
    parser.add_argument("--samples-per-ui", type=int, default=250)
    parser.add_argument("--bit-period-ps", type=float, default=25.0)
    parser.add_argument("--rise-time-ps", type=float, default=5.0)
    parser.add_argument("--vpp", type=float, default=2.0)
    parser.add_argument("--kaiser-alpha", type=float, default=1.0)
    parser.add_argument("--truncation-ns", type=float, default=10.0)
    parser.add_argument("--tolerance", type=float, default=0.0002)
    parser.add_argument("--block-bits", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--voltage-bins", type=int, default=512)
    parser.add_argument("--max-traces", type=int, default=512)
    parser.add_argument("--no-plots", action="store_true", help="Write numerical outputs without importing Matplotlib")
    return parser


def main(argv=None):
    args = argument_parser().parse_args(argv)
    if args.voltage_bins < 2 or args.max_traces < 1 or args.block_bits < 1:
        raise ValueError("Need at least two voltage bins and positive trace/block sizes")
    run_analysis(args)


if __name__ == "__main__":
    main()
