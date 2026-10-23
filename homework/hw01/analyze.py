#!/usr/bin/env python3
"""Reproduce HW1 Kaiser-window analysis and estimate tones from the data.

The input has a header and two columns: time in ns and voltage in V.  From the
repository root, run ``python homework/hw01/analyze.py``. Input/output defaults
are relative to the script, so it also works when invoked by an absolute path.
The fit searches the entire sampled band; ``--max-frequency-ghz`` changes only
the figure limits.

Kaiser alpha follows the course convention beta = pi * alpha.  No built-in
Kaiser window, known tone frequencies, amplitudes, or component count is used.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import rfft, rfftfreq


DEFAULT_ALPHAS = (0.1, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0)


def validate_samples(times_s, values):
    """Return finite 1-D samples and dt; reject nonuniform/invalid sampling."""
    t = np.asarray(times_s, dtype=float)
    x = np.asarray(values, dtype=float)
    if t.ndim != 1 or x.ndim != 1 or t.size != x.size or t.size < 8:
        raise ValueError("Time and signal must be equal-length 1-D arrays of at least 8 samples.")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(x)):
        raise ValueError("Time and signal must contain only finite values.")
    steps = np.diff(t)
    dt = float(np.median(steps))
    if dt <= 0 or np.any(steps <= 0):
        raise ValueError("Sample times must be strictly increasing.")
    # Also account for cancellation when a nonzero absolute time origin is used.
    atol = max(abs(dt) * 1e-8, np.finfo(float).eps * max(abs(t).max(), abs(dt)) * 16)
    if not np.allclose(steps, dt, rtol=1e-8, atol=atol):
        raise ValueError("FFT analysis requires uniformly spaced sample times.")
    return t, x, dt


def load_signal(path, time_unit="ns"):
    """Load the assignment's header + two-column text format."""
    data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError("Input must contain exactly two columns: time and voltage.")
    scale = {"ns": 1e-9, "s": 1.0}[time_unit]
    return validate_samples(data[:, 0] * scale, data[:, 1])


def kaiser_window(n_samples, alpha):
    """Evaluate the assigned Kaiser equation explicitly (beta = pi * alpha)."""
    if n_samples < 2:
        raise ValueError("A window requires at least two samples.")
    if not np.isfinite(alpha) or alpha < 0 or alpha > 50:
        raise ValueError("Alpha must be finite and between 0 and 50.")
    n = np.arange(n_samples, dtype=float)
    radius = np.sqrt(np.maximum(0.0, 1.0 - (2 * n / (n_samples - 1) - 1) ** 2))
    return np.i0(np.pi * alpha * radius) / np.i0(np.pi * alpha)


def amplitude_spectrum(values, dt_s, alpha, fft_size=None):
    """Return frequency, one-sided amplitude in V, and coherent window gain.

    Interior positive-frequency bins combine a real signal's two spectral
    halves. DC, and Nyquist for an even FFT length, must not be doubled.
    Dividing by sum(window) corrects coherent gain, not overlap/scalloping.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)):
        raise ValueError("Signal must be a finite 1-D array.")
    if not np.isfinite(dt_s) or dt_s <= 0:
        raise ValueError("Sampling interval must be positive and finite.")
    m = x.size if fft_size is None else int(fft_size)
    if m < x.size:
        raise ValueError("FFT size must not truncate the recorded samples.")
    w = kaiser_window(x.size, alpha)
    a = np.abs(rfft(x * w, n=m)) / np.sum(w)
    if m % 2 == 0:
        a[1:-1] *= 2
    else:
        a[1:] *= 2
    return rfftfreq(m, dt_s), a, float(np.mean(w))


def _basis(u, tau):
    phase = 2 * np.pi * np.outer(tau, u)
    return np.column_stack((np.ones(tau.size), np.sin(phase), np.cos(phase)))


def _project(u, tau, x):
    design = _basis(u, tau)
    coefficient = np.linalg.lstsq(design, x, rcond=None)[0]
    return design @ coefficient, coefficient


def _candidate_frequencies(residual, dt, duration, existing_u):
    """Find residual peaks with both narrow and low-sidelobe windows."""
    m = 1 << (max(4096, 32 * residual.size) - 1).bit_length()
    candidates = []
    for alpha in (0.1, 1.0, 5.0):
        freq, amplitude, _ = amplitude_spectrum(residual, dt, alpha, m)
        peaks = np.flatnonzero((amplitude[1:-1] > amplitude[:-2]) &
                               (amplitude[1:-1] > amplitude[2:])) + 1
        peaks = peaks[(peaks > 0) & (freq[peaks] < 0.5 / dt)]
        strongest = peaks[np.argsort(amplitude[peaks])[::-1]][:3]
        for i in strongest:
            u = float(freq[i] * duration)
            if all(abs(u - old) > 0.15 for old in (*existing_u, *candidates)):
                candidates.append(u)
    return candidates


def _refine_frequencies(initial, tau, signal):
    """Bounded damped Gauss-Newton for the variable-projection residual.

    Differentiate the least-squares projection itself, including its residual
    term. Frequencies use DFT-bin units, avoiding poorly scaled Hz derivatives.
    Damping and a one-bin step limit keep spectral initializations local.
    """
    u = np.asarray(initial, dtype=float).copy()
    damping = 1e-3
    for _ in range(100):
        design = _basis(u, tau)
        inverse = np.linalg.pinv(design)
        coeff = inverse @ signal
        residual = design @ coeff - signal
        phase = 2 * np.pi * np.outer(tau, u)
        dsin = 2 * np.pi * tau[:, None] * np.cos(phase)
        dcos = -2 * np.pi * tau[:, None] * np.sin(phase)
        count = u.size
        raw = dsin * coeff[1:1 + count] + dcos * coeff[1 + count:]
        correction = np.zeros((1 + 2 * count, count))
        for j in range(count):
            correction[1 + j, j] = dsin[:, j] @ residual
            correction[1 + count + j, j] = dcos[:, j] @ residual
        jacobian = raw - design @ (inverse @ raw) - inverse.T @ correction
        hessian = jacobian.T @ jacobian
        gradient = jacobian.T @ residual
        diagonal = np.maximum(np.diag(hessian), 1e-16)
        old_error = float(residual @ residual)
        accepted = False
        for _ in range(16):
            step = np.linalg.lstsq(hessian + damping * np.diag(diagonal), -gradient, rcond=None)[0]
            largest = np.max(np.abs(step))
            if largest > 1:
                step /= largest
            proposed_u = np.clip(u + step, 1e-8, tau.size / 2 - 1e-8)
            proposed, _ = _project(proposed_u, tau, signal)
            new_error = float(np.sum((proposed - signal) ** 2))
            if new_error < old_error:
                u = proposed_u
                damping = max(damping / 3, 1e-12)
                accepted = True
                break
            damping *= 10
        if not accepted or np.max(np.abs(step)) < 1e-12:
            break
    return np.sort(u)


def fit_tones(times_s, values, max_components=6, relative_tolerance=1e-10):
    """Estimate count, frequencies, amplitudes and phases without answer seeds.

    Add frequencies found in a residual spectrum. At each order, try several
    candidates from alpha 0.1, 1, and 5; refine all frequencies together with
    bounded damped Gauss-Newton. The offset and sine/cosine coefficients
    are solved by linear least squares at each frequency trial (variable
    projection). Choose the lowest BIC, counting 3 parameters per tone plus
    the offset. A declared relative residual tolerance stops fitting floating
    point roundoff in essentially noiseless data. The model is for stationary
    sinusoidal components; model selection is not a universal tone detector.
    """
    t, x, dt = validate_samples(times_s, values)
    if not isinstance(max_components, int) or not 1 <= max_components <= 20:
        raise ValueError("max_components must be an integer between 1 and 20.")
    if not 0 < relative_tolerance < 1:
        raise ValueError("relative_tolerance must lie strictly between 0 and 1.")
    duration = x.size * dt
    tau = (t - t[0]) / duration
    scale = max(float(np.sqrt(np.mean(x * x))), np.finfo(float).tiny)
    normalized = x / scale
    u = np.empty(0)
    prediction, coefficient = _project(u, tau, normalized)
    records = []
    models = []
    stop_reason = "maximum component count reached"
    for order in range(max_components + 1):
        residual = normalized - prediction
        rss = float(np.dot(residual, residual))
        physical_residual = x - prediction * scale
        physical_rss = float(physical_residual @ physical_residual)
        floor = x.size * relative_tolerance**2
        bic = float(x.size * np.log(max(rss, floor) / x.size) + (3 * order + 1) * np.log(x.size))
        records.append({"q": order, "rss_v2": physical_rss, "bic": bic,
                        "rmse_v": float(np.sqrt(physical_rss / x.size))})
        models.append((u.copy(), coefficient.copy(), prediction.copy()))
        if np.sqrt(rss / x.size) <= relative_tolerance:
            stop_reason = "relative numerical residual tolerance reached"
            break
        if order == max_components:
            break
        # Two increasing orders guard against spending work on noise peaks.
        if order >= 2 and records[-1]["bic"] > records[-2]["bic"] > records[-3]["bic"]:
            stop_reason = "BIC increased for two successive model orders"
            break
        candidates = _candidate_frequencies(residual, dt, duration, u)
        if not candidates:
            stop_reason = "no additional interior spectral peaks"
            break
        best = None
        for candidate in candidates:
            initial = np.sort(np.append(u, candidate))
            refined = _refine_frequencies(initial, tau, normalized)
            proposed, coeff = _project(refined, tau, normalized)
            error = float(np.sum((proposed - normalized) ** 2))
            if best is None or error < best[0]:
                best = (error, refined, proposed, coeff)
        _, u, prediction, coefficient = best
    selected = int(np.argmin([record["bic"] for record in records]))
    u, coefficient, prediction = models[selected]
    prediction = prediction * scale
    coefficient = coefficient * scale
    q = u.size
    components = []
    for i, frequency in enumerate(u / duration):
        sine, cosine = coefficient[1 + i], coefficient[1 + q + i]
        amplitude = float(np.hypot(sine, cosine))
        components.append({"index": i + 1, "frequency_hz": float(frequency),
                           "frequency_mhz": float(frequency / 1e6),
                           "amplitude_v": amplitude,
                           "phase_rad": float(np.arctan2(cosine, sine)),
                           "amplitude_dbv": float(20 * np.log10(max(amplitude, np.finfo(float).tiny)))})
    residual = x - prediction
    result = {"q": int(q), "offset_v": float(coefficient[0]),
              "time_origin_s": float(t[0]),
              "rmse_v": float(np.sqrt(np.mean(residual**2))),
              "max_abs_residual_v": float(np.max(np.abs(residual))),
              "relative_rmse": float(np.sqrt(np.mean(residual**2)) / scale),
              "bic": records[selected]["bic"], "components": components,
              "order_selection": records, "stop_reason": stop_reason}
    return result, prediction


def _db(amplitude):
    return 20 * np.log10(np.maximum(amplitude, 1e-16))


def _save(fig, path):
    fig.savefig(path, dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def create_figures(t, x, fit, prediction, output, alphas, fft_size, max_frequency_ghz):
    """Save deterministic, headless figures and return their relative names."""
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True,
                         "grid.alpha": 0.22, "figure.constrained_layout.use": True})
    times = (t - t[0]) * 1e9
    dt = float(np.median(np.diff(t)))
    files = {}
    fig, ax = plt.subplots(figsize=(9.2, 3.5))
    ax.plot(times, x, color="#24599b", linewidth=1.15)
    ax.set(xlabel="Time from first sample (ns)", ylabel="Voltage (V)",
           title=f"Recorded signal: {x.size} samples")
    files["time_domain"] = "time_domain.png"
    _save(fig, output / files["time_domain"])

    fig, ax = plt.subplots(figsize=(9.2, 3.7))
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.9, len(alphas)))
    for alpha, color in zip(alphas, colors):
        ax.plot(times, kaiser_window(x.size, alpha), label=f"alpha = {alpha:g}", color=color)
    ax.set(xlabel="Time from first sample (ns)", ylabel="Window value", title="Assigned Kaiser windows")
    ax.legend(ncol=4, fontsize=8)
    files["windows"] = "kaiser_windows.png"
    _save(fig, output / files["windows"])

    fig_all, ax_all = plt.subplots(figsize=(9.2, 4.5))
    spectrum_files = []
    gains = []
    for alpha, color in zip(alphas, colors):
        frequency, amplitude, gain = amplitude_spectrum(x, dt, alpha, fft_size)
        mask = frequency <= min(max_frequency_ghz * 1e9, 0.5 / dt)
        fg, db = frequency[mask] / 1e9, _db(amplitude[mask])
        ax_all.plot(fg, db, color=color, label=f"alpha = {alpha:g}", linewidth=1.0)
        fig, ax = plt.subplots(figsize=(9.2, 4.3))
        ax.plot(fg, db, color="#24599b", linewidth=1.1)
        ax.set(xlabel="Frequency (GHz)", ylabel="Amplitude (dBV; reference 1 V)",
               xlim=(0, min(max_frequency_ghz, 0.5 / dt / 1e9)), ylim=(-150, 5),
               title=f"Kaiser alpha = {alpha:g}; coherent gain = {gain:.6f}")
        for component in fit["components"]:
            ax.axvline(component["frequency_hz"] / 1e9, color="#555555", alpha=0.4,
                       linewidth=0.7, linestyle="--")
        tag = f"{alpha:g}".replace(".", "p")
        name = f"spectrum_alpha_{tag}.png"
        spectrum_files.append(name)
        gains.append({"alpha": float(alpha), "coherent_gain": gain})
        _save(fig, output / name)
    ax_all.set(xlabel="Frequency (GHz)", ylabel="Amplitude (dBV; reference 1 V)",
               xlim=(0, min(max_frequency_ghz, 0.5 / dt / 1e9)), ylim=(-150, 5),
               title="Window tradeoff: main-lobe width and sidelobe suppression")
    ax_all.legend(ncol=4, fontsize=8, loc="lower left")
    files["spectra"] = "spectra.png"
    _save(fig_all, output / files["spectra"])
    files["individual_spectra"] = spectrum_files

    fig, (ax, residual_ax) = plt.subplots(2, 1, figsize=(9.2, 6.0), sharex=True,
                                         gridspec_kw={"height_ratios": (2, 1)})
    ax.plot(times, x, label="Recorded", linewidth=1.8, color="#24599b", alpha=0.75)
    ax.plot(times, prediction, label="Estimated-tone reconstruction", linewidth=1.0,
            color="#de7c30", linestyle="--")
    ax.set(ylabel="Voltage (V)", title=f"Data-derived reconstruction: Q = {fit['q']}")
    ax.legend(fontsize=9)
    residual_ax.plot(times, x - prediction, linewidth=0.8, color="#6d3987")
    residual_ax.set(xlabel="Time from first sample (ns)", ylabel="Residual (V)",
                    title=f"Recorded minus reconstructed; RMSE = {fit['rmse_v']:.3e} V")
    residual_ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    files["reconstruction"] = "reconstruction.png"
    _save(fig, output / files["reconstruction"])
    return files, gains


def main(argv=None):
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=here / "instructions" / "HW1_Kaiser.txt")
    parser.add_argument("--output", type=Path, default=here / "results")
    parser.add_argument("--time-unit", choices=("ns", "s"), default="ns")
    parser.add_argument("--alphas", type=float, nargs="+", default=list(DEFAULT_ALPHAS))
    parser.add_argument("--fft-size", type=int, default=2**20 - 1)
    parser.add_argument("--max-frequency-ghz", type=float, default=2.0)
    parser.add_argument("--max-components", type=int, default=6)
    args = parser.parse_args(argv)
    if args.max_frequency_ghz <= 0 or not np.isfinite(args.max_frequency_ghz):
        parser.error("--max-frequency-ghz must be positive and finite")
    for alpha in args.alphas:
        kaiser_window(2, alpha)
    t, x, dt = load_signal(args.input, args.time_unit)
    if args.fft_size < x.size:
        parser.error("--fft-size must be at least the sample count")
    fit, prediction = fit_tones(t, x, args.max_components)
    args.output.mkdir(parents=True, exist_ok=True)
    figures, gains = create_figures(t, x, fit, prediction, args.output,
                                   args.alphas, args.fft_size, args.max_frequency_ghz)
    fields = ["index", "frequency_hz", "frequency_mhz", "amplitude_v", "phase_rad", "amplitude_dbv"]
    with (args.output / "parameters.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(fit["components"])
    summary = {
        "input_file": args.input.name,
        "sampling": {"n_samples": int(x.size), "dt_s": dt, "fs_hz": 1 / dt,
                     "nyquist_hz": 0.5 / dt, "sample_span_s": float(t[-1] - t[0]),
                     "record_duration_s": x.size * dt, "unpadded_df_hz": 1 / (x.size * dt),
                     "fft_size": args.fft_size, "plotted_df_hz": 1 / (args.fft_size * dt)},
        "analysis": {"alphas": args.alphas, "window_convention": "beta = pi * alpha",
                     "window_gains": gains, "max_components": args.max_components,
                     "numerical_relative_tolerance": 1e-10,
                     "fit_method": "Residual Kaiser spectral candidates; variable projection with bounded damped Gauss-Newton frequency refinement; BIC order selection",
                     "fit_search_band_hz": [0.0, 0.5 / dt],
                     "display_max_frequency_ghz": args.max_frequency_ghz,
                     "amplitude_reference": "dBV = 20 log10(amplitude / 1 V)",
                     "model": "offset + sum A_i sin(2 pi f_i (t - time_origin) + phase_i)"},
        "fit": fit,
        "artifacts": {**figures, "parameters": "parameters.csv", "summary": "summary.json"},
    }
    with (args.output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(f"N = {x.size}; dt = {dt:.12g} s; fs = {1 / dt:.12g} Hz")
    print(f"Selected Q = {fit['q']}; reconstruction RMSE = {fit['rmse_v']:.6e} V")
    for component in fit["components"]:
        print(f"  {component['frequency_mhz']:.9f} MHz: {component['amplitude_v']:.12g} V")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
