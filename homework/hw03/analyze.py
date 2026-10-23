"""NCU EE8091 HW3: bounded-memory overlap-add simulation and ADS comparison.

Only NumPy is needed for numerical work. Matplotlib is needed for figures.
The original assignment files are read-only inputs. No identifying information
is read from a student report, or written into generated artifacts.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
from pathlib import Path
import platform
import time

import numpy as np


def integer_points(duration_s: float, dt_s: float, label: str) -> int:
    """Use rounded integer counts and reject a duration not on the sample grid."""
    ratio = duration_s / dt_s
    count = round(ratio)
    if count < 1 or not math.isclose(ratio, count, rel_tol=0, abs_tol=1e-6):
        raise ValueError(f"{label} must be a positive integer multiple of dt")
    return count


def next_power_of_two(n: int) -> int:
    if n < 1:
        raise ValueError("FFT size must be positive")
    return 1 << (n - 1).bit_length()


def read_s21(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read standard Touchstone 1.x two-port RI, MA, or DB data (S21 column)."""
    options = "GHZ S MA R 50".split()
    values = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            options = line[1:].upper().split()
        elif line.startswith("["):
            raise ValueError("This reader supports Touchstone 1.x .s2p files")
        else:
            values.extend(float(v) for v in line.split())
    if len(options) < 3 or options[1] != "S" or len(values) % 9:
        raise ValueError("Expected two-port S parameter data with 9 values per point")
    units = {"HZ": 1, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}
    if options[0] not in units:
        raise ValueError("Unrecognized Touchstone frequency unit")
    data = np.asarray(values, dtype=float).reshape(-1, 9)
    freq = data[:, 0] * units[options[0]]
    first, second = data[:, 3], data[:, 4]
    if options[2] == "RI":
        s21 = first + 1j * second
    elif options[2] == "MA":
        s21 = first * np.exp(1j * np.deg2rad(second))
    elif options[2] == "DB":
        s21 = 10 ** (first / 20) * np.exp(1j * np.deg2rad(second))
    else:
        raise ValueError("Unrecognized Touchstone complex data format")
    if len(freq) < 2 or freq[0] != 0 or np.any(np.diff(freq) <= 0):
        raise ValueError("Frequency data must start at DC and increase")
    df = freq[1] - freq[0]
    if not np.allclose(np.diff(freq), df, rtol=1e-9, atol=1e-3):
        raise ValueError("Frequency samples must be uniformly spaced")
    return freq, s21


def sampled_impulse(
    freq: np.ndarray, s21: np.ndarray, dt_s: float,
    tail_s: float, alpha: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Convert Hermitian frequency data to a dimensionless discrete FIR.

    The continuous impulse is h_cont = h_discrete / dt. Therefore ordinary
    discrete convolution of voltage samples with h_discrete already includes
    the integration dt; adding another dt would give an incorrect voltage.
    irfft supplies conjugate negative frequencies and handles DC exactly once.
    This equals the lecture's twice-real one-sided IFFT with half-weight DC.
    Unmeasured frequencies above fmax are set to zero, not extrapolated.
    """
    if alpha < 0:
        raise ValueError("Kaiser alpha cannot be negative")
    df = float(freq[1] - freq[0])
    period = 1 / df
    n_period = integer_points(period, dt_s, "1/df")
    if n_period % 2:
        raise ValueError("Use an even point count for the sampled impulse period")
    if freq[-1] >= 1 / (2 * dt_s):
        raise ValueError("Sampling Nyquist frequency must exceed measured fmax")
    tail_points = integer_points(tail_s, dt_s, "tail duration")
    if tail_points >= n_period:
        raise ValueError("Tail duration must be smaller than 1/df")
    window = np.i0(np.pi * alpha * np.sqrt(np.maximum(0, 1 - (freq / freq[-1]) ** 2)))
    window /= np.i0(np.pi * alpha)
    positive_spectrum = np.zeros(n_period // 2 + 1, dtype=complex)
    positive_spectrum[: len(s21)] = window * s21
    full_h = np.fft.irfft(positive_spectrum, n=n_period)
    h = full_h[: tail_points + 1].copy()
    peak = float(np.max(np.abs(full_h)))
    # The last part of a period wraps the band-limited negative-time response.
    # Evaluate the positive-time tail only, ending before the final 20%.
    tail_end = max(tail_points + 1, int(0.8 * n_period))
    positive_tail = full_h[tail_points + 1 : tail_end]
    tail_peak = float(np.max(np.abs(positive_tail))) if positive_tail.size else None
    metadata = {
        "frequency_step_hz": df, "frequency_max_hz": float(freq[-1]),
        "frequency_samples": len(freq), "impulse_period_s": period,
        "impulse_period_points": n_period, "kaiser_alpha": alpha,
        "kaiser_beta": float(np.pi * alpha), "fir_points": len(h),
        "tail_duration_s": tail_s,
        "peak_delay_s": float(np.argmax(np.abs(full_h)) * dt_s),
        "peak_continuous_impulse_per_s": peak / dt_s,
        "channel_dc_gain": float(s21[0].real),
        "retained_discrete_dc_gain": float(h.sum()),
        "dc_gain_difference": float(h.sum() - s21[0].real),
        "positive_tail_check_start_s": (tail_points + 1) * dt_s,
        "positive_tail_check_end_s": tail_end * dt_s,
        "positive_tail_peak_relative_db":
            float(20 * np.log10(tail_peak / peak)) if tail_peak else None,
        "tail_target_db": -80.0,
        "tail_target_met_on_checked_interval":
            bool(tail_peak is not None and tail_peak / peak < 1e-4),
        "scaling": "h_discrete = dt*h_continuous; no extra dt in discrete convolution",
        "truncation_note": "No DC renormalization; causal truncation discards the band-limited periodic negative-time wrap. Tail bound applies only to the stated positive-time interval.",
    }
    return h, full_h, metadata


def read_bits(path: Path) -> np.ndarray:
    raw = b"".join(path.read_bytes().split())
    values = np.frombuffer(raw, dtype=np.uint8)
    if values.size == 0 or np.any((values != ord("0")) & (values != ord("1"))):
        raise ValueError("Bit pattern must contain only 0, 1, and whitespace")
    return (values == ord("1")).copy()


def input_section(
    bits: np.ndarray, samples_per_ui: int, rise_points: int,
    previous_bit: bool = False, amplitude_v: float = 1.0,
) -> np.ndarray:
    """Create incident NRZ samples, including the true preceding block bit.

    Transition starts at the UI boundary. At sample j of the transition,
    v = previous_level + (new_level-previous_level)*j/rise_points.
    The exact new level is reached at sample rise_points (20 ps by default).
    """
    if samples_per_ui < 1 or not 1 <= rise_points <= samples_per_ui:
        raise ValueError("Rise point count must be between 1 and samples_per_ui")
    bits = np.asarray(bits, dtype=bool)
    result = np.repeat(bits.astype(np.float32) * amplitude_v, samples_per_ui)
    if not bits.size:
        return result
    previous = np.concatenate(([bool(previous_bit)], bits[:-1]))
    changed = np.flatnonzero(previous != bits)
    ramp = np.arange(rise_points, dtype=np.float32) / rise_points
    for index in changed:
        old = float(previous[index]) * amplitude_v
        new = float(bits[index]) * amplitude_v
        start = index * samples_per_ui
        result[start : start + rise_points] = old + (new - old) * ramp
    return result


def simulate(
    bits: np.ndarray, h: np.ndarray, samples_per_ui: int,
    rise_points: int, block_bits: int = 512,
    output_path: Path | None = None, progress: bool = False,
) -> tuple[np.ndarray, dict]:
    """Linear convolution by FFT overlap-add, including any final partial block.

    FFTs and individual block outputs use float64; the accumulated waveform
    is stored in float32 after each addition.
    A disk-backed output keeps resident memory bounded by the block and FIR.
    No full-length time or input waveform is allocated.
    """
    if len(bits) == 0 or len(h) == 0 or block_bits < 1:
        raise ValueError("Bits, FIR, and block size must be nonempty/positive")
    total_points = len(bits) * samples_per_ui + len(h) - 1
    if output_path is None:
        y = np.zeros(total_points, dtype=np.float32)
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        y = np.memmap(output_path, dtype=np.float32, mode="w+", shape=(total_points,))
        y[:] = 0
    fft_points = next_power_of_two(block_bits * samples_per_ui + len(h) - 1)
    hf = np.fft.rfft(h, n=fft_points)
    n_blocks = math.ceil(len(bits) / block_bits)
    tic = time.perf_counter()
    for k, start_bit in enumerate(range(0, len(bits), block_bits)):
        end_bit = min(start_bit + block_bits, len(bits))
        x = input_section(
            bits[start_bit:end_bit], samples_per_ui, rise_points,
            bool(bits[start_bit - 1]) if start_bit else False,
        )
        count = len(x) + len(h) - 1
        section = np.fft.irfft(np.fft.rfft(x, n=fft_points) * hf, n=fft_points)[:count]
        start = start_bit * samples_per_ui
        y[start : start + count] += section
        if progress and (k + 1 == n_blocks or (k + 1) % max(1, n_blocks // 10) == 0):
            print(f"Convolution: {k + 1}/{n_blocks} blocks", flush=True)
    if isinstance(y, np.memmap):
        y.flush()
    elapsed = time.perf_counter() - tic
    metadata = {
        "bit_count": len(bits), "processed_bit_count": len(bits),
        "block_bits": block_bits, "block_count": n_blocks,
        "final_block_bits": len(bits) - (n_blocks - 1) * block_bits,
        "samples_per_ui": samples_per_ui, "rise_points": rise_points,
        "fft_points": fft_points,
        "maximum_linear_output_points_per_block": block_bits * samples_per_ui + len(h) - 1,
        "fft_padding_valid": bool(fft_points >= block_bits * samples_per_ui + len(h) - 1),
        "output_points": total_points, "output_dtype": str(y.dtype),
        "output_bytes": int(y.nbytes), "convolution_seconds": elapsed,
    }
    return y, metadata


def direct_convolution_check() -> dict:
    """An independent whole-waveform convolution checks the small-block solver."""
    bits = np.array([1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1], dtype=bool)
    h = np.array([0.0, 0.25, 0.5, 0.2, 0.05])
    x = input_section(bits, 10, 2)
    expected = np.convolve(x.astype(float), h)
    actual, metadata = simulate(bits, h, 10, 2, block_bits=4)
    error = float(np.max(np.abs(expected - actual)))
    return {"bit_count": len(bits), "block_bits": 4,
            "final_block_bits": metadata["final_block_bits"],
            "max_absolute_error_v": error, "tolerance_v": 2e-7,
            "passed": bool(error < 2e-7)}


def compare_ads(
    path: Path, y: np.ndarray, dt_s: float, windows: list[tuple[float, float]],
    chunk_rows: int = 200_000, valid_end_s: float | None = None,
) -> tuple[dict, list[dict]]:
    """Stream plain/gzipped ADS data; compare at ADS timestamps by interpolation.

    ADS timestamps are used directly, so its nonuniform first step is retained.
    No data-dependent time shift, voltage offset, or amplitude fit is applied.
    Samples after the supplied bit pattern ends are excluded: ADS continues
    through 105 us, while the known input ends at 104.8576 us. The later input
    bits are unavailable, so a zero-padded synthetic continuation is invalid.
    """
    sums = {"count": 0, "squared": 0.0, "absolute": 0.0, "signed": 0.0,
            "max": 0.0, "first": None, "last": None, "min_dt": math.inf,
            "max_dt": 0.0, "previous_t": None, "file_count": 0,
            "file_first": None, "file_last": None, "excluded_count": 0}
    comparison_end = valid_end_s if valid_end_s is not None else (len(y) - 1) * dt_s
    if comparison_end > (len(y) - 1) * dt_s + 1e-15:
        raise ValueError("Comparison interval exceeds computed output")
    slices = [{"start_s": lo, "end_s": hi, "time_s": [], "ads_v": [], "computed_v": []}
              for lo, hi in windows]
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        header = next(handle)
        if not header.lower().strip().startswith("time"):
            raise ValueError("Expected ADS header: time vout")
        while True:
            rows = []
            for _ in range(chunk_rows):
                line = handle.readline()
                if not line:
                    break
                rows.append(line)
            if not rows:
                break
            block = np.loadtxt(rows, dtype=float, ndmin=2)
            if block.shape[1] != 2 or not np.all(np.isfinite(block)):
                raise ValueError("ADS rows must contain finite time and voltage")
            t, reference = block[:, 0], block[:, 1]
            diffs = np.diff(t)
            if sums["previous_t"] is not None:
                diffs = np.concatenate(([t[0] - sums["previous_t"]], diffs))
            if np.any(diffs <= 0):
                raise ValueError("ADS times must increase strictly")
            if diffs.size:
                sums["min_dt"] = min(sums["min_dt"], float(diffs.min()))
                sums["max_dt"] = max(sums["max_dt"], float(diffs.max()))
            sums["previous_t"] = float(t[-1])
            sums["file_count"] += len(t)
            if sums["file_first"] is None:
                sums["file_first"] = float(t[0])
            sums["file_last"] = float(t[-1])
            if t[0] < 0:
                raise ValueError("ADS timestamps must be nonnegative")
            valid = t < comparison_end
            sums["excluded_count"] += int(np.sum(~valid))
            t, reference = t[valid], reference[valid]
            if len(t) == 0:
                continue
            index = t / dt_s
            lower = np.minimum(np.floor(index).astype(np.int64), len(y) - 2)
            weight = index - lower
            predicted = y[lower].astype(float) * (1 - weight) + y[lower + 1] * weight
            error = predicted - reference
            sums["count"] += len(t)
            sums["squared"] += float(np.dot(error, error))
            sums["absolute"] += float(np.abs(error).sum())
            sums["signed"] += float(error.sum())
            sums["max"] = max(sums["max"], float(np.abs(error).max()))
            if sums["first"] is None:
                sums["first"] = float(t[0])
            sums["last"] = float(t[-1])
            for window in slices:
                selection = (t >= window["start_s"]) & (t < window["end_s"])
                if np.any(selection):
                    for name, array in (("time_s", t), ("ads_v", reference), ("computed_v", predicted)):
                        window[name].append(array[selection].copy())
    if sums["count"] == 0:
        raise ValueError("ADS file has no samples")
    count = sums["count"]
    overall = {"sample_count": count, "start_s": sums["first"], "end_s": sums["last"],
               "file_sample_count": sums["file_count"],
               "file_start_s": sums["file_first"], "file_end_s": sums["file_last"],
               "excluded_sample_count": sums["excluded_count"],
               "known_input_end_s": comparison_end,
               "exclusion_reason": "ADS extends beyond supplied bits; compare only 0 <= t < known input duration",
               "rmse_v": math.sqrt(sums["squared"] / count),
               "mean_absolute_error_v": sums["absolute"] / count,
               "mean_signed_error_v": sums["signed"] / count,
               "max_absolute_error_v": sums["max"],
               "min_ads_step_s": sums["min_dt"], "max_ads_step_s": sums["max_dt"],
               "alignment": "Direct ADS timestamps, linear interpolation only; no fitted shift, gain, or offset"}
    for window in slices:
        for name in ("time_s", "ads_v", "computed_v"):
            window[name] = np.concatenate(window[name]) if window[name] else np.array([])
        error = window["computed_v"] - window["ads_v"]
        if error.size == 0:
            raise ValueError("A requested comparison window has no ADS samples")
        window["sample_count"] = len(error)
        window["rmse_v"] = float(np.sqrt(np.mean(error ** 2)))
        window["max_absolute_error_v"] = float(np.max(np.abs(error)))
        window["mean_signed_error_v"] = float(np.mean(error))
    return overall, slices


def figures(output: Path, full_h: np.ndarray, h: np.ndarray, dt_s: float,
            bits: np.ndarray, samples_per_ui: int, rise_points: int, windows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "figure.dpi": 140, "axes.spines.top": False,
                         "axes.spines.right": False})
    t = np.arange(len(full_h)) * dt_s
    normalized_db = 20 * np.log10(np.maximum(np.abs(full_h) / np.max(np.abs(full_h)), 1e-12))
    fig, ax = plt.subplots(figsize=(9, 3.8), layout="constrained")
    ax.plot(t * 1e9, normalized_db, linewidth=0.65)
    ax.axvline((len(h) - 1) * dt_s * 1e9, color="#c65c12", linestyle="--", label="Retained tail: 40 ns" if len(h)*dt_s<40.01e-9 and len(h)*dt_s>39.99e-9 else "Tail cutoff")
    ax.axhline(-80, color="#6b6b6b", linestyle=":", label="-80 dB relative to peak")
    ax.set(xlim=(0, t[-1] * 1e9), ylim=(-140, 5), xlabel="Time (ns)",
           ylabel="Impulse magnitude / peak (dB)", title="Band-limited impulse response over one 1/df period")
    ax.grid(alpha=0.2); ax.legend(loc="upper right")
    fig.savefig(output / "impulse_response.png"); plt.close(fig)
    x = input_section(bits[:20], samples_per_ui, rise_points)
    tx = np.arange(len(x)) * dt_s
    fig, ax = plt.subplots(figsize=(9, 3.3), layout="constrained")
    ax.plot(tx * 1e9, x, color="#16718c", linewidth=1.3)
    ax.set(xlabel="Time (ns)", ylabel="Incident voltage (V)",
           title="First 20 bits: 10 Gbps NRZ with 20 ps rise / fall time", ylim=(-0.08, 1.08))
    ax.grid(alpha=0.2); fig.savefig(output / "input_waveform.png"); plt.close(fig)
    for i, window in enumerate(windows, 1):
        ts = (window["time_s"] - window["start_s"]) * 1e9
        fig, (ax, error_ax) = plt.subplots(2, 1, figsize=(9, 4.4), sharex=True,
                                           gridspec_kw={"height_ratios": [3, 1]}, layout="constrained")
        ax.plot(ts, window["computed_v"], color="#116b8c", linewidth=1.1, label="Overlap-add")
        ax.plot(ts, window["ads_v"], color="#d15a2a", linewidth=0.9, linestyle="--", label="ADS")
        ax.set(ylabel="Voltage (V)", title=f"10 ns window starting at {window['start_s']*1e6:.4f} us")
        ax.legend(loc="upper right", ncol=2); ax.grid(alpha=0.2)
        error_ax.plot(ts, (window["computed_v"] - window["ads_v"]) * 1e3, linewidth=0.8, color="#5e4c8d")
        error_ax.set(xlabel="Time from window start (ns)", ylabel="Error (mV)", xlim=(0, 10))
        error_ax.grid(alpha=0.2)
        fig.savefig(output / f"comparison_{i}.png"); plt.close(fig)


def verification_study(
    bits: np.ndarray, freq: np.ndarray, s21: np.ndarray,
    windows: list[dict], base_dt: float, base_tail: float, base_alpha: float,
) -> dict:
    """Refine dt, window strength, and T2 on the exact same four ADS windows.

    Each window is recomputed with its full preceding FIR history, taken from
    the supplied bits. Local calculation equals the whole-sequence finite-FIR
    result on the window. Nothing is fitted to the reference.
    """
    variants = [
        ("half_sample_interval", base_dt / 2, base_tail, base_alpha),
        ("no_kaiser_taper", base_dt, base_tail, 0.0),
        ("weak_kaiser_taper", base_dt, base_tail, 0.1),
        ("shorter_tail", base_dt, base_tail / 2, base_alpha),
        ("longer_tail", base_dt, min(base_tail * 1.5, .7 / (freq[1] - freq[0])), base_alpha),
    ]
    results = []
    ui, rise = 100e-12, 20e-12
    for label, dt, tail, alpha in variants:
        n_ui = integer_points(ui, dt, "UI")
        n_rise = integer_points(rise, dt, "rise")
        h, _, impulse_info = sampled_impulse(freq, s21, dt, tail, alpha)
        errors, differences, per_window = [], [], []
        tic = time.perf_counter()
        for window in windows:
            # Include one extra UI before the full finite impulse history.
            first_bit = max(0, math.floor((window["start_s"] - tail) / ui) - 1)
            last_bit = min(len(bits), math.ceil(window["end_s"] / ui) + 1)
            # simulate assumes prior level zero, so explicitly generate each
            # local waveform with its true prior bit before direct FFT OLA.
            local_x = input_section(bits[first_bit:last_bit], n_ui, n_rise,
                                    bool(bits[first_bit - 1]) if first_bit else False)
            fft_n = next_power_of_two(len(local_x) + len(h) - 1)
            local_y = np.fft.irfft(np.fft.rfft(local_x, fft_n) * np.fft.rfft(h, fft_n), fft_n)
            index = (window["time_s"] - first_bit * ui) / dt
            lower = np.floor(index).astype(np.int64)
            fraction = index - lower
            predicted = local_y[lower] * (1 - fraction) + local_y[lower + 1] * fraction
            error = predicted - window["ads_v"]
            difference = predicted - window["computed_v"]
            errors.append(error); differences.append(difference)
            per_window.append({"start_s": window["start_s"],
                               "rmse_v": float(np.sqrt(np.mean(error ** 2))),
                               "max_absolute_error_v": float(np.max(np.abs(error)))})
        error = np.concatenate(errors); difference = np.concatenate(differences)
        results.append({"variant": label, "dt_s": dt, "tail_s": tail, "alpha": alpha,
                        "ads_rmse_v": float(np.sqrt(np.mean(error ** 2))),
                        "ads_max_absolute_error_v": float(np.max(np.abs(error))),
                        "rms_change_from_baseline_v": float(np.sqrt(np.mean(difference ** 2))),
                        "max_change_from_baseline_v": float(np.max(np.abs(difference))),
                        "retained_dc_gain": impulse_info["retained_discrete_dc_gain"],
                        "compute_seconds": time.perf_counter() - tic,
                        "windows": per_window})
    baseline_error = np.concatenate([w["computed_v"] - w["ads_v"] for w in windows])
    return {"scope": "The four 10 ns windows, using all finite-FIR preceding input history",
            "reference_samples": len(baseline_error),
            "baseline_ads_rmse_v": float(np.sqrt(np.mean(baseline_error ** 2))),
            "baseline_ads_max_absolute_error_v": float(np.max(np.abs(baseline_error))),
            "no_fitted_parameters": True, "variants": results,
            "interpretation": "Use the changes to assess numerical sensitivity; these window checks do not establish a unique cause of the residual ADS disagreement."}


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve().parent
    parser.add_argument("--data-dir", type=Path, default=here / "instructions")
    parser.add_argument("--output-dir", type=Path, default=here / "results")
    parser.add_argument("--samples-per-ui", type=int, default=50, help="Default 50 => dt=2 ps")
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--tail-ns", type=float, default=40)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--keep-waveform", action="store_true", help="Keep the large float32 output file")
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--verification-study", action="store_true",
                        help="Compare half dt, alternate tapers, and tail lengths on the four windows")
    args = parser.parse_args(argv)
    if args.samples_per_ui < 1:
        parser.error("--samples-per-ui must be positive")
    start = time.perf_counter()
    output = args.output_dir.resolve(); output.mkdir(parents=True, exist_ok=True)
    data = args.data_dir.resolve()
    bits_path = data / "HW3_bitpattern_1048576.txt"
    s2p_path = data / "HW3_mline_channel.s2p"
    ads_path = data / "HW3_Rx_waveform_ADS.txt"
    if not ads_path.exists():
        ads_path = data / "HW3_Rx_waveform_ADS.txt.gz"
    bits = read_bits(bits_path)
    bit_rate, ui, rise = 10e9, 100e-12, 20e-12
    dt = ui / args.samples_per_ui
    rise_points = integer_points(rise, dt, "rise time")
    freq, s21 = read_s21(s2p_path)
    h, full_h, impulse_metadata = sampled_impulse(freq, s21, dt, args.tail_ns * 1e-9, args.alpha)
    check = direct_convolution_check()
    if not check["passed"]:
        raise RuntimeError("Direct-convolution validation failed")
    scratch = output / "waveform.float32"
    y, simulation_metadata = simulate(bits, h, args.samples_per_ui, rise_points,
                                      args.block_bits, scratch, progress=True)
    duration = len(bits) * ui
    starts = [0.0, 8e-6, 50e-6, duration - 10e-9]
    if duration < 50.01e-6:
        starts = list(np.linspace(0, max(0, duration - 10e-9), 4))
    windows = [(float(lo), float(lo + 10e-9)) for lo in starts]
    ads_tic = time.perf_counter()
    overall, slices = compare_ads(ads_path, y, dt, windows, valid_end_s=duration)
    ads_seconds = time.perf_counter() - ads_tic
    study = verification_study(bits, freq, s21, slices, dt, args.tail_ns * 1e-9, args.alpha) if args.verification_study else None
    if not args.no_figures:
        figures(output, full_h, h, dt, bits, args.samples_per_ui, rise_points, slices)
    serial_windows = [{k: v for k, v in window.items() if k not in ("time_s", "ads_v", "computed_v")}
                      for window in slices]
    # Compact reusable ADS/solver points, only the requested four short windows.
    for index, window in enumerate(slices, 1):
        np.savetxt(output / f"comparison_{index}.csv",
                   np.column_stack([window["time_s"], window["ads_v"], window["computed_v"]]),
                   delimiter=",", header="time_s,ads_v,overlap_add_v", comments="", fmt="%.12e")
    summary = {
        "course": "NCU EE8091", "homework": "HW3 sectional convolution",
        "method": "Finite causal FIR from Hermitian S21 followed by FFT overlap-add",
        "parameters": {"bit_rate_bps": bit_rate, "ui_s": ui, "dt_s": dt,
                       "rise_fall_s": rise, "incident_low_v": 0.0, "incident_high_v": 1.0,
                       "line_length_m": 0.2, "source_ohm": 50, "load_ohm": 50},
        "input": {"bits": len(bits), "ones": int(bits.sum()),
                  "zeros": int(len(bits) - bits.sum()), "duration_s": duration,
                  "last_20_bits": "".join("1" if b else "0" for b in bits[-20:])},
        "impulse": impulse_metadata, "simulation": simulation_metadata,
        "direct_convolution_check": check, "ads_comparison": overall,
        "comparison_windows": serial_windows,
        "verification_study": study,
        "timing": {"ads_streaming_seconds": ads_seconds,
                   "total_seconds": time.perf_counter() - start},
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "platform": platform.system(), "processor": platform.machine()},
        "files": {"bits": bits_path.name, "channel": s2p_path.name, "ads": ads_path.name,
                  "waveform_saved": bool(args.keep_waveform)},
        "limitations": [
            "S21 is measured only through 40 GHz; higher frequencies are zero-filled and a Kaiser window is applied.",
            "The finite causal FIR omits the periodic negative-time wrap and the tail after T2; the DC gain difference is reported without correction.",
            "ADS comparison includes finite input, sampling, windowing, truncation, and simulator-model differences; it is not a proof of zero error.",
            "ADS continues to 105 us, beyond the supplied 104.8576 us bit pattern. Samples outside the known-input interval are excluded, not extrapolated.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    del y
    if not args.keep_waveform:
        scratch.unlink()
    print(json.dumps({"bits": len(bits), "convolution_seconds": simulation_metadata["convolution_seconds"],
                      "ads_rmse_v": overall["rmse_v"], "ads_max_absolute_error_v": overall["max_absolute_error_v"],
                      "output": str(output)}, indent=2))
    return summary


if __name__ == "__main__":
    main()
