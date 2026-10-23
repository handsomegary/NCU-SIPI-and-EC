"""Reproducible, anonymous HW4 analysis of finite-rise-time NRZ signals.

Run ``python src/hw4.py`` from the homework directory. Paths are relative to
this file, so the same command also works from another working directory.
Only the supplied bit streams are read; all generated files go to results/.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from numbers import Integral
from pathlib import Path
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import rfft, rfftfreq


@dataclass(frozen=True)
class Config:
    bit_duration: float = 100e-12
    rise_time: float = 20e-12
    dt: float = 2e-12
    vpp: float = 1.0
    z0: float = 50.0
    block_bits: int = 2000
    nfft: int = 2**19
    kaiser_beta: float = float(np.pi)
    averages: tuple[int, ...] = (10, 100)

    @property
    def samples_per_bit(self) -> int:
        return round(self.bit_duration / self.dt)

    @property
    def rise_samples(self) -> int:
        return round(self.rise_time / self.dt)

    @property
    def block_samples(self) -> int:
        return self.block_bits * self.samples_per_bit

    def validate(self) -> None:
        positive = (self.bit_duration, self.rise_time, self.dt, self.vpp, self.z0)
        if not all(np.isfinite(x) and x > 0 for x in positive):
            raise ValueError("Time, voltage, and impedance parameters must be finite and positive.")
        if self.rise_time > self.bit_duration:
            raise ValueError("Rise time must not exceed one bit duration.")
        for name, value in (("bit duration", self.bit_duration), ("rise time", self.rise_time)):
            ratio = value / self.dt
            if round(ratio) < 1 or not np.isclose(ratio, round(ratio), rtol=0, atol=1e-9):
                raise ValueError(f"The {name} must be an integer multiple of the sample interval.")
        integers = (self.block_bits, self.nfft, *self.averages)
        if not all(isinstance(value, Integral) and not isinstance(value, bool) for value in integers):
            raise ValueError("Block length, FFT length, and averaging counts must be integers.")
        if self.block_bits < 1 or self.nfft < self.block_samples:
            raise ValueError("Block length must be positive and nfft must be at least the block sample count.")
        if not np.isfinite(self.kaiser_beta) or self.kaiser_beta < 0:
            raise ValueError("Kaiser beta must be finite and nonnegative.")
        if not self.averages or any(k < 1 for k in self.averages):
            raise ValueError("At least one positive averaging count is required.")
        if len(set(self.averages)) != len(self.averages):
            raise ValueError("Averaging counts must be distinct.")


def load_bits(path: Path | str) -> np.ndarray:
    """Read ASCII zero/one data, allowing whitespace but no other characters."""
    raw = Path(path).read_text(encoding="utf-8-sig")
    compact = "".join(raw.split())
    if not compact or set(compact) - {"0", "1"}:
        raise ValueError(f"{Path(path).name} must contain only binary digits and whitespace.")
    return np.frombuffer(compact.encode("ascii"), dtype=np.uint8).copy() - ord("0")


def build_nrz_waveform(
    bits: np.ndarray,
    previous_bit: int = 0,
    samples_per_bit: int = 50,
    rise_samples: int = 10,
    vpp: float = 1.0,
) -> np.ndarray:
    """Construct bipolar NRZ, with a linear ramp at each transition.

    The sample at a bit boundary retains the preceding level; the new level
    is reached after ``rise_samples`` sample intervals. Pass the actual
    preceding bit when constructing successive blocks to keep transitions
    continuous. The first bit of the complete record assumes preceding zero.
    """
    bits = np.asarray(bits)
    if bits.ndim != 1 or bits.size == 0 or not np.all((bits == 0) | (bits == 1)):
        raise ValueError("bits must be a nonempty one-dimensional binary array.")
    if not all(isinstance(value, Integral) and not isinstance(value, bool) for value in (samples_per_bit, rise_samples)):
        raise ValueError("Bit and ramp sample counts must be integers.")
    if previous_bit not in (0, 1) or samples_per_bit < 1 or not 0 <= rise_samples <= samples_per_bit:
        raise ValueError("Invalid predecessor or ramp sample count.")
    if not np.isfinite(vpp) or vpp <= 0:
        raise ValueError("vpp must be finite and positive.")
    levels = (bits.astype(np.float64) - 0.5) * vpp
    waveform = np.repeat(levels, samples_per_bit)
    if rise_samples == 0:
        return waveform
    preceding = np.concatenate(([float(previous_bit) * vpp - vpp / 2], levels[:-1]))
    transitions = np.flatnonzero(levels != preceding)
    offsets = np.arange(rise_samples)
    ramps = preceding[transitions, None] + (
        levels[transitions] - preceding[transitions]
    )[:, None] * offsets[None, :] / rise_samples
    waveform[transitions[:, None] * samples_per_bit + offsets[None, :]] = ramps
    return waveform


def kaiser_window(length: int, beta: float = float(np.pi)) -> np.ndarray:
    """Symmetric Kaiser window; beta = pi * alpha, so beta=pi means alpha=1."""
    if not isinstance(length, Integral) or isinstance(length, bool) or length < 1 or not np.isfinite(beta) or beta < 0:
        raise ValueError("Window length must be positive and beta nonnegative.")
    if length == 1:
        return np.ones(1)
    coordinate = 2 * np.arange(length, dtype=np.float64) / (length - 1) - 1
    return np.i0(beta * np.sqrt(np.maximum(0.0, 1 - coordinate**2))) / np.i0(beta)


def one_sided_periodogram(
    waveform: np.ndarray,
    dt: float,
    z0: float = 50.0,
    nfft: int | None = None,
    window: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return frequency (Hz) and one-sided PSD (W/Hz), with no detrending.

    With U = mean(window**2), S2 = dt**2 |FFT(x*w)|**2 / (T*z0*U).
    Interior positive-frequency bins double; DC and an even FFT's Nyquist
    bin do not. Summing S1*df therefore exactly recovers the window-weighted
    mean-square voltage divided by impedance (including DC).
    """
    waveform = np.asarray(waveform, dtype=np.float64)
    if waveform.ndim != 1 or waveform.size == 0 or not np.all(np.isfinite(waveform)):
        raise ValueError("waveform must be a nonempty, finite one-dimensional array.")
    if not np.isfinite(dt) or dt <= 0 or not np.isfinite(z0) or z0 <= 0:
        raise ValueError("dt and z0 must be finite and positive.")
    nfft = waveform.size if nfft is None else nfft
    if not isinstance(nfft, Integral) or isinstance(nfft, bool) or nfft < waveform.size:
        raise ValueError("nfft must not truncate the waveform.")
    window = np.ones(waveform.size) if window is None else np.asarray(window, dtype=np.float64)
    if window.shape != waveform.shape or not np.all(np.isfinite(window)) or np.sum(window**2) <= 0:
        raise ValueError("window must match the waveform and have positive finite energy.")
    spectrum = rfft(waveform * window, n=nfft)
    psd = np.abs(spectrum) ** 2 * dt / (z0 * np.sum(window**2))
    psd[1:-1 if nfft % 2 == 0 else None] *= 2
    return rfftfreq(nfft, d=dt), psd


def theoretical_random_nrz_psd(
    frequency_hz: np.ndarray,
    bit_duration: float = 100e-12,
    rise_time: float = 20e-12,
    vpp: float = 1.0,
    z0: float = 50.0,
) -> np.ndarray:
    """One-sided continuous-time PSD for independent, equiprobable NRZ bits.

    The finite linear rise is modeled as rectangular smoothing of width Tr.
    This is a reference envelope, not a theory for coded or periodic data.
    The DC sample uses the undoubled two-sided value, consistent with rFFT.
    """
    f = np.asarray(frequency_hz, dtype=np.float64)
    if np.any(f < 0) or not np.all(np.isfinite(f)):
        raise ValueError("The one-sided frequency axis must be finite and nonnegative.")
    if not all(np.isfinite(value) and value > 0 for value in (bit_duration, vpp, z0)):
        raise ValueError("Bit duration, voltage swing, and impedance must be finite and positive.")
    if not np.isfinite(rise_time) or not 0 <= rise_time <= bit_duration:
        raise ValueError("Rise time must be finite and between zero and one bit duration.")
    two_sided = vpp**2 * bit_duration / (4 * z0) * np.sinc(f * bit_duration)**2 * np.sinc(f * rise_time)**2
    return two_sided * np.where(f == 0, 1.0, 2.0)


def sampled_random_nrz_psd(frequency_hz: np.ndarray, config: Config) -> np.ndarray:
    """Exact sampled-model reference, including discrete hold and ramp filters.

    Away from the removable zeros this equals the continuous reference times
    sinc(f*dt)**(-4). It allows the small discretization error of the plotted
    continuous-time reference to be quantified without changing the model.
    """
    f = np.asarray(frequency_hz, dtype=np.float64)
    denominator = np.sinc(f * config.dt)
    hold = np.sinc(f * config.bit_duration) / denominator
    ramp = np.sinc(f * config.rise_time) / denominator
    return (config.vpp**2 * config.bit_duration / (4 * config.z0)
            * hold**2 * ramp**2 * np.where(f == 0, 1.0, 2.0))


def analyze_pattern(bits: np.ndarray) -> dict:
    """Summarize observable structure without claiming a specific line code."""
    bits = np.asarray(bits)
    if bits.ndim != 1 or bits.size == 0 or not np.all((bits == 0) | (bits == 1)):
        raise ValueError("bits must be a nonempty one-dimensional binary array.")
    bipolar = 2 * bits.astype(np.int64) - 1
    boundaries = np.flatnonzero(np.diff(bits)) + 1
    runs = np.diff(np.concatenate(([0], boundaries, [bits.size])))
    disparity = np.concatenate(([0], np.cumsum(bipolar)))
    p1 = float(np.mean(bits))
    centered = bits.astype(float) - p1
    variance = float(np.mean(centered**2))
    correlations = {
        str(lag): float(np.mean(centered[:-lag] * centered[lag:]) / variance) if variance else None
        for lag in (1, 2, 3, 4, 5, 8, 10, 16, 20, 32, 64, 100) if lag < bits.size
    }
    width = min(10, bits.size)
    period10_ones = [float(np.mean(bits[offset::10])) for offset in range(width)]
    full_words = bits[:bits.size // 10 * 10].reshape(-1, 10)
    word_disparities = 2 * np.sum(full_words, axis=1, dtype=np.int64) - 10
    disparity_values, disparity_counts = np.unique(word_disparities, return_counts=True)
    pair_values = 2 * full_words[:, 8] + full_words[:, 9] if full_words.size else np.array([], dtype=int)
    return {
        "total_bits": int(bits.size),
        "ones": int(bits.sum()),
        "ones_fraction": p1,
        "transition_fraction": float(boundaries.size / (bits.size - 1)) if bits.size > 1 else 0.0,
        "maximum_run_length": int(np.max(runs)),
        "run_length_counts": {str(n): int(np.sum(runs == n)) for n in np.unique(runs)},
        "cumulative_disparity_min": int(np.min(disparity)),
        "cumulative_disparity_max": int(np.max(disparity)),
        "cumulative_disparity_final": int(disparity[-1]),
        "bit_autocorrelation": correlations,
        "ones_fraction_by_position_modulo_10": period10_ones,
        "aligned_10bit_word_disparity_counts": {
            str(int(value)): int(count) for value, count in zip(disparity_values, disparity_counts)
        },
        "aligned_10bit_last_pair_counts": {
            format(value, "02b"): int(np.sum(pair_values == value)) for value in range(4)
        },
    }


def analyze_file(path: Path | str, config: Config) -> dict:
    """Compute requested K checkpoints from consecutive, nonoverlapping blocks."""
    config.validate()
    bits = load_bits(path)
    needed_bits = max(config.averages) * config.block_bits
    if bits.size < needed_bits:
        raise ValueError(f"{Path(path).name}: need {needed_bits} bits, found {bits.size}.")
    window = kaiser_window(config.block_samples, config.kaiser_beta)
    window_energy = float(np.sum(window**2))
    psd_sum = np.zeros(config.nfft // 2 + 1)
    frequency = rfftfreq(config.nfft, d=config.dt)
    sampled_reference = sampled_random_nrz_psd(frequency, config)
    reference_band = (frequency >= 0.1e9) & (frequency <= 9e9)
    raw_power_sum = weighted_power_sum = 0.0
    estimates, checkpoints = {}, {}
    for block in range(max(config.averages)):
        start = block * config.block_bits
        waveform = build_nrz_waveform(
            bits[start:start + config.block_bits], int(bits[start - 1]) if start else 0,
            config.samples_per_bit, config.rise_samples, config.vpp,
        )
        frequency, psd = one_sided_periodogram(waveform, config.dt, config.z0, config.nfft, window)
        psd_sum += psd
        raw_power_sum += float(np.mean(waveform**2) / config.z0)
        weighted_power_sum += float(np.sum((waveform * window)**2) / window_energy / config.z0)
        k = block + 1
        if k in config.averages:
            averaged = psd_sum / k
            estimates[k] = averaged.copy()
            frequency_step = 1 / (config.nfft * config.dt)
            integrated = float(np.sum(averaged) * frequency_step)
            weighted_power = weighted_power_sum / k
            checkpoints[str(k)] = {
                "averages": k,
                "bits_used": k * config.block_bits,
                "record_duration_s": k * config.block_bits * config.bit_duration,
                "mean_time_domain_power_w": raw_power_sum / k,
                "mean_window_weighted_power_w": weighted_power,
                "integrated_full_band_psd_w": integrated,
                "parseval_relative_error": abs(integrated - weighted_power) / weighted_power if weighted_power else 0.0,
                "power_0_to_50_GHz_w": float(np.sum(averaged[frequency <= 50e9]) * frequency_step),
                "dc_bin_psd_w_per_hz": float(averaged[0]),
            }
            if np.any(reference_band):
                ratio = averaged[reference_band] / sampled_reference[reference_band]
                checkpoints[str(k)]["sampled_random_reference_comparison"] = {
                    "frequency_band_hz": [0.1e9, 9e9],
                    "mean_psd_ratio": float(np.mean(ratio)),
                    "std_psd_ratio": float(np.std(ratio)),
                    "median_ratio_db": float(np.median(10 * np.log10(np.maximum(ratio, 1e-300)))),
                    "std_ratio_db": float(np.std(10 * np.log10(np.maximum(ratio, 1e-300)))),
                    "interpretation": "Descriptive only; agreement is expected only for the independent random pattern.",
                }
    return {"frequency_hz": frequency, "estimates": estimates, "statistics": analyze_pattern(bits), "checkpoints": checkpoints}


def to_dbm_per_hz(psd_w_hz: np.ndarray) -> np.ndarray:
    """Convert W/Hz to dBm/Hz; the tiny floor is for plotting only."""
    return 10 * np.log10(np.maximum(np.asarray(psd_w_hz) * 1e3, 1e-300))


def _spectrum_axes(ax, frequency, estimate, theory, pattern, k, config):
    mask = frequency <= 50e9
    ax.plot(frequency[mask] / 1e9, to_dbm_per_hz(estimate[mask]), lw=0.8, color="#1558a6", label="Averaged periodogram")
    ax.plot(frequency[mask] / 1e9, to_dbm_per_hz(theory[mask]), lw=1.4, color="#c7473d", label="Independent random NRZ reference")
    ax.set(xlim=(0, 50), ylim=(-155, -55), xlabel="Frequency (GHz)", ylabel="PSD (dBm/Hz)",
           title=f"Pattern {pattern} | K = {k} | T = {config.block_bits * config.bit_duration * 1e9:g} ns")
    ax.grid(alpha=0.28)
    ax.legend(loc="upper right", fontsize=8)


def generate_outputs(data_dir: Path, output_dir: Path, config: Config, write_csv: bool = False) -> dict:
    """Generate the report figures, a JSON audit, and the full-band spectra."""
    config.validate()
    started = time.perf_counter()
    paths = [data_dir / f"HW4_bit_pattern_{index}.txt" for index in (1, 2, 3)]
    missing = [path.name for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing input files: " + ", ".join(missing))
    results = [analyze_file(path, config) for path in paths]
    output_dir.mkdir(parents=True, exist_ok=True)
    frequency = results[0]["frequency_hz"]
    theory = theoretical_random_nrz_psd(frequency, config.bit_duration, config.rise_time, config.vpp, config.z0)
    sampled_theory = sampled_random_nrz_psd(frequency, config)
    window = kaiser_window(config.block_samples, config.kaiser_beta)
    duration = config.block_samples * config.dt
    enbw_hz = float(np.sum(window**2) / np.sum(window)**2 / config.dt)
    artifacts = []
    arrays = {"frequency_hz": frequency, "continuous_random_reference_w_per_hz": theory,
              "sampled_random_reference_w_per_hz": sampled_theory}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "savefig.dpi": 180})
    overview, axes = plt.subplots(3, len(config.averages), figsize=(7 * len(config.averages), 11), squeeze=False)
    patterns = []
    for index, result in enumerate(results, 1):
        pattern_summary = {"pattern": index, "input_file": paths[index - 1].name,
                           "statistics": result["statistics"], "estimates": result["checkpoints"]}
        for column, k in enumerate(config.averages):
            estimate = result["estimates"][k]
            arrays[f"pattern_{index}_K{k}_w_per_hz"] = estimate
            fig, ax = plt.subplots(figsize=(8, 5))
            _spectrum_axes(ax, frequency, estimate, theory, index, k, config)
            fig.tight_layout()
            filename = f"pattern_{index}_K{k}.png"
            fig.savefig(output_dir / filename)
            plt.close(fig)
            artifacts.append(filename)
            _spectrum_axes(axes[index - 1, column], frequency, estimate, theory, index, k, config)
        patterns.append(pattern_summary)
    overview.tight_layout()
    overview.savefig(output_dir / "psd_overview.png")
    plt.close(overview)
    artifacts.append("psd_overview.png")
    fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
    for index, (path, ax) in enumerate(zip(paths, axes), 1):
        preview = build_nrz_waveform(load_bits(path)[:60], 0, config.samples_per_bit, config.rise_samples, config.vpp)
        ax.plot(np.arange(preview.size) * config.dt / 1e-9, preview, color="#1558a6", lw=1)
        ax.set(ylabel="Voltage (V)", title=f"Pattern {index}: first 60 bits", ylim=(-0.65, 0.65))
        ax.grid(alpha=0.25)
    axes[-1].set(xlabel="Time (ns)", xlim=(0, 6))
    fig.tight_layout()
    fig.savefig(output_dir / "waveform_examples.png")
    plt.close(fig)
    artifacts.append("waveform_examples.png")
    np.savez_compressed(output_dir / "spectra.npz", **arrays)
    artifacts.append("spectra.npz")
    if write_csv:
        selected = frequency <= 50e9
        columns = [name for name in arrays if name != "frequency_hz"]
        values = np.column_stack([frequency[selected]] + [arrays[name][selected] for name in columns])
        np.savetxt(output_dir / "spectra_0_to_50_GHz.csv", values, delimiter=",",
                   header=",".join(["frequency_hz"] + columns), comments="", fmt="%.12e")
        artifacts.append("spectra_0_to_50_GHz.csv")
    correction_db = float(-40 * np.log10(np.sinc(50e9 * config.dt)))
    summary = {
        "title": "Spectrum of digital signals: reproducible HW4 analysis",
        "configuration": asdict(config),
        "derived_parameters": {
            "bit_rate_bps": 1 / config.bit_duration,
            "samples_per_bit": config.samples_per_bit,
            "rise_samples": config.rise_samples,
            "block_samples": config.block_samples,
            "block_duration_s": duration,
            "unwindowed_resolution_scale_hz": 1 / duration,
            "fft_grid_step_hz": 1 / (config.nfft * config.dt),
            "nyquist_hz": 1 / (2 * config.dt),
            "window_power_gain_U": float(np.mean(window**2)),
            "window_coherent_gain": float(np.mean(window)),
            "window_enbw_hz": enbw_hz,
            "window_enbw_in_unpadded_bins": enbw_hz * duration,
            "kaiser_alpha": config.kaiser_beta / np.pi,
            "sampled_vs_continuous_reference_at_50_GHz_db": correction_db,
            "independent_random_sampled_expected_total_power_w": float(config.vpp**2 / (4 * config.z0)
                * (1 - (config.rise_time / config.bit_duration) / 3
                   * (1 - 1 / config.rise_samples**2))),
        },
        "method": {
            "segments": "Consecutive nonoverlapping blocks; the same initial blocks are used for every K.",
            "first_predecessor": 0,
            "block_boundaries": "Use the actual preceding input bit; never reset at each block.",
            "detrending": "None; DC and deterministic components are retained.",
            "window": "Symmetric Kaiser window, beta = pi * alpha.",
            "normalization": "S2 = dt**2 * abs(FFT(x*w))**2 / (T * Z0 * U), U = mean(w**2).",
            "one_sided": "Double interior positive-frequency bins only; retain DC and even-FFT Nyquist without doubling.",
            "units": "Stored spectra are W/Hz; figures use 10*log10(PSD*1000) dBm/Hz.",
            "resolution": "Zero padding refines the plotted grid; independent resolution is set by block duration and window.",
            "reference": "The continuous sinc-squared reference assumes independent equiprobable bits, linear ramps, and no sampling.",
            "parseval_check": "Full-band sum(PSD)*df equals window-weighted time-domain power; trapezoidal endpoint weighting is inappropriate for this DFT-bin check.",
        },
        "patterns": patterns,
        "artifacts": artifacts + ["summary.json"],
        "elapsed_seconds": time.perf_counter() - started,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    base = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=base / "instructions")
    parser.add_argument("--output-dir", type=Path, default=base / "results")
    parser.add_argument("--averages", nargs="+", type=int, default=[10, 100], metavar="K")
    parser.add_argument("--block-bits", type=int, default=2000)
    parser.add_argument("--dt-ps", type=float, default=2.0)
    parser.add_argument("--nfft", type=int, default=2**19)
    parser.add_argument("--kaiser-beta", type=float, default=float(np.pi))
    parser.add_argument("--csv", action="store_true", help="Also export the displayed frequency range as CSV.")
    args = parser.parse_args(argv)
    config = Config(dt=args.dt_ps * 1e-12, block_bits=args.block_bits,
                    nfft=args.nfft, kaiser_beta=args.kaiser_beta,
                    averages=tuple(sorted(args.averages)))
    try:
        summary = generate_outputs(args.data_dir, args.output_dir, config, args.csv)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Generated {len(summary['artifacts'])} files in {args.output_dir}")
    print(f"T = {config.block_bits * config.bit_duration * 1e9:g} ns; "
          f"FFT grid = {summary['derived_parameters']['fft_grid_step_hz'] / 1e6:.6f} MHz; "
          f"window ENBW = {summary['derived_parameters']['window_enbw_hz'] / 1e6:.6f} MHz")
    print(f"Completed in {summary['elapsed_seconds']:.2f} seconds.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
