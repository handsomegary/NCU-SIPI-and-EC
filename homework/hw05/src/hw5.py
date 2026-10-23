"""Construct fixed-rate and clock-triggered, vertically normalized eye histograms.

Time-bin functions use picoseconds. CSV time is in seconds. Waveform voltages
are never interpolated: only differential-clock crossing times are interpolated.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import platform

import numpy as np


@dataclass(frozen=True)
class EyeConfig:
    ui_ps: float = 25.0
    bins_per_ui: int = 25
    vmin: float = -0.1
    vmax: float = 1.1
    dv: float = 0.01
    origin_ps: float = 0.0

    def __post_init__(self):
        values = (self.ui_ps, self.vmin, self.vmax, self.dv, self.origin_ps)
        if not all(np.isfinite(value) for value in values):
            raise ValueError("Configuration values must be finite.")
        if self.ui_ps <= 0 or self.dv <= 0 or self.vmax <= self.vmin:
            raise ValueError("UI and dv must be positive, and vmax must exceed vmin.")
        if not isinstance(self.bins_per_ui, (int, np.integer)) or self.bins_per_ui < 1:
            raise ValueError("bins_per_ui must be a positive integer.")
        count = (self.vmax - self.vmin) / self.dv
        if not np.isclose(count, round(count), rtol=0, atol=1e-9):
            raise ValueError("The voltage range must be an integer multiple of dv.")

    @property
    def dt_ps(self):
        return self.ui_ps / self.bins_per_ui

    @property
    def window_ps(self):
        return 2 * self.ui_ps

    @property
    def ntime(self):
        return 2 * self.bins_per_ui

    @property
    def nvolt(self):
        return int(round((self.vmax - self.vmin) / self.dv))


def _finite_1d(values, name):
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite one-dimensional array.")
    return result


def _snap_integer_roundoff(values):
    """Repair only floating-point roundoff near an integer bin boundary."""
    values = np.asarray(values, dtype=float)
    nearest = np.rint(values)
    tolerance = 32 * np.finfo(float).eps * np.maximum(np.abs(values), 1.0)
    return np.where(np.abs(values - nearest) <= tolerance, nearest, values)


def voltage_bin_indices(voltage, config):
    """Floor-based bins; final bin includes exact vmax. Reject out-of-range data."""
    voltage = _finite_1d(voltage, "voltage")
    if np.any(voltage < config.vmin) or np.any(voltage > config.vmax):
        raise ValueError("Voltage samples lie outside [vmin, vmax].")
    scaled = _snap_integer_roundoff((voltage - config.vmin) / config.dv)
    indices = np.floor(scaled).astype(np.int64)
    return np.minimum(indices, config.nvolt - 1)


def fixed_time_bin_indices(time_ps, config):
    """Map absolute time to [0, 2UI), with fixed origin and robust grid edges."""
    time_ps = _finite_1d(time_ps, "time_ps")
    ticks = _snap_integer_roundoff((time_ps - config.origin_ps) / config.dt_ps)
    # Taking floor before integer modulo avoids t % window returning nearly 2UI
    # at decimal CSV timestamps that should lie exactly on the time grid.
    return np.floor(ticks).astype(np.int64) % config.ntime


def positive_zero_crossings(time_ps, differential_clock):
    """Interpolate each transition clk[i-1] < 0 <= clk[i]."""
    time_ps = _finite_1d(time_ps, "time_ps")
    clock = _finite_1d(differential_clock, "differential_clock")
    if time_ps.size != clock.size or time_ps.size < 2:
        raise ValueError("Time and clock must have matching lengths of at least two.")
    if np.any(np.diff(time_ps) <= 0):
        raise ValueError("Time must increase strictly.")
    right = np.flatnonzero((clock[:-1] < 0) & (clock[1:] >= 0)) + 1
    left = right - 1
    fraction = -clock[left] / (clock[right] - clock[left])
    return time_ps[left] + fraction * (time_ps[right] - time_ps[left])


def triggered_time_bin_indices(time_ps, triggers_ps, config):
    """Return indices and validity mask using the latest crossing per sample.

    Times before the first trigger and ages >= 2UI are excluded. Triggered bins
    never wrap. A sample belongs to at most one trigger interval.
    """
    time_ps = _finite_1d(time_ps, "time_ps")
    triggers = _finite_1d(triggers_ps, "triggers_ps")
    if triggers.size and np.any(np.diff(triggers) <= 0):
        raise ValueError("Trigger times must increase strictly.")
    indices = np.full(time_ps.size, -1, dtype=np.int64)
    if not triggers.size:
        return indices, np.zeros(time_ps.size, dtype=bool)
    latest = np.searchsorted(triggers, time_ps, side="right") - 1
    available = latest >= 0
    age = time_ps[available] - triggers[latest[available]]
    ticks = _snap_integer_roundoff(age / config.dt_ps)
    indices[available] = np.floor(ticks).astype(np.int64)
    valid = available & (indices >= 0) & (indices < config.ntime)
    return indices, valid


def histogram(voltage_indices, time_indices, config):
    voltage_indices = np.asarray(voltage_indices, dtype=np.int64)
    time_indices = np.asarray(time_indices, dtype=np.int64)
    if voltage_indices.shape != time_indices.shape or voltage_indices.ndim != 1:
        raise ValueError("Histogram indices must be matching one-dimensional arrays.")
    if np.any(voltage_indices < 0) or np.any(voltage_indices >= config.nvolt):
        raise ValueError("Voltage-bin index outside histogram.")
    if np.any(time_indices < 0) or np.any(time_indices >= config.ntime):
        raise ValueError("Time-bin index outside histogram.")
    flat = voltage_indices * config.ntime + time_indices
    return np.bincount(flat, minlength=config.nvolt * config.ntime).reshape(
        config.nvolt, config.ntime
    )


def vertical_normalize(counts):
    """P(voltage bin | time bin); empty columns remain all zero, with no NaNs."""
    counts = np.asarray(counts)
    if counts.ndim != 2 or not np.all(np.isfinite(counts)) or np.any(counts < 0):
        raise ValueError("Histogram counts must be finite, nonnegative, and two-dimensional.")
    totals = counts.sum(axis=0)
    probability = np.divide(
        counts, totals[np.newaxis, :],
        out=np.zeros(counts.shape, dtype=float), where=totals[np.newaxis, :] > 0,
    )
    return probability, totals


def read_waveforms(path):
    with Path(path).open("r", encoding="utf-8-sig") as stream:
        header = stream.readline().strip().split(",")
    if header != ["Time", "Vo", "CLKp", "CLKn"]:
        raise ValueError("Expected CSV columns: Time,Vo,CLKp,CLKn.")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if data.shape[1] != 4 or data.shape[0] < 2 or not np.all(np.isfinite(data)):
        raise ValueError("CSV must contain at least two rows of four finite values.")
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError("CSV time values must increase strictly.")
    return data


def _normalization_metrics(probability, totals):
    occupied = totals > 0
    return {
        "total_samples": int(totals.sum()),
        "nonempty_time_bins": int(np.count_nonzero(occupied)),
        "empty_time_bins": int(np.count_nonzero(~occupied)),
        "minimum_column_count": int(totals.min()),
        "maximum_column_count": int(totals.max()),
        "maximum_column_sum_error": float(np.max(np.abs(
            probability[:, occupied].sum(axis=0) - 1.0
        ))) if np.any(occupied) else 0.0,
        "all_probabilities_finite": bool(np.all(np.isfinite(probability))),
    }


def analyze(data, config):
    """Analyze a validated four-column matrix; return arrays and summary metrics."""
    data = np.asarray(data, dtype=float)
    if data.ndim != 2 or data.shape[1] != 4 or data.shape[0] < 2:
        raise ValueError("Expected at least two rows of [Time_s, Vo, CLKp, CLKn].")
    if not np.all(np.isfinite(data)) or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError("Waveforms must be finite with strictly increasing time.")
    time_ps = data[:, 0] * 1e12
    voltage_indices = voltage_bin_indices(data[:, 1], config)
    fixed_indices = fixed_time_bin_indices(time_ps, config)
    triggers = positive_zero_crossings(time_ps, data[:, 2] - data[:, 3])
    triggered_indices, valid = triggered_time_bin_indices(time_ps, triggers, config)
    fixed_counts = histogram(voltage_indices, fixed_indices, config)
    triggered_counts = histogram(voltage_indices[valid], triggered_indices[valid], config)
    fixed_probability, fixed_totals = vertical_normalize(fixed_counts)
    triggered_probability, triggered_totals = vertical_normalize(triggered_counts)
    intervals = np.diff(triggers)
    phase_error = (
        triggers - (triggers[0] + np.arange(triggers.size) * config.window_ps)
        if triggers.size else np.array([])
    )
    spacing = np.diff(time_ps)
    before_first = time_ps < triggers[0] if triggers.size else np.ones(time_ps.size, bool)
    arrays = {
        "fixed_counts": fixed_counts,
        "triggered_counts": triggered_counts,
        "fixed_probability": fixed_probability,
        "triggered_probability": triggered_probability,
        "fixed_column_counts": fixed_totals,
        "triggered_column_counts": triggered_totals,
        "time_edges_ps": np.linspace(0, config.window_ps, config.ntime + 1),
        "voltage_edges_V": np.linspace(config.vmin, config.vmax, config.nvolt + 1),
        "trigger_times_ps": triggers,
    }
    summary = {
        "configuration": {
            **asdict(config), "dt_ps": config.dt_ps,
            "window_ps": config.window_ps, "time_bins": config.ntime,
            "voltage_bins": config.nvolt,
        },
        "input": {
            "samples": int(time_ps.size), "start_time_ps": float(time_ps[0]),
            "end_time_ps": float(time_ps[-1]),
            "duration_ps": float(time_ps[-1] - time_ps[0]),
            "minimum_voltage_V": float(data[:, 1].min()),
            "maximum_voltage_V": float(data[:, 1].max()),
            "median_sample_spacing_ps": float(np.median(spacing)),
            "minimum_sample_spacing_ps": float(spacing.min()),
            "maximum_sample_spacing_ps": float(spacing.max()),
        },
        "trigger": {
            "polarity": "positive differential zero crossing",
            "crossing_method": "linear interpolation between adjacent clock samples",
            "crossings": int(triggers.size),
            "first_crossing_ps": float(triggers[0]) if triggers.size else None,
            "last_crossing_ps": float(triggers[-1]) if triggers.size else None,
            "minimum_interval_ps": float(intervals.min()) if intervals.size else None,
            "maximum_interval_ps": float(intervals.max()) if intervals.size else None,
            "mean_interval_ps": float(intervals.mean()) if intervals.size else None,
            "phase_error_reference": "first crossing plus n * nominal 2UI",
            "minimum_phase_error_ps": float(phase_error.min()) if triggers.size else None,
            "maximum_phase_error_ps": float(phase_error.max()) if triggers.size else None,
            "peak_to_peak_phase_error_ps": float(np.ptp(phase_error)) if triggers.size else None,
            "samples_before_first_trigger": int(before_first.sum()),
            "samples_excluded_after_2UI": int((~valid & ~before_first).sum()),
        },
        "fixed_rate": _normalization_metrics(fixed_probability, fixed_totals),
        "triggered": _normalization_metrics(triggered_probability, triggered_totals),
        "probability_definition": "P(voltage bin | time bin), not a continuous probability density",
        "time_bin_convention": "[left, right); fixed origin 0 ps by default; triggered age in [0, 2UI)",
        "voltage_bin_convention": "floor((V-vmin)/dv); [left, right), final bin includes vmax",
        "trigger_assignment": "latest positive crossing; each measured sample is counted at most once",
        "voltage_interpolation": False,
    }
    return arrays, summary


def create_plots(data, arrays, config, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    maximum = max(arrays["fixed_probability"].max(), arrays["triggered_probability"].max())
    maximum = float(maximum) if maximum > 0 else 1.0
    extent = [0, config.window_ps, config.vmin, config.vmax]

    def draw_eye(axis, probability, title):
        image = axis.imshow(
            probability, origin="lower", extent=extent, aspect="auto",
            interpolation="nearest", cmap="magma", vmin=0, vmax=maximum,
        )
        axis.set(xlabel="Time within 2 UI (ps)", ylabel="Output voltage (V)", title=title)
        return image

    for key, filename, title in (
        ("fixed_probability", "fixed_rate_eye.png", "Fixed-rate statistical eye"),
        ("triggered_probability", "triggered_eye.png", "Triggered statistical eye"),
    ):
        figure, axis = plt.subplots(figsize=(8, 4.7), layout="constrained")
        image = draw_eye(axis, arrays[key], title)
        figure.colorbar(image, ax=axis, label="Conditional bin probability")
        figure.savefig(output_dir / filename, dpi=200, metadata={"Software": "HW5 eye analysis"})
        plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.6), layout="constrained")
    draw_eye(axes[0], arrays["fixed_probability"], "Fixed rate: absolute nominal clock")
    image = draw_eye(axes[1], arrays["triggered_probability"], "Triggered: latest differential crossing")
    figure.colorbar(image, ax=axes, label="Conditional bin probability", shrink=0.9)
    figure.savefig(output_dir / "eye_comparison.png", dpi=200, metadata={"Software": "HW5 eye analysis"})
    plt.close(figure)

    time_ps = data[:, 0] * 1e12
    triggers = arrays["trigger_times_ps"]
    start = triggers[0] if triggers.size else time_ps[0]
    stop = start + 6 * config.window_ps
    segment = (time_ps >= start - config.window_ps / 2) & (time_ps <= stop)
    figure, axes = plt.subplots(3, 1, figsize=(10, 8), layout="constrained")
    axes[0].plot(time_ps[segment] - start, (data[:, 2] - data[:, 3])[segment], color="#3568a9")
    axes[0].axhline(0, color="gray", linewidth=0.7)
    axes[0].set(ylabel="CLKp - CLKn (V)", title="Positive clock crossings define time zero for each interval")
    axes[1].plot(time_ps[segment] - start, data[segment, 1], color="#c46235")
    axes[1].set(xlabel="Time relative to first trigger (ps)", ylabel="Output voltage (V)")
    for trigger in triggers[(triggers >= start) & (triggers <= stop)]:
        for axis in axes[:2]:
            axis.axvline(trigger - start, color="#327955", alpha=0.5, linewidth=0.8)
    if triggers.size:
        phase_error = triggers - (triggers[0] + np.arange(triggers.size) * config.window_ps)
        axes[2].plot((triggers - triggers[0]) / 1000, phase_error, color="#7254a3", linewidth=1)
    axes[2].set(
        xlabel="Elapsed time since first trigger (ns)",
        ylabel="Clock phase error (ps)",
        title=f"Crossing time relative to nominal {config.window_ps:g} ps intervals",
    )
    figure.savefig(output_dir / "clock_alignment.png", dpi=200, metadata={"Software": "HW5 eye analysis"})
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(9, 2.7), layout="constrained")
    if triggers.size:
        axis.plot((triggers - triggers[0]) / 1000, phase_error, color="#7254a3", linewidth=1)
    axis.set(
        xlabel="Elapsed time since first trigger (ns)",
        ylabel="Clock phase error (ps)",
        title=f"Received crossings relative to nominal {config.window_ps:g} ps intervals",
    )
    figure.savefig(output_dir / "clock_phase.png", dpi=200, metadata={"Software": "HW5 eye analysis"})
    plt.close(figure)


def main(argv=None):
    base = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=base / "instructions" / "HW5_ADS_waveforms.csv")
    parser.add_argument("--output-dir", type=Path, default=base / "results")
    parser.add_argument("--ui-ps", type=float, default=25.0)
    parser.add_argument("--bins-per-ui", type=int, default=25)
    parser.add_argument("--vmin", type=float, default=-0.1)
    parser.add_argument("--vmax", type=float, default=1.1)
    parser.add_argument("--dv", type=float, default=0.01)
    parser.add_argument("--origin-ps", type=float, default=0.0)
    parser.add_argument("--no-plots", action="store_true", help="Only save numerical outputs.")
    args = parser.parse_args(argv)
    try:
        config = EyeConfig(args.ui_ps, args.bins_per_ui, args.vmin, args.vmax, args.dv, args.origin_ps)
        data = read_waveforms(args.input)
        arrays, summary = analyze(data, config)
        if summary["trigger"]["crossings"] == 0:
            raise ValueError("No positive differential-clock crossings were found.")
    except (ValueError, OSError) as error:
        parser.error(str(error))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary["input"]["filename"] = args.input.name
    summary["input"]["sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    summary["software"] = {"python": platform.python_version(), "numpy": np.__version__}
    np.savez_compressed(args.output_dir / "histograms.npz", **arrays)
    if not args.no_plots:
        create_plots(data, arrays, config, args.output_dir)
        import matplotlib
        summary["software"]["matplotlib"] = matplotlib.__version__
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n"
    )
    print(json.dumps({
        "samples": summary["input"]["samples"],
        "crossings": summary["trigger"]["crossings"],
        "fixed_count": summary["fixed_rate"]["total_samples"],
        "triggered_count": summary["triggered"]["total_samples"],
        "output": args.output_dir.name,
    }, indent=2))


if __name__ == "__main__":
    main()
