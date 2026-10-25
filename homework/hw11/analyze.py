"""Reproducible HW11 eight-channel crosstalk simulation.

H has axes (frequency, receiver, transmitter). The loaded channel includes
the external source/load RC components. Impulse coefficients have voltage
gain units (not volts/second), so discrete convolution uses no extra dt.
Finite-bandwidth pre-ringing is retained as explicit negative impulse lags.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import re
import time

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.special import i0


def load_touchstone(path):
    """Read full-matrix Touchstone 1.x RI/MA/DB; return Hz, S[f,i,j], Z0.

    Three-or-more-port files use row-major S11,S12,... ordering. Two-port
    Touchstone uses its special S11,S21,S12,S22 ordering. No reciprocity
    assumption or transpose is applied. Touchstone 2.x is rejected explicitly.
    """
    path = Path(path)
    match = re.search(r"\.s(\d+)p$", path.name, re.I)
    if not match:
        raise ValueError("The filename must end in .sNp")
    nport = int(match.group(1))
    options = ["ghz", "s", "ma", "r", "50"]
    numbers = []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            options = line[1:].lower().split()
            continue
        if line.startswith("["):
            raise ValueError("Only full-matrix Touchstone 1.x is supported")
        numbers.extend(float(v.replace("D", "E").replace("d", "e")) for v in line.split())
    width = 1 + 2 * nport * nport
    if not numbers or len(numbers) % width:
        raise ValueError("Incomplete full-matrix Touchstone frequency record")
    rows = np.asarray(numbers, float).reshape(-1, width)
    unit = next((v for v in options if v in {"hz", "khz", "mhz", "ghz"}), "ghz")
    fmt = next((v for v in options if v in {"ri", "ma", "db"}), "ma")
    if "s" not in options:
        raise ValueError("Only S parameters are supported")
    z0 = float(options[options.index("r") + 1]) if "r" in options else 50.0
    freq = rows[:, 0] * {"hz": 1, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}[unit]
    if not np.isfinite(rows).all() or z0 <= 0 or freq[0] < 0 or np.any(np.diff(freq) <= 0):
        raise ValueError("Finite increasing nonnegative frequencies and positive Z0 required")
    pair = rows[:, 1:].reshape(-1, nport * nport, 2)
    if fmt == "ri":
        values = pair[..., 0] + 1j * pair[..., 1]
    else:
        magnitude = 10 ** (pair[..., 0] / 20) if fmt == "db" else pair[..., 0]
        values = magnitude * np.exp(1j * np.deg2rad(pair[..., 1]))
    s = values.reshape(-1, nport, nport)
    if nport == 2:
        s = s.transpose(0, 2, 1).copy()
    return freq, s, z0


def solve_transfer(frequency, s, z0, rs=30.0, cs=1.5e-12, rl=65.0, cl=1e-12):
    """Unit-source MNA solution H[f,Rx,Tx], including external RC loads.

    Port currents point into the S network. The equations are Yv+i=J and
    (I-S)v-Z0(I+S)i=0. Odd numbered ports are Tx; even ports are Rx.
    """
    f, s = np.asarray(frequency, float), np.asarray(s, complex)
    if s.ndim != 3 or s.shape[0] != f.size or s.shape[1] != s.shape[2] or s.shape[1] % 2:
        raise ValueError("An even-port S matrix with axes (frequency,port,port) is required")
    if rs <= 0 or rl <= 0 or cs < 0 or cl < 0 or z0 <= 0:
        raise ValueError("Positive resistances/reference and nonnegative capacitances required")
    nport = s.shape[1]
    tx, rx = np.arange(0, nport, 2), np.arange(1, nport, 2)
    identity = np.eye(nport)
    base = np.zeros((2 * nport, 2 * nport), complex)
    base[tx, tx], base[rx, rx] = 1 / rs, 1 / rl
    base[:nport, nport:] = identity
    b = np.zeros((2 * nport, tx.size), complex)
    b[tx, np.arange(tx.size)] = 1 / rs
    h = np.empty((f.size, rx.size, tx.size), complex)
    for k, fk in enumerate(f):
        a = base.copy()
        a[tx, tx] += 2j * np.pi * fk * cs
        a[rx, rx] += 2j * np.pi * fk * cl
        a[nport:, :nport] = identity - s[k]
        a[nport:, nport:] = -z0 * (identity + s[k])
        h[k] = np.linalg.solve(a, b)[rx]
    return h


@dataclass
class ImpulseKernel:
    samples: np.ndarray
    first_lag: int
    dt: float
    diagnostics: dict = field(default_factory=dict)


def frequency_to_impulse(freq, h, dt=2e-12, kaiser_a=1.0, pre_time=5e-9, post_time=20e-9):
    """Real Hermitian inverse DFT, retaining lags [-pre_time,post_time].

    Kaiser W(f)=I0(pi*a*sqrt(1-(f/fmax)^2))/I0(pi*a). DC is not halved:
    irfft already accounts for conjugate negative frequencies. Frequencies
    above the measured band are zero. Negative lags represent the finite
    bandwidth approximation, not physically acausal propagation.
    """
    f, h = np.asarray(freq, float), np.asarray(h, complex)
    if f.ndim != 1 or f.size < 2 or h.shape[0] != f.size or h.ndim != 3 or dt <= 0:
        raise ValueError("Uniform frequency vector and H[f,Rx,Tx] required")
    df = f[1] - f[0]
    if f[0] != 0 or df <= 0 or not np.allclose(np.diff(f), df, rtol=1e-9, atol=1e-5):
        raise ValueError("Frequency samples must include DC with a uniform positive step")
    nfft = round(1 / (dt * df))
    if nfft < 2 or not np.isclose(nfft * dt * df, 1, rtol=1e-10):
        raise ValueError("dt and df must define an integral inverse-transform length")
    if f[-1] > 0.5 / dt or f.size > nfft // 2 + 1 or kaiser_a < 0:
        raise ValueError("Invalid band limit or Kaiser parameter")
    npre, npost = round(pre_time / dt), round(post_time / dt)
    if min(npre, npost) < 0 or npre >= nfft // 2 or npost >= nfft // 2:
        raise ValueError("Support must lie within the +/- half-period unambiguous lag range")
    window = i0(np.pi * kaiser_a * np.sqrt(np.maximum(0, 1 - (f / f[-1]) ** 2))) / i0(np.pi * kaiser_a)
    spectrum = np.zeros((nfft // 2 + 1, *h.shape[1:]), complex)
    spectrum[:f.size] = h * window[:, None, None]
    spectrum[0] = spectrum[0].real
    if nfft % 2 == 0:
        spectrum[-1] = spectrum[-1].real
    circular = irfft(spectrum, nfft, axis=0)
    indices = np.arange(-npre, npost + 1) % nfft
    samples = circular[indices]
    discarded = circular.copy()
    discarded[indices] = 0
    dc_error = samples.sum(axis=0) - h[0].real
    main_gain = np.abs(np.diag(h[0])) if h.shape[1] == h.shape[2] else np.max(np.abs(h[0]), axis=1)
    relative_dc_error = np.divide(np.abs(dc_error), main_gain[:, None],
                                  out=np.full_like(dc_error.real, np.nan), where=main_gain[:, None] > 0)
    diagnostics = {
        "dft_length": nfft, "dft_period_s": nfft * dt,
        "first_lag_s": -npre * dt, "last_lag_s": npost * dt,
        "coefficient_units": "dimensionless discrete-convolution weights",
        "pre_ringing_policy": "Retain negative lags explicitly and align linear convolution by first_lag.",
        "discarded_l1_by_path": np.abs(discarded).sum(axis=0).tolist(),
        "discarded_l2_by_path": np.sqrt(np.square(discarded).sum(axis=0)).tolist(),
        "discarded_peak_by_path": np.abs(discarded).max(axis=0).tolist(),
        "dc_sum_error_by_path": dc_error.tolist(),
        "dc_sum_max_abs_error": float(np.abs(dc_error).max()),
        "dc_sum_abs_error_relative_receiver_main_dc_gain_by_path":
            [[None if not np.isfinite(v) else float(v) for v in row] for row in relative_dc_error],
        "dc_sum_max_error_relative_receiver_main_dc_gain":
            float(np.nanmax(relative_dc_error)) if np.isfinite(relative_dc_error).any() else None,
        "negative_lag_energy_fraction": float(np.square(samples[:npre]).sum() / np.square(circular).sum()),
    }
    return ImpulseKernel(samples, -npre, dt, diagnostics)


def independent_bits(n_bits, n_tx, seed=2026):
    """Deterministic, independent pseudorandom streams using spawned sequences."""
    if n_bits < 1 or n_tx < 1:
        raise ValueError("Positive bit and transmitter counts required")
    return np.column_stack([np.random.default_rng(s).integers(0, 2, n_bits, dtype=np.uint8)
                            for s in np.random.SeedSequence(seed).spawn(n_tx)])


def nrz_waveform(bits, npui, rise_samples, delay_samples=0, start=0, stop=None, vpp=1.0, *, _validate=True):
    """Sample finite-rise 0..Vpp NRZ streams at absolute input sample indices.

    Bits have shape (bit,Tx), or a one-dimensional bit stream. The first ramp
    sample has the old value; after rise_samples the new level is reached.
    Before and after the provided bit stream the level is zero. Channel
    delays must be integer sample counts, which holds for this assignment.
    """
    bits = np.asarray(bits)
    one_dim = bits.ndim == 1
    if one_dim:
        bits = bits[:, None]
    if bits.ndim != 2 or npui < 1 or not 1 <= rise_samples <= npui or (_validate and not np.isin(bits, [0, 1]).all()):
        raise ValueError("Binary bit columns and 1 <= rise_samples <= npui required")
    delays = np.broadcast_to(np.asarray(delay_samples, int), (bits.shape[1],))
    if (delays < 0).any():
        raise ValueError("Nonnegative channel delays required")
    stop = bits.shape[0] * npui + int(delays.max()) if stop is None else int(stop)
    if stop < start:
        raise ValueError("Invalid sample range")
    positions = np.arange(start, stop)[:, None] - delays[None]
    index, phase = np.divmod(positions, npui)
    valid = (index >= 0) & (index < bits.shape[0])
    valid_prev = (index > 0) & (index <= bits.shape[0])
    columns = np.arange(bits.shape[1])[None]
    current = np.where(valid, bits[np.clip(index, 0, bits.shape[0] - 1), columns], 0).astype(float)
    previous = np.where(valid_prev, bits[np.clip(index - 1, 0, bits.shape[0] - 1), columns], 0)
    fraction = np.minimum(phase / rise_samples, 1.0)
    output = np.where(valid, vpp * (previous + (current - previous) * fraction), 0.0)
    return output[:, 0] if one_dim else output


class OverlapAdd:
    """Stateful real MIMO linear convolution; process then flush once."""
    def __init__(self, kernel, block_size=32768, dtype=np.float64):
        coefficients = kernel.samples if isinstance(kernel, ImpulseKernel) else np.asarray(kernel)
        if coefficients.ndim == 1:
            coefficients = coefficients[:, None, None]
        if coefficients.ndim != 3 or len(coefficients) < 1 or block_size < 1:
            raise ValueError("Kernel shape (lag,Rx,Tx) and positive block size required")
        self.dtype, self.block_size = np.dtype(dtype), int(block_size)
        self.length, self.nout, self.nin = coefficients.shape
        self.nfft = next_fast_len(self.block_size + self.length - 1)
        self.spectrum = rfft(coefficients.astype(self.dtype), self.nfft, axis=0)
        self.tail = np.zeros((self.length - 1, self.nout), self.dtype)
        self.finished = False

    def process(self, x):
        x = np.asarray(x, dtype=self.dtype)
        if x.ndim == 1:
            x = x[:, None]
        if self.finished or x.ndim != 2 or x.shape[1] != self.nin or not 0 < len(x) <= self.block_size:
            raise ValueError("Invalid convolution block or convolution already flushed")
        transformed = rfft(x, self.nfft, axis=0)
        product = np.einsum("foi,fi->fo", self.spectrum, transformed, optimize=False)
        block = irfft(product, self.nfft, axis=0)[:len(x) + self.length - 1]
        block[:len(self.tail)] += self.tail
        output = block[:len(x)].copy()
        self.tail = block[len(x):].copy()
        return output

    def flush(self):
        if self.finished:
            raise ValueError("Convolution was already flushed")
        self.finished = True
        return self.tail.copy()


def overlap_add_mimo(x, kernel, block_size=8192):
    """Full raw linear convolution (first output has kernel.first_lag time).

    A 2D result (sample,Rx) is returned, including the full trailing response.
    Input and kernel are not shifted or silently truncated.
    """
    x = np.asarray(x)
    processor = OverlapAdd(kernel, block_size)
    chunks = [processor.process(x[k:k + block_size]) for k in range(0, len(x), block_size)]
    return np.concatenate(chunks + [processor.flush()], axis=0)


class EyeHistogram:
    """Streaming fixed-rate eye; counts axes are (output,voltage,time phase)."""
    def __init__(self, npui, n_outputs, voltage_min=-0.25, voltage_max=1.25, voltage_bins=3000):
        if npui < 1 or n_outputs < 1 or voltage_bins < 2 or voltage_max <= voltage_min:
            raise ValueError("Invalid eye grid")
        self.npui, self.n_outputs = int(npui), int(n_outputs)
        self.edges = np.linspace(voltage_min, voltage_max, voltage_bins + 1)
        self.counts = np.zeros((n_outputs, voltage_bins, npui), np.int64)
        self.outside = np.zeros(n_outputs, np.int64)

    def update(self, y, start_sample=0):
        y = np.asarray(y)
        if y.ndim == 1:
            y = y[:, None]
        if y.ndim != 2 or y.shape[1] != self.n_outputs or not np.isfinite(y).all():
            raise ValueError("Finite eye data of shape (sample,output) required")
        phase = (np.arange(len(y)) + int(start_sample)) % self.npui
        for channel in range(self.n_outputs):
            voltage = np.floor((y[:, channel] - self.edges[0]) / (self.edges[1] - self.edges[0])).astype(int)
            valid = (voltage >= 0) & (voltage < self.counts.shape[1])
            self.outside[channel] += np.count_nonzero(~valid)
            flat = voltage[valid] * self.npui + phase[valid]
            self.counts[channel] += np.bincount(flat, minlength=self.counts[channel].size).reshape(self.counts.shape[1:])

    def probabilities(self):
        totals = self.counts.sum(axis=1, keepdims=True)
        return np.divide(self.counts, totals, out=np.zeros_like(self.counts, float), where=totals != 0)

    def metrics(self, thresholds):
        """Conditional 0.1%/99.9% eye opening; an empirical finite-sample metric.

        At each phase split by the provided decision voltage. Estimate the
        99.9th percentile of the lower cluster and 0.1th percentile of the
        upper cluster. Select the largest gap, requiring both clusters to
        contain at least 20% of phase samples. Resolution is one voltage bin.
        This is not a BER prediction or a guaranteed worst-case eye opening.
        """
        centers = (self.edges[:-1] + self.edges[1:]) / 2
        thresholds = np.broadcast_to(thresholds, (self.n_outputs,))
        result = []
        for channel, threshold in enumerate(thresholds):
            candidates = []
            lower = centers < threshold
            for phase in range(self.npui):
                count = self.counts[channel, :, phase]
                total = count.sum()
                low, high = count[lower], count[~lower]
                if total == 0 or min(low.sum(), high.sum()) < 0.2 * total:
                    continue
                low_idx = min(np.searchsorted(np.cumsum(low), 0.999 * low.sum()), low.size - 1)
                high_idx = min(np.searchsorted(np.cumsum(high), 0.001 * high.sum()), high.size - 1)
                low_v, high_v = centers[lower][low_idx], centers[~lower][high_idx]
                candidates.append((high_v - low_v, phase, low_v, high_v))
            opening, phase, low_v, high_v = max(candidates) if candidates else (0.0, 0, None, None)
            result.append({"opening_v": float(opening), "best_phase_ui": phase / self.npui,
                           "lower_cluster_99_9_percentile_v": None if low_v is None else float(low_v),
                           "upper_cluster_0_1_percentile_v": None if high_v is None else float(high_v),
                           "threshold_v": float(threshold), "samples": int(self.counts[channel].sum()),
                           "outside_voltage_range": int(self.outside[channel]),
                           "minimum_samples_per_phase": int(self.counts[channel].sum(axis=0).min())})
        return result


def _slice_kernel(kernel, receivers, transmitters):
    return ImpulseKernel(kernel.samples[:, receivers][:, :, transmitters], kernel.first_lag, kernel.dt)


def simulate_case(bits, kernel, npui, rise_samples, delays, n_bits, guard_ui, receivers, transmitters,
                  block_size=32768, snapshot_start=145e-9, snapshot_duration=4e-9):
    """Stream a guarded case, retaining exactly n_bits full UIs per receiver."""
    selected = _slice_kernel(kernel, receivers, transmitters)
    processor = OverlapAdd(selected, block_size, dtype=np.float32)
    eye = EyeHistogram(npui, len(receivers))
    begin = guard_ui * npui - kernel.first_lag
    end = begin + n_bits * npui
    # Raw linear-convolution index r has physical sample r+first_lag.
    # Generate enough post-roll to include negative-lag pre-ringing at end.
    n_input = end
    snap_begin = int(round(snapshot_start / kernel.dt))
    snap_end = snap_begin + int(round(snapshot_duration / kernel.dt))
    snapshot, snapshot_indices = [], []
    sum_y = np.zeros(len(receivers))
    minimum = np.full(len(receivers), np.inf)
    maximum = np.full(len(receivers), -np.inf)
    selected_bits = bits[:, transmitters]
    for start in range(0, n_input, block_size):
        stop = min(start + block_size, n_input)
        # Avoid validating/rebuilding all bit streams in every short block.
        x = nrz_waveform(selected_bits, npui, rise_samples, np.asarray(delays)[transmitters], start, stop, _validate=False)
        y = processor.process(x)
        left, right = max(start, begin), min(stop, end)
        if right <= left:
            continue
        values = y[left - start:right - start]
        relative = left - begin
        eye.update(values, relative)
        sum_y += values.sum(axis=0, dtype=float)
        minimum = np.minimum(minimum, values.min(axis=0))
        maximum = np.maximum(maximum, values.max(axis=0))
        sl, sr = max(relative, snap_begin), min(right - begin, snap_end)
        if sr > sl:
            snapshot.append(values[sl - relative:sr - relative].copy())
            snapshot_indices.append(np.arange(sl, sr))
    if eye.outside.any():
        raise RuntimeError("Eye voltage grid clipped samples; expand its range before publishing")
    extra = {"mean_v": (sum_y / (n_bits * npui)).tolist(), "minimum_v": minimum.tolist(),
             "maximum_v": maximum.tolist(), "retained_samples_per_receiver": n_bits * npui,
             "snapshot_start_s": snapshot_start}
    values = np.concatenate(snapshot) if snapshot else np.empty((0, len(receivers)), np.float32)
    times = np.concatenate(snapshot_indices) * kernel.dt if snapshot_indices else np.empty(0)
    return eye, times, values, extra


def _convergence_check(freq, h, kernel, bits, npui, rise_samples, delays, guard_ui, block_size, kaiser_a):
    """All 64 paths and eight active streams checked against larger support."""
    half_period = 0.5 / (freq[1] - freq[0]) - kernel.dt
    larger = frequency_to_impulse(freq, h, kernel.dt, kaiser_a,
                                  min(2 * abs(kernel.first_lag) * kernel.dt, half_period),
                                  min(2 * (len(kernel.samples) + kernel.first_lag - 1) * kernel.dt, half_period))
    n_sample = min(2048, bits.shape[0] - 2 * guard_ui) * npui
    x = nrz_waveform(bits[:], npui, rise_samples, delays, stop=(2 * guard_ui) * npui + n_sample)
    y = overlap_add_mimo(x, kernel, block_size)
    reference = overlap_add_mimo(x, larger, block_size)
    physical_start = guard_ui * npui
    retained = y[physical_start - kernel.first_lag:physical_start - kernel.first_lag + n_sample]
    expected = reference[physical_start - larger.first_lag:physical_start - larger.first_lag + n_sample]
    delta = retained - expected
    main_dc_gain = np.abs(np.diag(h[0]))
    return {"reference_pre_time_s": -larger.first_lag * kernel.dt,
            "reference_post_time_s": (len(larger.samples) + larger.first_lag - 1) * kernel.dt,
            "test_bits": n_sample // npui, "all_eight_transmitters_active": True,
            "maximum_abs_waveform_difference_v": float(np.abs(delta).max()),
            "rms_waveform_difference_v": float(np.sqrt(np.mean(delta ** 2))),
            "maximum_abs_by_receiver_v": np.abs(delta).max(axis=0).tolist(),
            "rms_by_receiver_v": np.sqrt(np.mean(delta ** 2, axis=0)).tolist(),
            "maximum_abs_by_receiver_relative_main_dc_gain": (np.abs(delta).max(axis=0) / main_dc_gain).tolist(),
            "reference_dc_sum_max_abs_error": larger.diagnostics["dc_sum_max_abs_error"]}


def _sampling_check(freq, h, kernel, bits, npui, rise_samples, delays, guard_ui, block_size, kaiser_a):
    """Short all-active record at twice sample rate, compared on common times."""
    fine = frequency_to_impulse(freq, h, kernel.dt / 2, kaiser_a,
                                 -kernel.first_lag * kernel.dt,
                                 (len(kernel.samples) + kernel.first_lag - 1) * kernel.dt)
    n_sample = 512 * npui
    total = 2 * guard_ui * npui + n_sample
    coarse_x = nrz_waveform(bits, npui, rise_samples, delays, stop=total)
    fine_x = nrz_waveform(bits, 2 * npui, 2 * rise_samples, 2 * delays, stop=2 * total)
    coarse_y = overlap_add_mimo(coarse_x, kernel, block_size)
    fine_y = overlap_add_mimo(fine_x, fine, 2 * block_size)
    coarse_start = guard_ui * npui - kernel.first_lag
    fine_start = 2 * guard_ui * npui - fine.first_lag
    delta = coarse_y[coarse_start:coarse_start + n_sample] - fine_y[fine_start:fine_start + 2 * n_sample:2]
    return {"coarse_samples_per_ui": npui, "fine_samples_per_ui": 2 * npui,
            "test_bits": 512, "all_eight_transmitters_active": True,
            "maximum_abs_waveform_difference_v": float(np.abs(delta).max()),
            "rms_waveform_difference_v": float(np.sqrt(np.mean(delta ** 2))),
            "maximum_abs_by_receiver_v": np.abs(delta).max(axis=0).tolist(),
            "maximum_abs_by_receiver_relative_main_dc_gain": (np.abs(delta).max(axis=0) / np.abs(np.diag(h[0]))).tolist(),
            "rms_by_receiver_v": np.sqrt(np.mean(delta ** 2, axis=0)).tolist()}


def _plots(output, freq, s, eyes, times, snapshot, kernel):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "figure.dpi": 130})
    fig, axs = plt.subplots(4, 4, figsize=(13, 10), sharex=True)
    for j, ax in enumerate(axs.flat):
        ax.plot(freq / 1e9, 20 * np.log10(np.maximum(np.abs(s[:, 0, j]), 1e-15)))
        ax.set_title(f"S(1,{j + 1})")
        ax.set_xlabel("Frequency (GHz)")
        ax.set_ylabel("Magnitude (dB)")
    fig.suptitle("All sixteen S parameters in output-port row 1")
    fig.tight_layout()
    fig.savefig(figures / "s_parameters.png")
    plt.close(fig)

    # Use the same complete observed voltage range for every comparison.
    # Histogram bounds are independently checked against clipping upstream.
    occupied = np.zeros(next(iter(eyes.values())).counts.shape[1], dtype=bool)
    for eye in eyes.values():
        occupied |= eye.counts.sum(axis=(0, 2)) > 0
    grid = next(iter(eyes.values())).edges
    first, last = np.flatnonzero(occupied)[[0, -1]]
    pad = max(.02, .02 * (grid[last + 1] - grid[first]))
    voltage_limits = (max(grid[0], grid[first] - pad), min(grid[-1], grid[last + 1] + pad))

    def draw_eye(ax, eye, channel, title):
        probability = eye.probabilities()[channel]
        image = ax.imshow(np.log10(np.maximum(np.hstack([probability, probability]), 1e-8)),
                          origin="lower", extent=[0, 2, eye.edges[0], eye.edges[-1]],
                          aspect="auto", vmin=-6, vmax=-1, cmap="viridis", interpolation="nearest")
        ax.set(title=title, xlabel="Time (UI)", ylabel="Voltage (V)", ylim=voltage_limits)
        return image

    all_on = eyes["all_on"]
    fig, axs = plt.subplots(2, 4, figsize=(14, 7), sharex=True, sharey=True)
    for channel, ax in enumerate(axs.flat):
        im = draw_eye(ax, all_on, channel, f"Rx{channel}: all Tx on")
    fig.subplots_adjust(left=.06, right=0.87, bottom=.09, top=.93, wspace=0.25, hspace=0.35)
    colorbar_axis = fig.add_axes([.90, .17, .015, .65])
    fig.colorbar(im, cax=colorbar_axis, label="log10 probability per voltage bin")
    fig.savefig(figures / "eyes_all_on.png")
    plt.close(fig)
    for case, channel, name, title in [
        ("only_tx0", 0, "eye_rx0_only_tx0", "Rx0: only Tx0 on"),
        ("all_on", 0, "eye_rx0_all_on", "Rx0: all Tx on"),
        ("only_tx4", 3, "eye_rx4_only_tx4", "Rx4: only Tx4 on"),
        ("all_on", 4, "eye_rx4_all_on", "Rx4: all Tx on"),
    ]:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        im = draw_eye(ax, eyes[case], channel, title)
        fig.colorbar(im, ax=ax, label="log10 probability per voltage bin")
        fig.tight_layout()
        fig.savefig(figures / (name + ".png"))
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(12, 4))
    for j, rx in enumerate([1, 2, 3, 4]):
        scale = 1 if rx == 4 else 5
        ax.plot(times * 1e9, snapshot[:, j] * scale, label=f"Rx{rx}" + (" (FEXT x5)" if scale == 5 else " (main)"), lw=0.9)
    ax.set(xlabel="Retained steady-state time (ns)", ylabel="Voltage (V)",
           title="Only Tx4 on: four-nanosecond receiver snapshot", xlim=(times[0] * 1e9, times[0] * 1e9 + 4))
    ax.legend(ncol=4)
    fig.tight_layout()
    fig.savefig(figures / "tx4_only_snapshot.png")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4))
    lag = (np.arange(len(kernel.samples)) + kernel.first_lag) * kernel.dt
    for rx in range(8):
        for tx in range(8):
            ax.plot(lag / 1e-9, 20 * np.log10(np.maximum(np.abs(kernel.samples[:, rx, tx]), 1e-15)), alpha=.2, lw=.5)
    ax.set(xlabel="Signed lag (ns)", ylabel="Discrete coefficient magnitude (dB)",
           title="All 64 loaded paths: retained finite-bandwidth impulse coefficients", ylim=(-160, 0))
    fig.tight_layout()
    fig.savefig(figures / "impulse_all_paths.png")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    base = Path(__file__).resolve().parent
    parser.add_argument("--input", type=Path, default=base / "instructions" / "HW11_crosstalk_8lines.s16p")
    parser.add_argument("--output", type=Path, default=base / "results")
    parser.add_argument("--bits", type=int, default=2 ** 18, help="Retained random bits per transmitter (default 262144)")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--samples-per-ui", type=int, default=50)
    parser.add_argument("--block-size", type=int, default=32768)
    parser.add_argument("--kaiser-a", type=float, default=1.0)
    parser.add_argument("--pre-ns", type=float, default=5.0)
    parser.add_argument("--post-ns", type=float, default=20.0)
    args = parser.parse_args(argv)
    if args.bits < 2048 or args.samples_per_ui < 10 or args.samples_per_ui % 10:
        parser.error("Use at least 2048 bits and samples-per-ui divisible by 10")
    started = time.perf_counter()
    args.output.mkdir(parents=True, exist_ok=True)
    print("Reading 16-port channel and solving all 64 loaded transfer paths...", flush=True)
    frequency, s, z0 = load_touchstone(args.input)
    if s.shape[1] != 16:
        raise ValueError("HW11 requires the supplied sixteen-port network")
    h = solve_transfer(frequency, s, z0)
    tb, tr, vpp = 100e-12, 20e-12, 1.0
    npui, dt = args.samples_per_ui, tb / args.samples_per_ui
    delays = np.asarray([0, .5, 0, .5, 0, .5, 0, .5]) * tb
    delay_samples = np.rint(delays / dt).astype(int)
    rise_samples = round(tr / dt)
    kernel = frequency_to_impulse(frequency, h, dt, args.kaiser_a, args.pre_ns * 1e-9, args.post_ns * 1e-9)
    # Guard at least 100 ns at both ends; every retained sample sees real random
    # pre/post history, including the larger-support convergence reference.
    guard_ui = max(1000, int(np.ceil((args.post_ns + args.pre_ns) * 1e-9 / tb)) + 2)
    bits = independent_bits(args.bits + 2 * guard_ui + 2, 8, args.seed)
    print("Checking all-path impulse support convergence...", flush=True)
    convergence = _convergence_check(frequency, h, kernel, bits, npui, rise_samples, delay_samples, guard_ui, args.block_size, args.kaiser_a)
    sampling = _sampling_check(frequency, h, kernel, bits, npui, rise_samples, delay_samples, guard_ui, args.block_size, args.kaiser_a)
    eyes, extras = {}, {}
    for name, receivers, transmitters in [
        ("all_on", list(range(8)), list(range(8))),
        ("only_tx0", [0], [0]),
        ("only_tx4", [1, 2, 3, 4], [4]),
    ]:
        print(f"Simulating {name}: {args.bits} retained bits...", flush=True)
        eye, times, waveform, extra = simulate_case(bits, kernel, npui, rise_samples, delay_samples,
                                                   args.bits, guard_ui, receivers, transmitters, args.block_size)
        eyes[name], extras[name] = eye, extra
        if name == "only_tx4":
            snapshot_times, snapshot = times, waveform
    threshold = np.diag(h[0].real) / 2
    case_metrics = {
        "all_on": {f"rx{i}": metric for i, metric in enumerate(eyes["all_on"].metrics(threshold))},
        "only_tx0": {"rx0": eyes["only_tx0"].metrics(threshold[0])[0]},
        "only_tx4": {f"rx{i}": metric for i, metric in zip([1, 2, 3, 4], eyes["only_tx4"].metrics(threshold[[1, 2, 3, 4]]))},
    }
    _plots(args.output, frequency, s, eyes, snapshot_times, snapshot, kernel)
    np.savez_compressed(args.output / "channel.npz", frequency_hz=frequency, s=s, h=h,
                        impulse_coefficients=kernel.samples.astype(np.float32), first_lag=kernel.first_lag, dt_s=dt)
    np.savez_compressed(args.output / "simulation.npz", bits=bits, snapshot_times_s=snapshot_times,
                        tx4_only_snapshot_v=snapshot, voltage_edges_v=eyes["all_on"].edges,
                        all_on_counts=eyes["all_on"].counts, only_tx0_counts=eyes["only_tx0"].counts,
                        only_tx4_counts=eyes["only_tx4"].counts, only_tx4_receivers=np.array([1, 2, 3, 4]))
    metrics = {
        "configuration": {"bits": args.bits, "samples_per_ui": npui, "dt_s": dt, "Tb_s": tb,
                          "Tr_s": tr, "delays_s": delays.tolist(), "seed": args.seed, "vpp_v": vpp,
                          "rs_ohm": 30, "cs_f": 1.5e-12, "rl_ohm": 65, "cl_f": 1e-12,
                          "kaiser_a": args.kaiser_a, "pre_time_s": args.pre_ns * 1e-9,
                          "post_time_s": args.post_ns * 1e-9, "guard_ui": guard_ui,
                          "block_size": args.block_size, "simulation_precision": "float32 FFT and accumulation; diagnostics float64",
                          "ffe_taps": [1.0], "signaling": "single-ended NRZ"},
        "input": {"sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  "frequency_count": len(frequency), "df_hz": float(frequency[1] - frequency[0]),
                  "fmax_hz": float(frequency[-1]), "z0_ohm": z0,
                  "maximum_reciprocity_difference": float(np.abs(s - s.transpose(0, 2, 1)).max())},
        "transfer_dc_by_path": h[0].real.tolist(), "impulse_diagnostics": kernel.diagnostics,
        "convergence": convergence, "sampling_convergence": sampling, "cases": case_metrics, "waveform_statistics": extras,
        "eye_metric_definition": "Largest across phases of upper-cluster 0.1th percentile minus lower-cluster 99.9th percentile; split at half main-path DC gain. Each cluster requires >=20% samples. Empirical finite-sample metric, not BER or worst-case.",
        "eye_voltage_bin_v": float(eyes["all_on"].edges[1] - eyes["all_on"].edges[0]),
        "eye_probability_definition": "Probability mass per voltage bin, independently normalized at each phase.",
        "verification": {"ads_performed": False, "independent_tests": "Run python -m unittest discover -s homework/hw11/tests -v from repository root."},
        "reproducibility": {"command": f"python homework/hw11/analyze.py --bits {args.bits} --seed {args.seed} --samples-per-ui {npui} --block-size {args.block_size} --kaiser-a {args.kaiser_a:g} --pre-ns {args.pre_ns:g} --post-ns {args.post_ns:g}",
                            "analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                            "python_version": platform.python_version(), "numpy_version": np.__version__,
                            "scipy_version": version("scipy"), "matplotlib_version": version("matplotlib"),
                            "random_streams": "numpy.random.SeedSequence(seed).spawn(8), default_rng per channel"},
        "elapsed_seconds": time.perf_counter() - started,
    }
    (args.output / "metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Done in {metrics['elapsed_seconds']:.1f} s. Results: {args.output}", flush=True)
    return metrics


if __name__ == "__main__":
    main()
