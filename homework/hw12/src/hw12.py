#!/usr/bin/env python3
"""HW12: bounded-memory jitter simulation and basic analysis.

Run from any working directory, for example::

    python homework/hw12/src/hw12.py --output homework/hw12/results

The default evaluates exactly 2**20 steady-state bits at 200 samples/UI.
Waveforms are generated and filtered in blocks; long waveforms are never saved.
The requested horizontal eye slices and an independent, interpolated crossing
histogram are both reported. An eye slice is voltage occupancy, so its weighting
can differ from the event-count weighting of the crossing histogram.

Jitter shifts the START of each linear transition; the simulation rounds
shifts to the nearest sample (ties to even), avoiding truncation toward zero.
The reusable waveform function also accepts continuous shifts. Physically
overlapping transitions cause an explicit error instead of silent corruption.
The discrete channel kernel uses the reference single-sided Kaiser convention
(a=0.5) and is truncated at 40 ns. Its continuous-time impulse response is the
discrete kernel divided by dt. All RMS figures are centered RMS (sigma).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
from scipy.fft import ifft, irfft, next_fast_len, rfft
from scipy.special import i0


UI = 25e-12
RISE = 5e-12
VPP = 1.0
SJ_PEAK = 3e-12
SJ_FREQUENCY = 1e6
RJ_SIGMA = 1e-12
DV = 1e-3
TRUNCATION = 40e-9
KAISER_A = 0.5
CASES = (
    ("case1_no_tx_jitter", "No Tx jitter", 0.0, 0.0),
    ("case2_sj", "SJ only", SJ_PEAK, 0.0),
    ("case3_rj", "RJ only", 0.0, RJ_SIGMA),
    ("case4_sj_rj", "SJ + RJ", SJ_PEAK, RJ_SIGMA),
)


def read_touchstone_2port(path: str | Path):
    """Return f[Hz], S[f, receiving_port, incident_port], and real Z0.

    Touchstone 1.x two-port data are ordered S11,S21,S12,S22. MA, DB,
    and RI representations and Hz/kHz/MHz/GHz units are supported. In-line
    comments and wrapped data lines are accepted. No reciprocal assumption
    or transpose is applied.
    """
    unit, representation, z0 = "ghz", "ma", 50.0
    tokens: list[float] = []
    with Path(path).open(encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.split("!", 1)[0].strip()
            if not line:
                continue
            if line.startswith("#"):
                options = line[1:].lower().split()
                if not options or options[0] not in {"hz", "khz", "mhz", "ghz"}:
                    raise ValueError("Unsupported Touchstone frequency unit")
                unit = options[0]
                if "s" not in options:
                    raise ValueError("This reader requires S parameters")
                representation = next((x for x in options if x in {"ma", "db", "ri"}), "ma")
                if "r" in options:
                    z0 = float(options[options.index("r") + 1])
            elif line.startswith("["):
                raise ValueError("Use the supplied Touchstone 1.x two-port file")
            else:
                tokens.extend(float(x.replace("D", "E").replace("d", "e")) for x in line.split())
    if not tokens or len(tokens) % 9:
        raise ValueError("Incomplete two-port Touchstone record (expected 9 numbers)")
    data = np.asarray(tokens, dtype=float).reshape(-1, 9)
    f = data[:, 0] * {"hz": 1.0, "khz": 1e3, "mhz": 1e6, "ghz": 1e9}[unit]
    pairs = data[:, 1:].reshape(-1, 4, 2)
    if representation == "ri":
        values = pairs[:, :, 0] + 1j * pairs[:, :, 1]
    else:
        magnitude = pairs[:, :, 0] if representation == "ma" else 10 ** (pairs[:, :, 0] / 20)
        values = magnitude * np.exp(1j * np.deg2rad(pairs[:, :, 1]))
    s = np.empty((f.size, 2, 2), dtype=complex)
    s[:, 0, 0], s[:, 1, 0] = values[:, 0], values[:, 1]
    s[:, 0, 1], s[:, 1, 1] = values[:, 2], values[:, 3]
    if z0 <= 0 or not np.all(np.isfinite(s)) or not np.all(np.isfinite(f)):
        raise ValueError("Nonfinite data or invalid reference impedance")
    if f.size < 2 or f[0] != 0 or np.any(np.diff(f) <= 0):
        raise ValueError("Frequency data must start at DC and increase")
    if not np.allclose(np.diff(f), f[1] - f[0], rtol=1e-9, atol=1e-3):
        raise ValueError("Uniform frequency spacing is required")
    return f, s, z0


def channel_transfer(f, s, z0, rs=20.0, cs=0.3e-12, rl=85.0, cl=0.3e-12):
    """Loaded transfer Vload / Videal_source using a 4-by-4 MNA solve.

    Currents point into the two-port: (I-S)V - Z0(I+S)Iport = 0.
    The source is one volt; Vpp is applied only during waveform generation.
    """
    f, s = np.asarray(f), np.asarray(s)
    if s.shape != (f.size, 2, 2):
        raise ValueError("Expected S shape (frequency_count, 2, 2)")
    if min(rs, rl, z0) <= 0 or min(cs, cl) < 0:
        raise ValueError("Invalid circuit component")
    identity = np.eye(2)
    a = np.zeros((f.size, 4, 4), dtype=complex)
    a[:, 0, 0] = 1 / rs + 2j * np.pi * f * cs
    a[:, 1, 1] = 1 / rl + 2j * np.pi * f * cl
    a[:, :2, 2:] = identity
    a[:, 2:, :2] = identity - s
    a[:, 2:, 2:] = -z0 * (identity + s)
    rhs = np.zeros((f.size, 4, 1), dtype=complex)
    rhs[:, 0, 0] = 1 / rs
    return np.linalg.solve(a, rhs)[:, 1, 0]


def build_fir(f, h_frequency, dt, truncation=TRUNCATION, kaiser_a=KAISER_A):
    """Build a causal discrete convolution kernel, with correct DC/Nyquist.

    Frequencies above the supplied maximum are zero-filled. The sampled
    impulse has no extra 1/dt scaling: it directly convolves voltage samples.
    """
    f, h_frequency = np.asarray(f, float), np.asarray(h_frequency, complex)
    if dt <= 0 or truncation <= 0 or kaiser_a < 0 or f.size < 2:
        raise ValueError("Invalid FIR grid or truncation")
    if f.shape != h_frequency.shape or f[0] != 0:
        raise ValueError("H must have the same DC-inclusive grid as f")
    df = f[1] - f[0]
    if df <= 0 or not np.allclose(np.diff(f), df, rtol=1e-9, atol=1e-3):
        raise ValueError("Uniform positive df is required")
    m = int(round(1 / (dt * df)))
    if not np.isclose(m * dt * df, 1, rtol=1e-9):
        raise ValueError("dt and df must define an integer DFT period")
    if f[-1] > 0.5 / dt * (1 + 1e-12):
        raise ValueError("Input frequency exceeds the output Nyquist frequency")
    ntaps = int(round(truncation / dt))
    if ntaps < 1 or ntaps >= m:
        raise ValueError("Truncation must lie strictly inside the DFT period")
    w = i0(np.pi * kaiser_a * np.sqrt(np.maximum(0, 1 - (f / f[-1]) ** 2))) / i0(np.pi * kaiser_a)
    w[0] *= 0.5
    if m % 2 == 0 and f.size == m // 2 + 1:
        w[-1] *= 0.5
    full = 2 * np.real(ifft(h_frequency * w, n=m))
    return np.ascontiguousarray(full[:ntaps])


def generate_waveform_segment(bits, jitter, dt, ui, rise, start_sample, stop_sample,
                              bit_start=0, vpp=VPP):
    """Sample a piecewise-linear unipolar NRZ segment on the global n*dt grid.

    bits[k] starts its transition at (bit_start+k)*ui+jitter[k]. Before the
    first supplied bit the level is zero. Supply guard bits on both sides of
    a requested segment when using negative jitter or nonzero initial levels.
    """
    bits, jitter = np.asarray(bits), np.asarray(jitter, dtype=float)
    if bits.ndim != 1 or bits.size != jitter.size or not bits.size:
        raise ValueError("bits and jitter must be nonempty equal-length vectors")
    if dt <= 0 or ui <= 0 or rise <= 0 or rise > ui or stop_sample < start_sample:
        raise ValueError("Invalid waveform grid")
    if not np.all((bits == 0) | (bits == 1)) or not np.all(np.isfinite(jitter)):
        raise ValueError("Bits must be binary and jitter finite")
    previous = np.r_[0, bits[:-1]]
    transition = np.flatnonzero(bits != previous)
    if transition.size == 0:
        return np.zeros(stop_sample - start_sample, dtype=float)
    starts = (bit_start + transition) * (ui / dt) + jitter[transition] / dt
    # Exact grid-aligned ramps remain exact even at very large sample indices.
    if np.allclose(starts, np.rint(starts), rtol=0, atol=1e-6):
        starts = np.rint(starts)
    rise_samples = rise / dt
    if abs(rise_samples - round(rise_samples)) < 1e-10:
        rise_samples = float(round(rise_samples))
    if np.any(np.diff(starts) < rise_samples * (1 - 1e-10)):
        raise ValueError("Jitter causes physically overlapping transition ramps")
    t = np.arange(start_sample, stop_sample, dtype=float)
    event = np.searchsorted(starts, t, side="right") - 1
    result = np.zeros(t.size, dtype=float)
    valid = event >= 0
    selected = event[valid]
    old = previous[transition[selected]].astype(float)
    delta = bits[transition[selected]].astype(float) - old
    fraction = np.clip((t[valid] - starts[selected]) / rise_samples, 0, 1)
    result[valid] = vpp * (old + delta * fraction)
    return result


class _FIRFilter:
    def __init__(self, h):
        self.h = np.asarray(h, dtype=float)
        if self.h.ndim != 1 or not self.h.size:
            raise ValueError("FIR kernel must be nonempty")
        self.history = np.zeros(self.h.size - 1)
        self.spectra = {}

    def process(self, block):
        block = np.asarray(block, dtype=float)
        if not block.size:
            return block.copy()
        combined = np.concatenate((self.history, block))
        size = next_fast_len(combined.size, real=True)
        if size not in self.spectra:
            self.spectra[size] = rfft(self.h, size)
        transformed = irfft(rfft(combined, size) * self.spectra[size], size)
        offset = self.h.size - 1
        result = transformed[offset:offset + block.size].copy()
        if offset:
            self.history = combined[-offset:].copy()
        return result


def stream_fir(blocks: Iterable[np.ndarray], h) -> Iterator[np.ndarray]:
    """Causal FFT overlap-save, equivalent to np.convolve(x,h)[:len(x)]."""
    filt = _FIRFilter(h)
    for block in blocks:
        yield filt.process(block)


def accumulate_eye(eye, waveform, start_sample, points_per_ui, vmin, dv):
    """Add raw counts to one UI; caller controls voltage bounds."""
    waveform = np.asarray(waveform)
    if eye.ndim != 2 or eye.shape[1] != points_per_ui or dv <= 0:
        raise ValueError("Invalid eye shape or voltage grid")
    voltage_index = np.floor((waveform - vmin) / dv).astype(np.int64)
    if np.any(voltage_index < 0) or np.any(voltage_index >= eye.shape[0]):
        raise ValueError("Waveform outside eye voltage bounds")
    phase = (np.arange(waveform.size, dtype=np.int64) + int(start_sample)) % points_per_ui
    counts = np.bincount(voltage_index * points_per_ui + phase, minlength=eye.size)
    eye += counts.reshape(eye.shape)


class _EyeAccumulator:
    """Grow only the voltage grid, keeping a fixed dv-aligned global origin."""
    def __init__(self, points_per_ui, dv):
        self.p = points_per_ui
        self.dv = dv
        self.lo = -200
        self.counts = np.zeros((1601, points_per_ui), dtype=np.int64)

    @property
    def vmin(self):
        return self.lo * self.dv

    def add(self, waveform, start_sample):
        if not waveform.size:
            return
        required_lo = math.floor(float(np.min(waveform)) / self.dv) - 1
        required_hi = math.floor(float(np.max(waveform)) / self.dv) + 2
        old_hi = self.lo + self.counts.shape[0]
        lo, hi = min(self.lo, required_lo), max(old_hi, required_hi)
        if lo != self.lo or hi != old_hi:
            counts = np.zeros((hi - lo, self.p), dtype=np.int64)
            counts[self.lo - lo:old_hi - lo] = self.counts
            self.counts, self.lo = counts, lo
        accumulate_eye(self.counts, waveform, start_sample, self.p, self.vmin, self.dv)


def _distribution(raw, dt, reference_index=None):
    raw = np.asarray(raw, dtype=float)
    if raw.ndim != 1 or raw.size == 0 or not np.all(np.isfinite(raw)) or np.any(raw < 0):
        raise ValueError("Invalid distribution counts")
    total = float(raw.sum())
    if total == 0:
        raise ValueError("No observations in the selected threshold row")
    reference_index = int(np.argmax(raw)) if reference_index is None else int(reference_index)
    reference_index %= raw.size
    centered = np.roll(raw, raw.size // 2 - reference_index)
    pmf = centered / total
    t = (np.arange(raw.size) - raw.size // 2) * dt
    mean = float(np.sum(t * pmf))
    variance = float(np.sum((t - mean) ** 2 * pmf))
    return {"raw_slice": raw.copy(), "pmf": pmf, "pdf": pmf / dt,
            "time_seconds": t, "reference_index": reference_index,
            "mean_seconds": mean, "variance_seconds2": variance,
            "rms_seconds": math.sqrt(variance), "observation_count": total}


def analyze_eye_slice(eye, vmin, dv, threshold, dt, reference_index=None):
    """Normalize a one-UI horizontal eye row, preserving exact zero cells.

    Each time column has the same number of evaluated bits, so raw counts
    and vertically normalized eye probabilities give the same final PMF.
    """
    eye = np.asarray(eye)
    row = int(math.floor((float(np.real(threshold)) - vmin) / dv))
    if eye.ndim != 2 or not 0 <= row < eye.shape[0]:
        raise ValueError("Threshold is outside the eye")
    result = _distribution(eye[row], dt, reference_index)
    result["threshold_volts"] = float(np.real(threshold))
    result["voltage_bin_index"] = row
    result["voltage_bin_edges"] = [vmin + row * dv, vmin + (row + 1) * dv]
    return result


class _PatternBlocks:
    """Separate, chunk-invariant random streams with two guard bits/side."""
    def __init__(self, seed):
        self.bit_rng = np.random.default_rng(np.random.SeedSequence([seed, 11]))
        self.rj_rng = np.random.default_rng(np.random.SeedSequence([seed, 29]))
        self.bits = self.bit_rng.integers(0, 2, size=4, dtype=np.int64)
        self.normals = self.rj_rng.normal(size=4)
        self.position = 0

    def take(self, count):
        self.bits = np.r_[self.bits, self.bit_rng.integers(0, 2, size=count, dtype=np.int64)]
        self.normals = np.r_[self.normals, self.rj_rng.normal(size=count)]
        indices = np.arange(self.position - 2, self.position + count + 2)
        bits, normals = self.bits, self.normals
        self.bits, self.normals = bits[count:].copy(), normals[count:].copy()
        self.position += count
        return indices, bits, normals


class _Crossings:
    def __init__(self, threshold, p, dt, valid_start, valid_stop):
        self.threshold, self.p, self.dt = threshold, p, dt
        self.valid_start, self.valid_stop = valid_start, valid_stop
        self.counts = np.zeros(p, dtype=np.int64)
        self.previous = None

    def add(self, waveform, start_sample):
        if not waveform.size:
            return
        if self.previous is None:
            values = waveform
            first = start_sample
        else:
            values = np.r_[self.previous, waveform]
            first = start_sample - 1
        left, right = values[:-1], values[1:]
        crossing = np.flatnonzero(((left < self.threshold) & (right >= self.threshold)) |
                                  ((left > self.threshold) & (right <= self.threshold)))
        fraction = (self.threshold - left[crossing]) / (right[crossing] - left[crossing])
        sample_position = first + crossing + fraction
        keep = (sample_position >= self.valid_start) & (sample_position < self.valid_stop)
        phase = np.remainder(sample_position[keep], self.p)
        rounded = np.rint(phase)
        phase = np.where(np.abs(phase - rounded) < 1e-6, rounded, phase)
        bins = np.floor(phase).astype(int) % self.p
        self.counts += np.bincount(bins, minlength=self.p)
        self.previous = float(waveform[-1])


def _json_stats(result):
    return {key: value for key, value in result.items()
            if key not in {"raw_slice", "pmf", "pdf", "time_seconds"}}


def _simulate_case(bits, p, block_bits, seed, h, sj_amplitude, rj_sigma, threshold_rx, warmup_bits):
    dt = UI / p
    start_valid = warmup_bits * p
    stop_valid = (warmup_bits + bits) * p
    pattern, filt = _PatternBlocks(seed), _FIRFilter(h)
    eyes = {"tx": _EyeAccumulator(p, DV), "rx": _EyeAccumulator(p, DV)}
    crossings = {"tx": _Crossings(VPP / 2, p, dt, start_valid, stop_valid),
                 "rx": _Crossings(threshold_rx, p, dt, start_valid, stop_valid)}
    excerpt = {"tx": [], "rx": []}
    jitter_sums = {"requested": np.zeros(3), "applied": np.zeros(3)}
    excerpt_limit = start_valid + 40 * p
    total_bits = warmup_bits + bits
    elapsed_start, last_message = time.monotonic(), time.monotonic()
    for first_bit in range(0, total_bits, block_bits):
        count = min(block_bits, total_bits - first_bit)
        indices, block_pattern, normals = pattern.take(count)
        requested_jitter = sj_amplitude * np.sin(2 * np.pi * SJ_FREQUENCY * indices * UI) + rj_sigma * normals
        jitter = np.rint(requested_jitter / dt) * dt
        if np.max(np.abs(jitter)) >= UI:
            raise ValueError("Jitter exceeds the waveform generator guard interval")
        first_sample, last_sample = first_bit * p, (first_bit + count) * p
        x = generate_waveform_segment(block_pattern, jitter, dt, UI, RISE,
                                      first_sample, last_sample, int(indices[0]), VPP)
        y = filt.process(x)
        central_indices = indices[2:count + 2]
        actual_edges = ((block_pattern[2:count + 2] != block_pattern[1:count + 1]) &
                        (central_indices >= warmup_bits) & (central_indices < warmup_bits + bits))
        for label, values in (("requested", requested_jitter), ("applied", jitter)):
            selected = values[2:count + 2][actual_edges]
            jitter_sums[label] += [selected.size, selected.sum(), np.sum(selected**2)]
        lo, hi = max(start_valid, first_sample), min(stop_valid, last_sample)
        for side, waveform in (("tx", x), ("rx", y)):
            crossings[side].add(waveform, first_sample)
            if hi > lo:
                eyes[side].add(waveform[lo - first_sample:hi - first_sample], lo)
            take_hi = min(hi, excerpt_limit)
            if take_hi > lo:
                excerpt[side].append(waveform[lo - first_sample:take_hi - first_sample].copy())
        if time.monotonic() - last_message > 20:
            print(f"  {(first_bit + count) / total_bits:6.1%} ({time.monotonic()-elapsed_start:.1f} s)", flush=True)
            last_message = time.monotonic()
    jitter_stats = {}
    for label, (count, sum_x, sum_x2) in jitter_sums.items():
        mean = sum_x / count if count else 0.0
        variance = max(0.0, sum_x2 / count - mean**2) if count else 0.0
        jitter_stats[label] = {"observation_count": int(count), "mean_seconds": float(mean),
                               "variance_seconds2": float(variance), "rms_seconds": math.sqrt(variance)}
    return (eyes, crossings,
            {key: np.concatenate(value) if value else np.array([]) for key, value in excerpt.items()},
            jitter_stats)


def _save_plots(output, f, s, h_frequency, h, dt, results, eyes, crossings, excerpts, parameters):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ["#3366aa", "#df8f24", "#32956f", "#a54e9b"]
    half_ui_ps = parameters["ui_seconds"] / 2e-12
    plt.rcParams.update({"font.size": 13, "legend.fontsize": 11,
                         "axes.grid": True, "grid.alpha": 0.2})
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    axs[0].plot(f / 1e9, 20 * np.log10(np.maximum(np.abs(s[:, 1, 0]), 1e-16)), label="S21")
    axs[0].plot(f / 1e9, 20 * np.log10(np.maximum(np.abs(h_frequency), 1e-16)), label="Loaded H")
    axs[0].set(xlabel="Frequency (GHz)", ylabel="Magnitude (dB)", title="Channel and loaded circuit")
    axs[0].legend()
    axs[1].plot(np.arange(h.size) * dt / 1e-9, h / dt)
    axs[1].set(xlabel="Time (ns)", ylabel="Impulse response (1/s)", title="Kaiser a=0.5; retained 40 ns")
    fig.savefig(output / "channel_diagnostics.png", dpi=180)
    plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for ax, (key, name, _, _) in zip(axs.flat, CASES):
        eye = eyes[key]["rx"]
        normalized = eye.counts / parameters["bits"]
        image = np.ma.masked_where(normalized == 0, normalized)
        im = ax.imshow(np.ma.log10(np.ma.hstack((image, image))), origin="lower", aspect="auto",
                       extent=[0, 2, eye.vmin, eye.vmin + eye.counts.shape[0] * DV],
                       vmin=-6, vmax=-1, cmap="viridis", interpolation="nearest")
        ax.set(xlabel="Time (UI)", ylabel="Voltage (V)", title=name,
               ylim=(-0.15, max(0.95, parameters["threshold_rx_volts"] * 2.1)))
        ax.grid(False)
    fig.colorbar(im, ax=axs, label="log10 probability per voltage bin")
    fig.savefig(output / "rx_eyes.png", dpi=180)
    plt.close(fig)
    for quantity, filename, ylabel in (("pmf", "tx_rx_pmf.png", "Probability per time bin"),
                                       ("pdf", "tx_rx_pdf.png", "Density (1/ps)")):
        fig, axs = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
        for side, ax in zip(("tx", "rx"), axs):
            for color, (key, name, _, _) in zip(colors, CASES):
                dist = results[key][side]
                values = dist[quantity] * (1e-12 if quantity == "pdf" else 1)
                ax.semilogy(dist["time_seconds"] / 1e-12, np.ma.masked_equal(values, 0),
                            label=name, color=color, markersize=6,
                            marker="o" if np.count_nonzero(values) == 1 else None)
            ax.set(xlabel="Time relative to baseline reference (ps)", ylabel=ylabel,
                   title=f"{side.upper()} horizontal eye slice", xlim=(-half_ui_ps, half_ui_ps))
            ax.legend()
        fig.savefig(output / filename, dpi=180)
        plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for side, ax in zip(("tx", "rx"), axs):
        for color, (key, name, _, _) in zip(colors, CASES):
            dist = crossings[key][side]
            ax.semilogy(dist["time_seconds"] / 1e-12, np.ma.masked_equal(dist["pmf"], 0),
                        label=name, color=color, markersize=6,
                        marker="o" if np.count_nonzero(dist["pmf"]) == 1 else None)
        ax.set(xlabel="Time relative to baseline crossing reference (ps)", ylabel="Crossing probability per bin",
               title=f"{side.upper()} interpolated crossing events", xlim=(-half_ui_ps, half_ui_ps))
        ax.legend()
    fig.savefig(output / "crossing_diagnostics.png", dpi=180)
    plt.close(fig)
    if excerpts is None:
        return
    fig, axs = plt.subplots(4, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
    for ax, color, (key, name, _, _) in zip(axs, colors, CASES):
        waveform = excerpts[key]["rx"]
        ax.plot(np.arange(waveform.size) * dt / 1e-12, waveform, color=color, lw=0.8)
        ax.set(ylabel="Rx voltage (V)", title=name)
    axs[-1].set_xlabel("Time from evaluation start (ps)")
    fig.savefig(output / "waveform_excerpt.png", dpi=180)
    plt.close(fig)


def run(channel, output, bits=2**20, points_per_ui=200, block_bits=4096, seed=0, plots=True):
    """Run all four cases and write reproducible numerical/figure artifacts."""
    if min(bits, points_per_ui, block_bits) < 1 or seed < 0:
        raise ValueError("Positive bit/grid/block counts and nonnegative seed are required")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    f, s, z0 = read_touchstone_2port(channel)
    h_frequency = channel_transfer(f, s, z0)
    if abs(h_frequency[0].imag) > 1e-10 or h_frequency[0].real <= 0:
        raise ValueError("Expected positive real DC loaded gain")
    dt = UI / points_per_ui
    h = build_fir(f, h_frequency, dt)
    warmup_bits = math.ceil((h.size - 1) / points_per_ui) + 2
    threshold_rx = VPP * float(h_frequency[0].real) / 2
    parameters = {"bits": bits, "points_per_ui": points_per_ui, "block_bits": block_bits,
                  "seed": seed, "ui_seconds": UI, "rise_seconds": RISE, "dt_seconds": dt,
                  "dv_volts": DV, "vpp_volts": VPP, "sj_peak_seconds": SJ_PEAK,
                  "sj_frequency_hz": SJ_FREQUENCY, "rj_sigma_seconds": RJ_SIGMA,
                  "rs_ohms": 20.0, "rl_ohms": 85.0, "cs_farads": 0.3e-12, "cl_farads": 0.3e-12,
                  "kaiser_a": KAISER_A, "truncation_seconds": TRUNCATION,
                  "warmup_bits": warmup_bits, "threshold_tx_volts": VPP / 2,
                  "threshold_rx_volts": threshold_rx, "evaluated_sj_cycles": bits * UI * SJ_FREQUENCY,
                  "jitter_model": "nearest-sample displacement of linear-ramp start; np.rint ties-to-even; no clipping",
                  "quantization_error_bound_seconds": dt / 2,
                  "guard_bits_each_side": 2, "requested_bits": bits, "analyzed_bits": bits,
                  "rng": "PCG64, independent SeedSequence([seed,11]) bits and [seed,29] normal jitter",
                  "evaluation": "exact bit count after FIR-length warmup; no padded zero tail",
                  "rms_definition": "standard deviation (RMS after subtracting mean)"}
    channel_info = {"file": Path(channel).name, "frequency_samples": int(f.size),
                    "df_hz": float(f[1]-f[0]), "maximum_frequency_hz": float(f[-1]),
                    "reference_impedance_ohms": float(z0), "dc_gain": float(h_frequency[0].real),
                    "fir_taps": int(h.size), "fir_dc_gain": float(h.sum()),
                    "fir_dc_relative_error": float((h.sum()-h_frequency[0].real)/h_frequency[0].real),
                    "fir_last_1ns_l1_fraction": float(np.sum(abs(h[-max(1,round(1e-9/dt)):])) / np.sum(abs(h))),
                    "s_parameter_orientation": "S[receiving_port,incident_port]; no transpose"}
    results, all_eyes, all_crossings, excerpts = {}, {}, {}, {}
    eye_references, crossing_references = {}, {}
    checks = {"exact_evaluated_bits": True, "finite_outputs": True, "pmfs_sum_to_one": True,
              "pdfs_integrate_to_one": True, "zero_probabilities_preserved": True}
    for key, name, sj_amplitude, rj_sigma in CASES:
        print(f"Simulating {name}: {bits:,} evaluated bits, {points_per_ui} points/UI", flush=True)
        eyes, crossings, excerpt, jitter_stats = _simulate_case(bits, points_per_ui, block_bits, seed, h,
                                                  sj_amplitude, rj_sigma, threshold_rx, warmup_bits)
        case_result, case_crossing = {"name": name, "input_jitter": jitter_stats}, {}
        save = {}
        for side, threshold in (("tx", VPP/2), ("rx", threshold_rx)):
            eye = eyes[side]
            if not np.all(eye.counts.sum(axis=0) == bits):
                raise AssertionError("Eye columns do not contain the exact evaluated bit count")
            result = analyze_eye_slice(eye.counts, eye.vmin, DV, threshold, dt, eye_references.get(side))
            crossing_result = _distribution(crossings[side].counts, dt, crossing_references.get(side))
            eye_references.setdefault(side, result["reference_index"])
            crossing_references.setdefault(side, crossing_result["reference_index"])
            case_result[side] = result
            case_result[f"crossing_{side}"] = _json_stats(crossing_result)
            case_crossing[side] = crossing_result
            save[f"{side}_counts"] = eye.counts
            save[f"{side}_vmin"] = np.array(eye.vmin)
            save[f"{side}_crossing_counts"] = crossings[side].counts
            for dist in (result, crossing_result):
                if not np.isclose(dist["pmf"].sum(), 1, rtol=0, atol=1e-12):
                    raise AssertionError("PMF normalization failed")
                if not np.isclose(np.sum(dist["pdf"])*dt, 1, rtol=0, atol=1e-12):
                    raise AssertionError("PDF normalization failed")
        np.savez_compressed(output / f"{key}_eyes.npz", **save, dv=np.array(DV), dt=np.array(dt),
                            points_per_ui=np.array(points_per_ui), evaluated_bits=np.array(bits))
        results[key], all_eyes[key], all_crossings[key], excerpts[key] = case_result, eyes, case_crossing, excerpt
        print(f"  eye-slice sigma Tx={case_result['tx']['rms_seconds']/1e-12:.6f} ps; "
              f"Rx={case_result['rx']['rms_seconds']/1e-12:.6f} ps", flush=True)
    for filename, source in (("jitter_distributions.csv", results), ("crossing_distributions.csv", all_crossings)):
        with (output / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["case", "side", "time_seconds", "pmf", "pdf_per_second"])
            for key, _, _, _ in CASES:
                for side in ("tx", "rx"):
                    dist = source[key][side]
                    writer.writerows((key, side, float(t), float(p), float(d))
                                     for t, p, d in zip(dist["time_seconds"], dist["pmf"], dist["pdf"]))
    summary_cases = {}
    for key, name, sj_amplitude, rj_sigma in CASES:
        summary_cases[key] = {"name": name, "tx": _json_stats(results[key]["tx"]),
                              "rx": _json_stats(results[key]["rx"]),
                              "crossing_tx": results[key]["crossing_tx"],
                              "crossing_rx": results[key]["crossing_rx"],
                              "input_jitter": results[key]["input_jitter"],
                              "sj_peak_seconds": sj_amplitude, "rj_sigma_seconds": rj_sigma,
                              "ideal_tx_sigma_seconds": math.sqrt(sj_amplitude**2/2+rj_sigma**2)}
    summary = {"parameters": parameters, "channel": channel_info, "checks": checks, "cases": summary_cases,
               "notes": ["Horizontal eye slices are voltage occupancy, not exact crossing-event distributions.",
                         "Interpolated crossing diagnostics count actual rising and falling threshold events equally.",
                         "Time references are circular baseline modes, held fixed across cases separately at Tx/Rx.",
                         "A one-UI histogram cannot represent unbounded Gaussian tails without folding; current sigma is far below half UI.",
                         "Ideal SJ sigma assumes full phase coverage. Finite runs, especially smoke tests, need not match it.",
                         "Jitter is rounded to nearest sample, rather than truncated toward zero; continuum theory is a reference, not an exact finite-grid prediction.",
                         "Zero cells are preserved; plotting masks them without adding artificial probability.",
                         "Variance is centered. Fixed channel propagation delay is excluded from centered RMS.",
                         "The source assignment Q4 prints Tx4 on the RHS of an Rx relation; report treats this as a typo."]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    if plots:
        _save_plots(output, f, s, h_frequency, h, dt, results, all_eyes, all_crossings, excerpts, parameters)
    return summary


def render_existing(channel, output):
    """Re-render numerical artifacts without rerunning the simulation.

    The original waveform excerpt figure is left in place because full
    waveforms are intentionally not persisted. All statistical/channel plots
    are reconstructed from the saved integer counts and CSV distributions.
    """
    output = Path(output)
    summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
    parameters = summary["parameters"]
    dt = parameters["dt_seconds"]
    f, s, z0 = read_touchstone_2port(channel)
    h_frequency = channel_transfer(f, s, z0)
    h = build_fir(f, h_frequency, dt, parameters["truncation_seconds"], parameters["kaiser_a"])
    eyes = {}
    for key, _, _, _ in CASES:
        eyes[key] = {}
        with np.load(output / f"{key}_eyes.npz") as saved:
            for side in ("tx", "rx"):
                eye = _EyeAccumulator(parameters["points_per_ui"], parameters["dv_volts"])
                eye.counts = saved[f"{side}_counts"].copy()
                eye.lo = round(float(saved[f"{side}_vmin"]) / eye.dv)
                eyes[key][side] = eye
    distributions = []
    for filename in ("jitter_distributions.csv", "crossing_distributions.csv"):
        grouped = {key: {side: [] for side in ("tx", "rx")} for key, _, _, _ in CASES}
        with (output / filename).open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                grouped[row["case"]][row["side"]].append(
                    [float(row["time_seconds"]), float(row["pmf"]), float(row["pdf_per_second"])])
        mapped = {}
        for key, _, _, _ in CASES:
            mapped[key] = {}
            for side in ("tx", "rx"):
                data = np.asarray(grouped[key][side])
                if data.shape != (parameters["points_per_ui"], 3):
                    raise ValueError(f"Unexpected distribution length in {filename}: {key}/{side}")
                mapped[key][side] = {"time_seconds": data[:, 0], "pmf": data[:, 1], "pdf": data[:, 2]}
        distributions.append(mapped)
    _save_plots(output, f, s, h_frequency, h, dt, distributions[0], eyes, distributions[1], None, parameters)
    print(f"Re-rendered five statistical/channel figures in {output}; original waveform excerpt retained.", flush=True)


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--channel", type=Path, default=root / "instructions" / "HW12_channel.s2p")
    parser.add_argument("--output", type=Path, default=root / "results")
    parser.add_argument("--bits", type=int, default=2**20)
    parser.add_argument("--points-per-ui", type=int, default=200)
    parser.add_argument("--block-bits", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-plots", action="store_true", help="Write numerical artifacts without rendering figures")
    parser.add_argument("--render-only", action="store_true", help="Re-render existing results without resimulating")
    args = parser.parse_args(argv)
    if args.render_only:
        render_existing(args.channel, args.output)
    else:
        run(args.channel, args.output, args.bits, args.points_per_ui, args.block_bits, args.seed, not args.no_plots)


if __name__ == "__main__":
    main()
