"""Reproducible NCU EE8091 HW2: Kaiser-window step responses and TDR.

Run from HW2 with ``python src/hw2.py``. To recompute the independent
ngspice reference, add ``--ngspice /path/to/ngspice``. All paths default
to this homework directory and all generated outputs are anonymous.

Fourier convention: h(t) = integral H(f) exp(+j*2*pi*f*t) df.
The assignment's F is twice the one-sided cutoff: F = 2*f_max.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq
from scipy.signal import StateSpace, lsim


ROOT = Path(__file__).resolve().parents[1]
Z0 = 50.0
L = 4.314e-9
C = 0.5301e-12
ONE_WAY_DELAY = 250e-12
SOURCE_RAMP = 50e-12
TARGET_RISE = 40e-12
ALPHAS = np.arange(5, dtype=float)


@dataclass(frozen=True)
class Touchstone:
    """Touchstone 1.x S-parameter data, shaped [frequency, row, column]."""

    frequency_hz: np.ndarray
    s: np.ndarray
    reference_ohm: float


def read_touchstone(path: Path | str) -> Touchstone:
    """Read one- or two-port Touchstone 1.x RI, MA, or DB data.

    Honors frequency units, ignores comment tails, and permits wrapped
    numeric records. For two ports the standard order is S11,S21,S12,S22.
    Noise data and Touchstone 2.x keywords are deliberately rejected.
    """
    path = Path(path)
    match = re.fullmatch(r"\.s(\d+)p", path.suffix.lower())
    if match is None:
        raise ValueError("Expected a .sNp Touchstone filename")
    ports = int(match.group(1))
    if ports not in {1, 2}:
        raise ValueError("This homework reader supports one- and two-port files")
    # Touchstone defaults when an option line is absent.
    units, parameter, representation, impedance = "GHZ", "S", "MA", 50.0
    numbers: list[float] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        content = line.partition("!")[0].strip()
        if not content:
            continue
        if content.startswith("["):
            raise ValueError("This reader supports Touchstone 1.x, not 2.x")
        if content.startswith("#"):
            options = content[1:].upper().split()
            if len(options) < 3:
                raise ValueError("Incomplete Touchstone option line")
            units, parameter, representation = options[:3]
            if "R" in options:
                impedance = float(options[options.index("R") + 1])
            continue
        numbers.extend(float(x.replace("D", "E").replace("d", "e")) for x in content.split())
    scales = {"HZ": 1.0, "KHZ": 1e3, "MHZ": 1e6, "GHZ": 1e9}
    if units not in scales or parameter != "S" or representation not in {"RI", "MA", "DB"}:
        raise ValueError("Unsupported Touchstone units, parameter type, or representation")
    width = 1 + 2 * ports * ports
    if not numbers or len(numbers) % width:
        raise ValueError("Incomplete full-matrix Touchstone data record")
    raw = np.asarray(numbers).reshape(-1, width)
    frequency = raw[:, 0] * scales[units]
    if np.any(np.diff(frequency) <= 0) or np.any(frequency < 0):
        raise ValueError("Frequencies must be nonnegative and strictly increasing")
    pairs = raw[:, 1:].reshape(-1, ports * ports, 2)
    first, second = pairs[..., 0], pairs[..., 1]
    if representation == "RI":
        values = first + 1j * second
    else:
        magnitude = 10 ** (first / 20) if representation == "DB" else first
        values = magnitude * np.exp(1j * np.deg2rad(second))
    # Column-major flattening of a matrix: S11,S21,S12,S22 for N=2.
    matrices = values.reshape(-1, ports, ports).transpose(0, 2, 1)
    return Touchstone(frequency, matrices, impedance)


def kaiser_spectrum(frequency_hz: np.ndarray, F: float, alpha: float) -> np.ndarray:
    """Continuous frequency window, equal to zero outside |f| <= F/2."""
    if F <= 0 or alpha < 0:
        raise ValueError("F must be positive and alpha nonnegative")
    x = 2 * np.asarray(frequency_hz, dtype=float) / F
    inside = np.abs(x) <= 1
    result = np.zeros_like(x)
    result[inside] = np.i0(np.pi * alpha * np.sqrt(np.maximum(0, 1 - x[inside] ** 2)))
    return result / np.i0(np.pi * alpha)


def kaiser_impulse(time_s: np.ndarray, F: float, alpha: float) -> np.ndarray:
    """Analytic inverse Fourier transform of the continuous Kaiser window.

    np.sinc(x)=sin(pi*x)/(pi*x); therefore the square root contains
    alpha**2, not (pi*alpha)**2. The hyperbolic branch is evaluated
    explicitly so the result is real without complex cancellation.
    """
    q = F * np.asarray(time_s, dtype=float)
    square = q * q - alpha * alpha
    out = np.ones_like(q)
    positive, negative = square > 1e-14, square < -1e-14
    out[positive] = np.sinc(np.sqrt(square[positive]))
    u = np.pi * np.sqrt(-square[negative])
    out[negative] = np.sinh(u) / u
    return F * out / np.i0(np.pi * alpha)


_GL_X, _GL_W = leggauss(256)
_QUAD_X = (_GL_X + 1) / 2
_QUAD_W = _GL_W / 2


def kaiser_step(q: np.ndarray | float, alpha: float) -> np.ndarray:
    """Zero-centered normalized Kaiser step U(q), q=F*t.

    U(q)=1/2 + integral_0^1 W(x)*sin(pi*q*x)/(pi*x) dx.
    This form includes the infinite past analytically, avoiding the
    baseline error from integrating only a short time interval.
    """
    q = np.asarray(q, dtype=float)
    flat = q.reshape(-1)
    window = np.i0(np.pi * alpha * np.sqrt(1 - _QUAD_X ** 2)) / np.i0(np.pi * alpha)
    result = np.empty_like(flat)
    for begin in range(0, flat.size, 4096):
        value = flat[begin:begin + 4096, None]
        integrand = value * np.sinc(value * _QUAD_X) * window
        result[begin:begin + 4096] = 0.5 + integrand @ _QUAD_W
    return result.reshape(q.shape)


def normalized_rise(alpha: float) -> float:
    """Return q90-q10 for the central 10%-90% transition."""
    positive_crossing = brentq(lambda q: float(kaiser_step(q, alpha)) - 0.9, 0, 20, xtol=1e-13)
    return 2 * positive_crossing  # U(-q)=1-U(q).


def alpha_for_rise(F: float, rise_s: float) -> float:
    target = F * rise_s
    if target < normalized_rise(0):
        raise ValueError("Requested rise is faster than the alpha=0 window")
    return brentq(lambda a: normalized_rise(a) - target, 0, 20, xtol=1e-12)


def sparameters_from_circuit(frequency_hz: np.ndarray) -> np.ndarray:
    """Independent ABCD solution of the assigned line + L-C-L ladder."""
    frequency_hz = np.asarray(frequency_hz)
    omega = 2 * np.pi * frequency_hz
    n = frequency_hz.size
    series = np.zeros((n, 2, 2), complex)
    series[:, 0, 0] = series[:, 1, 1] = 1
    series[:, 0, 1] = 1j * omega * L
    shunt = np.zeros_like(series)
    shunt[:, 0, 0] = shunt[:, 1, 1] = 1
    shunt[:, 1, 0] = 1j * omega * C
    ladder = series @ shunt @ series
    A, B = ladder[:, 0, 0], ladder[:, 0, 1]
    cc, D = ladder[:, 1, 0], ladder[:, 1, 1]
    denominator = A + B / Z0 + cc * Z0 + D
    s11 = (A + B / Z0 - cc * Z0 - D) / denominator
    s21 = 2 / denominator
    result = np.zeros_like(ladder)
    result[:, 0, 0] = s11 * np.exp(-2j * omega * ONE_WAY_DELAY)
    result[:, 1, 0] = result[:, 0, 1] = s21 * np.exp(-1j * omega * ONE_WAY_DELAY)
    result[:, 1, 1] = s11
    return result


def windowed_reflection(
    frequency_hz: np.ndarray,
    s11: np.ndarray,
    alpha: float,
    delay_s: float = SOURCE_RAMP / 2,
    n_fft: int = 2 ** 19,
    source: str = "kaiser",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Reconstruct continuous h and its step integral on an exact FFT grid.

    DC appears once; each nonzero positive-frequency sample has one
    conjugate negative-frequency partner. The last measured endpoint is
    half-weighted for the frequency-domain trapezoidal rule. Zero padding
    refines dt but does not extend bandwidth or the period 1/df.
    ``source='linear_ramp'`` substitutes sinc(f*T) for the Kaiser window.
    """
    frequency_hz = np.asarray(frequency_hz, dtype=float)
    s11 = np.asarray(s11, dtype=complex)
    df = frequency_hz[1] - frequency_hz[0]
    if frequency_hz[0] != 0 or not np.allclose(np.diff(frequency_hz), df, rtol=1e-10, atol=1e-4):
        raise ValueError("IFFT requires a uniform frequency grid starting at DC")
    if n_fft % 2 or n_fft <= 2 * (frequency_hz.size - 1):
        raise ValueError("Use an even FFT size greater than twice the highest measured index")
    if abs(s11[0].imag) > 1e-12:
        raise ValueError("Real impulse responses require a real DC point")
    F = 2 * frequency_hz[-1]
    if source == "kaiser":
        window = kaiser_spectrum(frequency_hz, F, alpha)
    elif source == "linear_ramp":
        window = np.sinc(frequency_hz * SOURCE_RAMP)
    else:
        raise ValueError("Unknown source shape")
    positive = s11 * window * np.exp(-2j * np.pi * frequency_hz * delay_s)
    positive[-1] *= 0.5
    full = np.zeros(n_fft, dtype=complex)
    full[:positive.size] = positive
    full[-positive.size + 1:] = positive[1:][::-1].conj()
    raw = np.fft.ifft(full) * n_fft * df
    impulse = np.fft.fftshift(raw.real)
    dt = 1 / (n_fft * df)
    time = (np.arange(n_fft) - n_fft // 2) * dt
    step = cumulative_trapezoid(impulse, dx=dt, initial=0)
    diagnostic = {
        "df_hz": float(df), "fmax_hz": float(frequency_hz[-1]),
        "F_hz": float(F), "n_fft": n_fft, "dt_s": float(dt),
        "period_s": float(1 / df),
        "dc_s11_real": float(s11[0].real),
        "max_ifft_imaginary": float(np.max(np.abs(raw.imag))),
        "max_ifft_imaginary_relative": float(np.max(np.abs(raw.imag)) / max(np.max(np.abs(raw.real)), 1)),
        "integrated_impulse": float(np.sum(impulse) * dt),
        "initial_step": float(step[0]), "final_step": float(step[-1]),
    }
    return time, impulse, step, diagnostic


def circuit_ramp_reflection(time_s: np.ndarray) -> np.ndarray:
    """Independent time-domain state-space reference, in volts.

    At the ladder input the source is 2*a in series with 50 ohm.
    States are [i_L1, v_C, i_L3]; reflected voltage is a-50*i_L1.
    Moving the reference plane through the matched line delays it by 2*TD.
    """
    time_s = np.asarray(time_s)
    t = time_s - time_s[0]
    incident = np.clip((time_s - 2 * ONE_WAY_DELAY) / SOURCE_RAMP, 0, 1)
    A = np.array([[-Z0 / L, -1 / L, 0], [1 / C, 0, -1 / C], [0, 1 / L, -Z0 / L]])
    B = np.array([[2 / L], [0], [0]])
    system = StateSpace(A, B, np.array([[-Z0, 0, 0]]), np.array([[1.0]]))
    _, reflected, _ = lsim(system, U=incident, T=t, interp=True)
    return reflected


def circuit_kaiser_reflection(time_s: np.ndarray, F: float, alpha: float) -> np.ndarray:
    """Independent differential-circuit response to the same Kaiser edge.

    Integrate the analytic impulse on a 0.2 ps grid with 20 ns prehistory;
    center the incident edge at U(0)=0.5. This calculation uses neither the
    supplied S-parameters nor an FFT. Voltage-scaled inductor states and
    nanosecond units keep the three circuit equations well scaled.
    """
    time_s = np.asarray(time_s)
    step_ns = 0.0002
    end_ns = max(4.0, float(np.max(time_s)) * 1e9 + 1)
    samples = int(np.ceil((end_ns + 20) / step_ns))
    dense_ns = -20 + np.arange(samples + 1) * step_ns
    impulse_per_ns = kaiser_impulse(dense_ns * 1e-9, F, alpha) * 1e-9
    incident = cumulative_trapezoid(impulse_per_ns, dx=step_ns, initial=0)
    incident += .5 - np.interp(0, dense_ns, incident)
    k = Z0 / L * 1e-9
    k2 = 1 / (Z0 * C) * 1e-9
    A = np.array([[-k, -k, 0], [k2, 0, -k2], [0, k, -k]])
    B = np.array([[2 * k], [0], [0]])
    system = StateSpace(A, B, np.array([[-1, 0, 0]]), np.array([[1.0]]))
    _, reflected, _ = lsim(system, U=incident, T=dense_ns - dense_ns[0], interp=True)
    delay_ns = (2 * ONE_WAY_DELAY + SOURCE_RAMP / 2) * 1e9
    return np.interp(time_s * 1e9 - delay_ns, dense_ns, reflected)


def run_ngspice(executable: str, simulation: Path, scratch_parent: Path | None = None) -> np.ndarray:
    """Run the saved anonymous standalone circuit, returning time and volts."""
    with tempfile.TemporaryDirectory(prefix="hw02-ngspice-", dir=scratch_parent or simulation) as directory:
        scratch = Path(directory)
        shutil.copy2(simulation / "chebyshev_reflection.cir", scratch)
        completed = subprocess.run(
            [executable, "-b", "chebyshev_reflection.cir"], cwd=scratch,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            encoding="utf-8", errors="replace", timeout=120,
        )
        # Output can include machine-specific install paths. Do not persist it.
        result = scratch / "spice_reference.txt"
        if completed.returncode or not result.exists():
            raise RuntimeError("ngspice failed; inspect the standalone circuit and your local installation")
        values = np.loadtxt(result)
    if values.ndim != 2 or values.shape[1] < 4:
        raise RuntimeError("Unexpected ngspice wrdata columns")
    # First column is time; the saved vectors are source, input port, reflection.
    return values[:, [0, 1, 2, 3]]


def error_metrics(observed: np.ndarray, reference: np.ndarray) -> dict:
    residual = np.asarray(observed) - np.asarray(reference)
    return {"rmse_v": float(np.sqrt(np.mean(residual ** 2))),
            "max_abs_error_v": float(np.max(np.abs(residual)))}


def save_csv(path: Path, columns: list[np.ndarray], names: list[str]) -> None:
    np.savetxt(path, np.column_stack(columns), delimiter=",", header=",".join(names), comments="", fmt="%.12e")


def make_figures(
    results: Path, data: Touchstone, rows: list[dict],
    time: np.ndarray, kaiser_response: np.ndarray, ramp_response: np.ndarray,
    spice: np.ndarray, state: np.ndarray, alpha: float,
    alpha_responses: dict[float, np.ndarray],
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25,
                         "figure.dpi": 150, "savefig.dpi": 200})
    destination = results / "figures"
    destination.mkdir(exist_ok=True)
    def save(figure, name):
        figure.tight_layout()
        figure.savefig(destination / (name + ".png"))
        figure.savefig(destination / (name + ".pdf"), metadata={"Author": "", "Creator": "", "Title": name})
        plt.close(figure)
    q = np.linspace(-5, 5, 4001)
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.5, 3.6))
    x = np.linspace(-1.15, 1.15, 2001)
    for a in ALPHAS:
        left.plot(x, kaiser_spectrum(x / 2, 1, a), label=rf"$\alpha={a:.0f}$")
        right.plot(q, kaiser_impulse(q, 1, a))
    left.set(xlabel=r"Normalized frequency $2f/F$", ylabel=r"$W(f)$", xlim=(-1.15, 1.15))
    left.legend(fontsize=8)
    right.set(xlabel=r"Normalized time $Ft$", ylabel=r"Normalized impulse $w(t)/F$", xlim=(-5, 5))
    save(fig, "part1_kaiser_window_impulse")
    fig, ax = plt.subplots(figsize=(7.0, 4.1))
    for a in ALPHAS:
        ax.plot(q, kaiser_step(q, a), label=rf"$\alpha={a:.0f}$")
    ax.axhline(.1, color="0.65", lw=.7, ls=":")
    ax.axhline(.9, color="0.65", lw=.7, ls=":")
    ax.set(xlabel=r"Normalized time $Ft$", ylabel=r"Kaiser step $u_W(t)$", ylim=(-.15, 1.15), xlim=(-5, 5))
    ax.legend(ncol=3, fontsize=9)
    save(fig, "part1_kaiser_steps")
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.plot([r["alpha"] for r in rows], [r["normalized_rise"] for r in rows], "o-")
    ax.set(xlabel=r"Kaiser parameter $\alpha$", ylabel=r"Normalized rise time $F t_{10-90}$", xticks=ALPHAS)
    save(fig, "part1_rise_time")
    fig, ax = plt.subplots(figsize=(7.0, 4.1))
    for row, col, label in [(0, 0, "$S_{11}$"), (1, 0, "$S_{21}$")]:
        magnitude = np.abs(data.s[:, row, col])
        ax.plot(data.frequency_hz / 1e9, 20 * np.log10(np.maximum(magnitude, 1e-12)), label=label)
    ax.set(xlabel="Frequency (GHz)", ylabel="Magnitude (dB)", xlim=(0, 20), ylim=(-80, 3))
    ax.legend()
    save(fig, "part2_sparameters")
    fig, (ax, bottom) = plt.subplots(2, 1, figsize=(7.1, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(time * 1e9, kaiser_response, label=rf"Kaiser IDFT ($\alpha={alpha:.6f}$)")
    ax.plot(time * 1e9, spice, "--", label="ngspice: input - source/2")
    ax.plot(time * 1e9, ramp_response, ":", label="Linear-ramp IDFT", lw=1.6)
    ax.set(ylabel="Reflected voltage (V)", xlim=(0, 3))
    ax.legend(fontsize=8)
    bottom.plot(time * 1e9, kaiser_response - spice, label="Kaiser - ngspice")
    bottom.plot(time * 1e9, ramp_response - spice, label="Ramp - ngspice")
    bottom.set(xlabel="Time (ns)", ylabel="Error (V)")
    bottom.legend(fontsize=8)
    save(fig, "part2_reflection_comparison")
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(time * 1e9, kaiser_response, label="Kaiser IDFT")
    ax.plot(time * 1e9, spice, "--", label="ngspice")
    ax.set(xlabel="Time (ns)", ylabel="Reflected voltage (V)", xlim=(.45, .75))
    ax.legend()
    save(fig, "part2_reflection_edge")
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for a, ax in zip([0.0, 4.0], axes):
        ax.plot(time * 1e9, alpha_responses[a], label=rf"Kaiser $\alpha={a:.0f}$")
        ax.plot(time * 1e9, kaiser_response, ":", label=rf"Matched $\alpha={alpha:.6f}$")
        ax.plot(time * 1e9, spice, "--", label="ngspice")
        ax.set(xlabel="Time (ns)", xlim=(0.35, 1.2))
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Reflected voltage (V)")
    save(fig, "part2_alpha_comparison")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ngspice", help="Executable for a fresh independent circuit simulation")
    parser.add_argument("--output", type=Path, help="Alternative output directory; default is HW2/results")
    parser.add_argument("--n-fft", type=int, default=2 ** 19)
    args = parser.parse_args(argv)
    results = args.output.resolve() if args.output else ROOT / "results"
    simulation = ROOT / "simulation"
    results.mkdir(parents=True, exist_ok=True)
    data = read_touchstone(ROOT / "data" / "chebyshev_lpf.s2p")
    frequency = data.frequency_hz
    F = 2 * frequency[-1]
    alpha = alpha_for_rise(F, TARGET_RISE)
    rows = [{"alpha": float(a), "normalized_rise": normalized_rise(a),
             "rise_ps_at_F_40GHz": normalized_rise(a) / F * 1e12} for a in ALPHAS]
    save_csv(results / "part1_rise_times.csv",
             [np.array([row[key] for row in rows]) for key in rows[0]], list(rows[0]))
    q = np.linspace(-5, 5, 4001)
    save_csv(results / "part1_steps.csv", [q] + [kaiser_step(q, a) for a in ALPHAS],
             ["normalized_time_Ft"] + [f"step_alpha_{int(a)}" for a in ALPHAS])
    t_fft, impulse, step, grid = windowed_reflection(frequency, data.s[:, 0, 0], alpha, n_fft=args.n_fft)
    _, _, ramp_step, _ = windowed_reflection(frequency, data.s[:, 0, 0], alpha, n_fft=args.n_fft, source="linear_ramp")
    _, _, coarse_step, _ = windowed_reflection(frequency, data.s[:, 0, 0], alpha, n_fft=args.n_fft // 2)
    t_coarse = (np.arange(args.n_fft // 2) - args.n_fft // 4) / (args.n_fft // 2 * grid["df_hz"])
    # A common 0.5 ps grid makes reference comparisons and CSV reuse simple.
    time = np.linspace(0, 3e-9, 6001)
    kaiser_response = np.interp(time, t_fft, step)
    ramp_response = np.interp(time, t_fft, ramp_step)
    alpha_responses = {}
    for a in (0.0, 4.0):
        _, _, response, _ = windowed_reflection(frequency, data.s[:, 0, 0], a, n_fft=args.n_fft)
        alpha_responses[a] = np.interp(time, t_fft, response)
    state = circuit_ramp_reflection(time)
    same_kaiser_state = circuit_kaiser_reflection(time, F, alpha)
    executable = args.ngspice or shutil.which("ngspice") or shutil.which("ngspice_con")
    if executable:
        raw_spice = run_ngspice(executable, simulation, results)
        reference_output = results / "spice_reference.csv" if args.output else simulation / "spice_reference.csv"
        save_csv(reference_output, [raw_spice[:, i] for i in range(4)],
                 ["time_s", "source_v", "input_port_v", "reflected_v"])
        spice_status = "recomputed"
    else:
        reference = simulation / "spice_reference.csv"
        if not reference.exists():
            raise RuntimeError("ngspice is unavailable and no saved independent reference exists")
        raw_spice = np.loadtxt(reference, delimiter=",", skiprows=1)
        spice_status = "saved_reference"
    spice = np.interp(time, raw_spice[:, 0], raw_spice[:, 3])
    save_csv(results / "part2_reflection.csv", [time, kaiser_response, ramp_response, spice, state, same_kaiser_state],
             ["time_s", "kaiser_idft_v", "linear_ramp_idft_v", "ngspice_reflected_v", "state_space_reflected_v", "same_kaiser_state_space_v"])
    modeled = sparameters_from_circuit(frequency)
    s11_error = float(np.max(np.abs(modeled[:, 0, 0] - data.s[:, 0, 0])))
    s21_error = np.abs(modeled[:, 1, 0] - data.s[:, 1, 0])
    anomaly = s21_error > 1e-4
    coarse = np.interp(time, t_coarse, coarse_step)
    convergence = error_metrics(kaiser_response, coarse)
    checks = {
        "kaiser_impulse_even": bool(np.allclose(kaiser_impulse(q, 1, alpha), kaiser_impulse(-q, 1, alpha), atol=1e-13)),
        "kaiser_step_symmetry": bool(np.allclose(kaiser_step(q, alpha) + kaiser_step(-q, alpha), 1, atol=1e-12)),
        "target_rise_within_0_001ps": bool(abs(normalized_rise(alpha) / F - TARGET_RISE) < 1e-15),
        "real_ifft": bool(grid["max_ifft_imaginary_relative"] < 1e-14),
        "dc_area_matches_s11": bool(abs(grid["integrated_impulse"] - data.s[0, 0, 0].real) < 1e-12),
        "fft_convergence_below_0_1mV": bool(convergence["max_abs_error_v"] < 1e-4),
        "ngspice_vs_state_below_0_1mV": bool(error_metrics(spice, state)["max_abs_error_v"] < 1e-4),
        "supplied_s11_matches_circuit_below_1ppm": bool(s11_error < 1e-6),
        "same_kaiser_state_vs_idft_below_0_1mV": bool(error_metrics(kaiser_response, same_kaiser_state)["max_abs_error_v"] < 1e-4),
    }
    metrics = {
        "course": "NCU EE8091", "homework": "HW2", "all_checks_passed": all(checks.values()),
        "part1": rows,
        "part2": {"alpha_for_40ps": float(alpha), "target_rise_ps": 40,
                  "actual_kaiser_rise_ps": float(normalized_rise(alpha) / F * 1e12),
                  "source_full_ramp_ps": SOURCE_RAMP * 1e12,
                  "source_ramp_10_90_ps": SOURCE_RAMP * .8 * 1e12,
                  "source_center_delay_ps": SOURCE_RAMP / 2 * 1e12,
                  "line_one_way_delay_ps": ONE_WAY_DELAY * 1e12,
                  "reflection_round_trip_delay_ps": 2 * ONE_WAY_DELAY * 1e12,
                  "spice_reference_status": spice_status,
                  "s11_circuit_max_complex_error": s11_error,
                  "s21_phase_anomaly_frequencies_ghz": (frequency[anomaly] / 1e9).tolist(),
                  "s21_other_points_max_complex_error": float(np.max(s21_error[~anomaly])),
                  "s21_max_magnitude_error": float(np.max(np.abs(np.abs(modeled[:, 1, 0]) - np.abs(data.s[:, 1, 0])))),
                  "kaiser_vs_ngspice": error_metrics(kaiser_response, spice),
                  "linear_ramp_idft_vs_ngspice": error_metrics(ramp_response, spice),
                  "state_space_vs_ngspice": error_metrics(state, spice),
                  "same_kaiser_state_vs_idft": error_metrics(same_kaiser_state, kaiser_response),
                  "same_kaiser_state_grid": {"dt_ps": .2, "prehistory_ns": 20, "uses_sparameters_or_fft": False},
                  "fft_size_doubled_convergence": convergence,
                  "alpha_comparison": [{"alpha": a, "rise_ps": normalized_rise(a) / F * 1e12,
                                        "peak_v": float(np.max(alpha_responses[a])),
                                        "minimum_v": float(np.min(alpha_responses[a])),
                                        "peak_time_ns": float(time[np.argmax(alpha_responses[a])] * 1e9),
                                        **error_metrics(alpha_responses[a], spice)} for a in (0.0, 4.0)],
                  "comparison_interval_ns": [0, 3], "comparison_samples": time.size,
                  "grid": grid},
        "checks": checks,
        "notes": [
            "Kaiser and piecewise-linear edges share the 40 ps 10-90 rise but differ in shape.",
            "Finite 20 GHz bandwidth and the window cause residual disagreement with broadband SPICE.",
            "The line TD is one way; input reflection accumulates twice that delay.",
            "Saved SPICE reference is generated independently from the standalone netlist.",
            "The supplied S11 agrees with the nominal circuit within export rounding.",
            "Supplied S21 has isolated 180-degree phase anomalies at 2,6,10,14,18 GHz; preserved as source data and not used in reflection.",
        ],
    }
    (results / "verification.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    save_csv(results / "part2_alpha_comparison.csv", [time, alpha_responses[0], kaiser_response, alpha_responses[4], spice],
             ["time_s", "kaiser_alpha0_v", "kaiser_matched_v", "kaiser_alpha4_v", "ngspice_reflected_v"])
    make_figures(results, data, rows, time, kaiser_response, ramp_response, spice, state, alpha, alpha_responses)
    print(json.dumps({"all_checks_passed": metrics["all_checks_passed"], "alpha": alpha,
                      "kaiser_vs_ngspice": metrics["part2"]["kaiser_vs_ngspice"],
                      "output_directory": "results"}, indent=2))
    if not metrics["all_checks_passed"]:
        raise SystemExit("One or more numerical checks failed; see results/verification.json")


if __name__ == "__main__":
    main()
