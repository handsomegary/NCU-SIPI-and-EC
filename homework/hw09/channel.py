"""Discrete channel utilities adapted from the verified HW8 implementation.

Only channel/pulse/convolution utilities are reused; HW9 applies no FFE.
"""
from pathlib import Path
import numpy as np
from scipy.fft import irfft, rfft, next_fast_len
from scipy.special import i0


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
