"""Reproducible NRZ/PAM4 comparison through a measured differential channel.

The differential source is a Thevenin source with 1 V total open-circuit
peak-to-peak swing; matched 50-ohm single-ended ports give H = Sdd21/2.
FFE operates on symbols before DAC transitions are placed on a jittered clock.
No transmitted symbols are used by the decision-feedback slicer.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.fft import irfft
from scipy.optimize import linprog
from scipy.signal import fftconvolve


SOURCE_VPP = 1.0
RISE_S = 10e-12
SJ_AMPLITUDE_S = 1.5e-12
SJ_FREQUENCY_HZ = 10e6
RJ_SIGMA_S = 1e-12
CTLE_PEAK_DB = 18.0
FFE_COUNT = 3
DFE_COUNT = 5
DESIGN_JITTER_S = 3e-12
DEFAULT_SEED = 2025


def read_touchstone(path):
    """Read row-major Touchstone 1.x four-port data, including MA/RI/DB."""
    options = ['ghz', 's', 'ma', 'r', '50']
    numbers = []
    for raw in Path(path).read_text(encoding='utf-8-sig').splitlines():
        line = raw.split('!', 1)[0].strip()
        if not line:
            continue
        if line.startswith('#'):
            options = line[1:].lower().split()
        elif line.startswith('['):
            raise ValueError('Only Touchstone 1.x four-port files are supported')
        else:
            numbers.extend(float(v.replace('D', 'E').replace('d', 'e')) for v in line.split())
    if not numbers or len(numbers) % 33:
        raise ValueError('Each four-port frequency point requires 33 numbers')
    if 's' not in options:
        raise ValueError('Scattering parameters are required')
    scales = {'hz': 1.0, 'khz': 1e3, 'mhz': 1e6, 'ghz': 1e9}
    unit = next((v for v in options if v in scales), 'ghz')
    fmt = next((v for v in options if v in ('ma', 'ri', 'db')), 'ma')
    z0 = float(options[options.index('r') + 1]) if 'r' in options else 50.0
    values = np.asarray(numbers, dtype=float).reshape(-1, 33)
    f = values[:, 0] * scales[unit]
    q = values[:, 1:].reshape(-1, 4, 4, 2)
    if fmt == 'ri':
        s = q[..., 0] + 1j * q[..., 1]
    else:
        magnitude = 10 ** (q[..., 0] / 20) if fmt == 'db' else q[..., 0]
        s = magnitude * np.exp(1j * np.deg2rad(q[..., 1]))
    if len(f) < 2 or not np.isfinite(values).all() or f[0] != 0 or not np.all(np.diff(f) > 0):
        raise ValueError('Finite, increasing frequency points including DC are required')
    if not np.allclose(np.diff(f), f[1] - f[0], rtol=1e-8, atol=1e-3):
        raise ValueError('A uniformly sampled frequency grid is required')
    return f, s, z0


def differential_forward(s):
    """Forward differential transmission for left (1,3), right (2,4)."""
    s = np.asarray(s)
    return (s[..., 1, 0] - s[..., 1, 2] - s[..., 3, 0] + s[..., 3, 2]) / 2


def gray_encode(bits):
    """00,01,11,10 map to ascending normalized PAM4 voltage levels."""
    bits = np.asarray(bits)
    if bits.ndim != 1 or bits.size % 2 or not np.isin(bits, [0, 1]).all():
        raise ValueError('PAM4 requires an even-length one-dimensional binary record')
    code = 2 * bits[0::2] + bits[1::2]
    return np.array([-1.0, -1 / 3, 1.0, 1 / 3])[code]


def gray_decode(symbols):
    symbols = np.asarray(symbols, dtype=float)
    levels = np.array([-1.0, -1 / 3, 1 / 3, 1.0])
    indices = np.argmin(np.abs(symbols[..., None] - levels), axis=-1)
    return np.array([[0, 0], [0, 1], [1, 1], [1, 0]], dtype=np.uint8)[indices].reshape(-1)


def ctle_peak_frequency(zero_hz, pole1_hz, pole2_hz):
    """Analytic continuous-frequency maximum of a one-zero/two-pole CTLE."""
    z, p1, p2 = map(float, (zero_hz, pole1_hz, pole2_hz))
    if min(z, p1, p2) <= 0:
        raise ValueError('CTLE corner frequencies must be positive')
    b, c = (z / p1) ** 2, (z / p2) ** 2
    if b + c >= 1:
        return 0.0
    x = -1 + np.sqrt(1 + (1 - b - c) / (b * c))
    return float(z * np.sqrt(max(0.0, x)))


def ctle_response(f, zero_hz, pole1_hz, pole2_hz, peak_db=CTLE_PEAK_DB):
    """Normalize the absolute continuous CTLE gain to the common peak cap."""
    def raw(freq):
        return (1 + 1j * np.asarray(freq) / zero_hz) / ((1 + 1j * np.asarray(freq) / pole1_hz) * (1 + 1j * np.asarray(freq) / pole2_hz))
    peak_frequency = ctle_peak_frequency(zero_hz, pole1_hz, pole2_hz)
    return raw(f) * 10 ** (peak_db / 20) / abs(raw(peak_frequency))


def impulse_response(f, transfer, dt_s, duration_s):
    """Real Hermitian IDFT, zero above measured band, causal finite kernel.

    The kernel is dimensionless and used in discrete convolution without an
    extra dt factor. No DC renormalization or synthetic high-frequency data
    are introduced. The full 1/df record is retained for truncation diagnostics.
    """
    df = float(f[1] - f[0])
    n = round(1 / (df * dt_s))
    if not np.isclose(n * df * dt_s, 1.0, rtol=1e-8):
        raise ValueError('Requested symbol/time grid is not integral in 1/df')
    if f[-1] > 0.5 / dt_s or len(f) > n // 2 + 1:
        raise ValueError('Time-domain Nyquist frequency must exceed the measured band')
    count = round(duration_s / dt_s)
    if not 0 < count < n:
        raise ValueError('Impulse truncation must lie strictly inside the IDFT period')
    spectrum = np.zeros(n // 2 + 1, dtype=complex)
    spectrum[:len(f)] = transfer
    spectrum[0] = spectrum[0].real
    full = irfft(spectrum, n=n)
    return full[:count].copy(), full


def source_pulse(samples_per_ui, dt_s, rise_s=RISE_S, vpp=SOURCE_VPP):
    """Single normalized +1 symbol; ramp duration is explicitly 0%-100%."""
    ui = samples_per_ui * dt_s
    t = np.arange(int(np.ceil((ui + rise_s) / dt_s)) + 1) * dt_s
    return vpp / 2 * (np.clip(t / rise_s, 0, 1) - np.clip((t - ui) / rise_s, 0, 1))


def jittered_boundaries(count, ui_s, seed, bit_stride=1):
    """Use common underlying bit-clock RJ draws for NRZ and PAM4.

    PAM4 selects every second bit boundary. Thus both modulations see the same
    absolute SJ law, and exactly the same RJ draw at a common bit boundary.
    """
    if count < 1 or bit_stride not in (1, 2):
        raise ValueError('Positive count and bit stride 1 or 2 are required')
    nominal = np.arange(count, dtype=float) * ui_s
    rj = np.random.default_rng(seed).normal(0, RJ_SIGMA_S, count * bit_stride)[::bit_stride]
    sj = SJ_AMPLITUDE_S * np.sin(2 * np.pi * SJ_FREQUENCY_HZ * nominal)
    return nominal + rj + sj, rj, sj


def source_waveform(levels_v, boundary_times_s, dt_s, rise_s, stop_s=None):
    """Sample causal linear ramps beginning at each jittered clock boundary."""
    levels = np.asarray(levels_v, dtype=float)
    boundary = np.asarray(boundary_times_s, dtype=float)
    if levels.ndim != 1 or boundary.shape != levels.shape or not len(levels):
        raise ValueError('One finite voltage and clock boundary per symbol are required')
    if not np.isfinite(levels).all() or not np.isfinite(boundary).all() or min(dt_s, rise_s) <= 0:
        raise ValueError('Finite data and positive time steps/rise time are required')
    if np.any(np.diff(boundary) <= rise_s):
        raise ValueError('Jittered clock produces overlapping source ramps')
    stop = boundary[-1] + rise_s + dt_s if stop_s is None else float(stop_s)
    t = np.arange(max(1, int(np.ceil(stop / dt_s)))) * dt_s
    xp = np.column_stack((boundary, boundary + rise_s)).ravel()
    yp = np.column_stack((np.r_[0.0, levels[:-1]], levels)).ravel()
    return np.interp(t, xp, yp, left=0.0, right=levels[-1])


def apply_ffe_pulse(pulse, taps, samples_per_ui):
    result = np.zeros(len(pulse) + (len(taps) - 1) * samples_per_ui)
    for j, tap in enumerate(taps):
        result[j * samples_per_ui:j * samples_per_ui + len(pulse)] += tap * pulse
    return result


def cursor_matrix(pulse, clock, samples_per_ui, n_taps=FFE_COUNT, pad_samples=0):
    first = int(np.ceil((-clock - pad_samples) / samples_per_ui))
    last = int(np.floor((len(pulse) - 1 + (n_taps - 1) * samples_per_ui - clock + pad_samples) / samples_per_ui))
    k = np.arange(first, last + 1)
    index = clock + (k[:, None] - np.arange(n_taps)[None, :]) * samples_per_ui
    matrix = np.zeros(index.shape)
    good = (index >= 0) & (index < len(pulse))
    matrix[good] = pulse[index[good]]
    return k, matrix


def design_ffe(pulse, clock, samples_per_ui, modulation, dt_s):
    """Conservative worst-pattern LP with a common +/-3 ps clock stress.

    For M equally spaced levels, the residual half-eye is main/(M-1)-ISI.
    Central postcursors 1..5 define fixed DFE coefficients for design. At an
    offset phase their difference from the central coefficient remains ISI.
    Tiny rows are bounded by column L1 sums, not silently discarded. The
    three PAM4 eye thresholds share one clock and one central main cursor.
    """
    m = 2 if modulation == 'nrz' else 4
    phase_shift = max(1, round(DESIGN_JITTER_S / dt_s))
    k, matrix = cursor_matrix(pulse, clock, samples_per_ui, pad_samples=phase_shift)
    main = matrix[k == 0][0]
    # A fixed +/-3 ps clock displacement is a deterministic design stress,
    # not a probabilistic bound on Gaussian jitter. Its purpose is to avoid
    # choosing an extremely tall but narrow noiseless central eye.
    phases = (-phase_shift, 0, phase_shift)
    groups = []
    for offset in phases:
        index = clock + offset + (k[:, None] - np.arange(FFE_COUNT)[None, :]) * samples_per_ui
        shifted = np.zeros(index.shape)
        good = (index >= 0) & (index < len(pulse))
        shifted[good] = pulse[index[good]]
        remaining = shifted.copy()
        remaining[k == 0] = 0
        post = (k >= 1) & (k <= DFE_COUNT)
        remaining[post] -= matrix[post]
        active = np.max(np.abs(remaining), axis=1) >= max(np.max(np.abs(matrix)) * 1e-4, 1e-12)
        a = remaining[active]
        tail = np.abs(remaining[~active]).sum(axis=0)
        main_phase = shifted[k == 0][0]
        margins = [main_phase] if m == 2 else [main_phase / 3, main_phase - 2 * main / 3, 2 * main / 3 - main_phase / 3]
        groups.append((a, tail, margins))
    nr = sum(len(a) for a, _, _ in groups)
    n_margin = sum(len(margins) for _, _, margins in groups)
    n_var = 2 * FFE_COUNT + nr + 1
    objective = np.zeros(n_var); objective[-1] = -1
    constraints = np.zeros((2 * FFE_COUNT + 1 + 2 * nr + n_margin, n_var))
    upper = np.zeros(len(constraints))
    for j in range(FFE_COUNT):
        constraints[j, j] = 1
        constraints[j, FFE_COUNT + j] = -1
        constraints[FFE_COUNT + j, j] = -1
        constraints[FFE_COUNT + j, FFE_COUNT + j] = -1
    constraints[2 * FFE_COUNT, FFE_COUNT:2 * FFE_COUNT] = 1
    upper[2 * FFE_COUNT] = 1
    row = 2 * FFE_COUNT + 1
    residual_column = 2 * FFE_COUNT
    for a, tail, margins in groups:
        begin = residual_column
        for residual in a:
            constraints[row, :FFE_COUNT] = residual
            constraints[row + 1, :FFE_COUNT] = -residual
            constraints[row:row + 2, residual_column] = -1
            residual_column += 1; row += 2
        for margin in margins:
            constraints[row, :FFE_COUNT] = -margin
            constraints[row, FFE_COUNT:2 * FFE_COUNT] = tail
            constraints[row, begin:residual_column] = 1
            constraints[row, -1] = 1
            row += 1
    fit = linprog(objective, A_ub=constraints, b_ub=upper,
                  bounds=[(None, None)] * FFE_COUNT + [(0, None)] * (FFE_COUNT + nr) + [(None, None)], method='highs')
    if not fit.success:
        raise RuntimeError('FFE linear program did not converge: ' + fit.message)
    taps = fit.x[:FFE_COUNT]
    if np.abs(taps).sum() < 1e-8:
        taps = np.array([0.0, 1.0, 0.0])
    else:
        taps /= np.abs(taps).sum()
    cursors = matrix @ taps
    main_value = float(cursors[k == 0][0])
    interference = float(np.abs(cursors[(k < 0) | (k > DFE_COUNT)]).sum())
    score = float(fit.x[-1])
    return taps, float(score), main_value, interference


def design_equalizers(f, sdd21, rate_gbps, modulation, samples_per_ui=64, impulse_ns=40.0):
    stride = 1 if modulation == 'nrz' else 2
    ui = stride / (rate_gbps * 1e9)
    dt = ui / samples_per_ui
    pulse = source_pulse(samples_per_ui, dt)
    h_channel, full_channel = impulse_response(f, sdd21 / 2, dt, impulse_ns * 1e-9)
    raw_pulse = fftconvolve(h_channel, pulse)
    best = None
    grid = [(z, p, 40.0) for p in (10.0, 20.0, 40.0, 60.0)
            for z in (0.5, 1.0, 2.0, 4.0, 8.0) if z < p]
    for zero, pole1, pole2 in grid:
        response = ctle_response(f, zero * 1e9, pole1 * 1e9, pole2 * 1e9)
        h_eq, _ = impulse_response(f, sdd21 / 2 * response, dt, impulse_ns * 1e-9)
        single = fftconvolve(h_eq, pulse)
        middle = int(np.argmax(single)) + samples_per_ui
        phases = np.unique(np.rint(middle + np.linspace(-0.5, 0.5, 9) * samples_per_ui).astype(int))
        for clock in phases:
            taps, score, main, isi = design_ffe(single, clock, samples_per_ui, modulation, dt)
            if best is None or score > best['score']:
                best = dict(score=score, main=main, isi=isi, clock=int(clock), ffe=taps,
                            ctle=response, zero=zero, pole1=pole1, pole2=pole2,
                            h=h_eq, sbr=single)
    center = best['clock']
    for clock in range(center - samples_per_ui // 8, center + samples_per_ui // 8 + 1):
        taps, score, main, isi = design_ffe(best['sbr'], clock, samples_per_ui, modulation, dt)
        if score > best['score']:
            best.update(score=score, main=main, isi=isi, clock=clock, ffe=taps)
    full_eq = impulse_response(f, sdd21 / 2 * best['ctle'], dt, impulse_ns * 1e-9)[1]
    sbr_ffe = apply_ffe_pulse(best['sbr'], best['ffe'], samples_per_ui)
    dfe = sbr_ffe[best['clock'] + np.arange(1, DFE_COUNT + 1) * samples_per_ui]
    thresholds = np.array([0.0]) if modulation == 'nrz' else best['main'] * np.array([-2 / 3, 0.0, 2 / 3])
    peak_f = ctle_peak_frequency(best['zero'] * 1e9, best['pole1'] * 1e9, best['pole2'] * 1e9)
    expected_dc = float((sdd21[0] * best['ctle'][0] / 2).real)
    settings = {
        'ctle': {'zero_ghz': best['zero'], 'poles_ghz': [best['pole1'], best['pole2']],
                 'peak_gain_db': CTLE_PEAK_DB, 'peak_frequency_ghz': peak_f / 1e9,
                 'dc_gain_db': float(20 * np.log10(abs(best['ctle'][0])))},
        'ffe': {'taps': best['ffe'].tolist(), 'l1_norm': float(np.abs(best['ffe']).sum()), 'tap_count': FFE_COUNT},
        'dfe': {'taps_v': dfe.tolist(), 'thresholds_v': thresholds.tolist(),
                'main_cursor_v': best['main'], 'tap_count': DFE_COUNT},
    }
    design = {'method': 'Noiseless single-symbol response and conservative worst-pattern LP with fixed DFE coefficients and fixed thresholds at nominal clock and +/-3 ps stress; payload is not used for fitting',
              'ctle_grid_candidates': len(grid), 'clock_sample': int(best['clock']),
              'clock_time_ps': best['clock'] * dt * 1e12,
              'pda_eye_height_v': 2 * (best['main'] / (1 if modulation == 'nrz' else 3) - best['isi']),
              'robust_design_half_margin_v': best['score'], 'design_clock_stress_ps': DESIGN_JITTER_S * 1e12,
              'uncancelled_isi_sum_v': best['isi'],
              'channel_kernel_dc': float(h_channel.sum()), 'channel_expected_dc': float(sdd21[0].real / 2),
              'equalized_kernel_dc': float(best['h'].sum()), 'equalized_expected_dc': expected_dc,
              'equalized_dc_relative_error': float((best['h'].sum() - expected_dc) / expected_dc),
              'channel_discarded_l1': float(np.abs(full_channel[len(h_channel):]).sum()),
              'equalized_discarded_l1': float(np.abs(full_eq[len(best['h']):]).sum())}
    return {'f': f, 'sdd21': sdd21, 'ui': ui, 'dt': dt, 'modulation': modulation,
            'rate_gbps': rate_gbps, 'samples_per_ui': samples_per_ui, 'impulse_ns': impulse_ns,
            'h_channel': h_channel, 'h_eq': best['h'], 'source_pulse': pulse,
            'sbr_channel': raw_pulse, 'sbr_eq': best['sbr'], 'sbr_ffe': sbr_ffe,
            'ctle': best['ctle'], 'ffe': best['ffe'], 'dfe': dfe,
            'main': best['main'], 'thresholds': thresholds, 'clock': int(best['clock']),
            'equalizers': settings, 'design': design}


def fixed_design(f, sdd21, existing, samples_per_ui, impulse_ns):
    """Refine discretization without retuning physical settings or clock."""
    ui = existing['ui']; dt = ui / samples_per_ui
    h_channel, full_channel = impulse_response(f, sdd21 / 2, dt, impulse_ns * 1e-9)
    h_eq, full_eq = impulse_response(f, sdd21 / 2 * existing['ctle'], dt, impulse_ns * 1e-9)
    pulse = source_pulse(samples_per_ui, dt)
    sbr = fftconvolve(h_eq, pulse)
    sbr_ffe = apply_ffe_pulse(sbr, existing['ffe'], samples_per_ui)
    clock = round(existing['clock'] * existing['dt'] / dt)
    model = dict(existing)
    # Freeze actual analog DFE voltages and thresholds during a numerical check.
    model.update(samples_per_ui=samples_per_ui, impulse_ns=impulse_ns, dt=dt,
                 h_channel=h_channel, h_eq=h_eq, source_pulse=pulse,
                 sbr_channel=fftconvolve(h_channel, pulse), sbr_eq=sbr,
                 sbr_ffe=sbr_ffe, clock=clock)
    return model


def load_design(f, sdd21, case, samples_per_ui, impulse_ns):
    """Rebuild a saved physical design without tuning it on a new record."""
    modulation = case['modulation']; rate = float(case['bit_rate_gbps'])
    ui = (1 if modulation == 'nrz' else 2) / (rate * 1e9)
    dt = ui / samples_per_ui
    settings = case['equalizers']
    c = settings['ctle']
    response = ctle_response(f, c['zero_ghz'] * 1e9, c['poles_ghz'][0] * 1e9,
                             c['poles_ghz'][1] * 1e9, c['peak_gain_db'])
    h_channel, _ = impulse_response(f, sdd21 / 2, dt, impulse_ns * 1e-9)
    h_eq, _ = impulse_response(f, sdd21 / 2 * response, dt, impulse_ns * 1e-9)
    pulse = source_pulse(samples_per_ui, dt)
    single = fftconvolve(h_eq, pulse)
    taps = np.array(settings['ffe']['taps'])
    clock = round(case['design']['clock_time_ps'] * 1e-12 / dt)
    return {'f': f, 'sdd21': sdd21, 'ui': ui, 'dt': dt, 'modulation': modulation,
            'rate_gbps': rate, 'samples_per_ui': samples_per_ui, 'impulse_ns': impulse_ns,
            'h_channel': h_channel, 'h_eq': h_eq, 'source_pulse': pulse,
            'sbr_channel': fftconvolve(h_channel, pulse), 'sbr_eq': single,
            'sbr_ffe': apply_ffe_pulse(single, taps, samples_per_ui),
            'ctle': response, 'ffe': taps, 'dfe': np.array(settings['dfe']['taps_v']),
            'main': settings['dfe']['main_cursor_v'], 'thresholds': np.array(settings['dfe']['thresholds_v']),
            'clock': clock, 'equalizers': settings, 'design': {**case['design'],
                'evaluation_clock_sample': clock, 'evaluation_clock_time_ps': clock * dt * 1e12,
                'evaluation_impulse_ns': impulse_ns, 'evaluation_samples_per_ui': samples_per_ui,
                'loaded_frozen_design': True}}


def dfe_decisions(samples_v, taps_v, main_cursor_v, modulation):
    """Recursive decisions only; normalized outputs feed back through V taps."""
    if modulation not in ('nrz', 'pam4') or main_cursor_v <= 0:
        raise ValueError('Use NRZ/PAM4 and a positive main cursor')
    samples = np.asarray(samples_v, dtype=float)
    taps = np.asarray(taps_v, dtype=float)
    decisions = np.zeros(len(samples), dtype=float)
    feedback = np.zeros(len(samples), dtype=float)
    corrected = np.empty(len(samples), dtype=float)
    threshold = main_cursor_v * 2 / 3
    for n, sample in enumerate(samples):
        count = min(n, len(taps))
        if count:
            feedback[n] = np.dot(taps[:count], decisions[n - count:n][::-1])
        value = sample - feedback[n]
        corrected[n] = value
        if modulation == 'nrz':
            decisions[n] = 1.0 if value >= 0 else -1.0
        else:
            decisions[n] = -1.0 if value < -threshold else (-1 / 3 if value < 0 else (1 / 3 if value < threshold else 1.0))
    return decisions, feedback, corrected


def contiguous_width(phase_ps, margin, center):
    """Contiguous strict-positive margin around clock, interpolated crossings."""
    if margin[center] <= 0:
        return 0.0, float(phase_ps[center]), float(phase_ps[center])
    left = right = center
    while left > 0 and margin[left - 1] > 0:
        left -= 1
    while right < len(margin) - 1 and margin[right + 1] > 0:
        right += 1
    def crossing(a, b):
        return float(phase_ps[a] - margin[a] * (phase_ps[b] - phase_ps[a]) / (margin[b] - margin[a]))
    lo = crossing(left - 1, left) if left > 0 else float(phase_ps[left])
    hi = crossing(right, right + 1) if right < len(margin) - 1 else float(phase_ps[right])
    return hi - lo, lo, hi


def measure_eye(received, truth_symbols, clock_sample, samples_per_ui, dt_s,
                start, count, feedback_v=None, main_cursor_v=None, histogram=True):
    """Conditional finite-record envelopes and fixed-threshold eye openings.

    DFE feedback for symbol n is held on [clock_n-UI/2, clock_n+UI/2).
    Its right endpoint belongs to the next symbol's hold. These holds are
    causal because all feedback uses decisions at least one UI earlier.
    """
    p = samples_per_ui
    offsets = np.arange(-p // 2, p // 2 + 1)
    idx = np.arange(start, start + count)
    truth = np.asarray(truth_symbols)[idx]
    levels = np.unique(truth)
    if len(levels) not in (2, 4):
        raise ValueError('The retained record must contain every modulation level')
    if main_cursor_v is None:
        main_cursor_v = 1.0
    thresholds = main_cursor_v * (levels[:-1] + levels[1:]) / 2
    minima = np.empty((len(levels), len(offsets)))
    maxima = np.empty_like(minima)
    centers = clock_sample + idx * p
    span = float(np.max(np.abs(received))) + (float(np.max(np.abs(feedback_v))) if feedback_v is not None else 0)
    span = max(span * 1.02, 1e-3)
    edges = np.linspace(-span, span, 601)
    hist = np.zeros((600, len(offsets)), dtype=np.uint32) if histogram else None
    for column, offset in enumerate(offsets):
        voltage = np.asarray(received)[centers + offset].copy()
        if feedback_v is not None:
            voltage -= np.asarray(feedback_v)[idx + int(offset == p // 2)]
        for i, level in enumerate(levels):
            selected = voltage[truth == level]
            minima[i, column] = np.min(selected)
            maxima[i, column] = np.max(selected)
        if histogram:
            hist[:, column] = np.histogram(voltage, bins=edges)[0]
    phase = offsets * dt_s * 1e12
    center = p // 2
    heights = minima[1:, center] - maxima[:-1, center]
    widths, crossings = [], []
    for i, threshold in enumerate(thresholds):
        margin = np.minimum(minima[i + 1] - threshold, threshold - maxima[i])
        width, lo, hi = contiguous_width(phase, margin, center)
        widths.append(width); crossings.append([lo, hi])
    ui_ps = p * dt_s * 1e12
    metrics = {'eye_heights_v': heights.tolist(), 'eye_widths_ps': widths,
               'eye_widths_ui': (np.asarray(widths) / ui_ps).tolist(),
               'worst_eye_height_v': float(np.min(heights)),
               'worst_eye_width_ps': float(np.min(widths)),
               'worst_eye_width_ui': float(np.min(widths) / ui_ps),
               'crossings_ps': crossings, 'measured_symbols': int(count),
               'level_counts': [int(np.count_nonzero(truth == level)) for level in levels],
               'thresholds_v': thresholds.tolist()}
    arrays = {'phase_ps': phase, 'level_min_v': minima, 'level_max_v': maxima,
              'levels': levels, 'voltage_edges_v': edges}
    if histogram:
        arrays['histogram'] = hist
    return metrics, arrays


def simulate_record(model, bits=65536, seed=DEFAULT_SEED, histogram=True):
    p = model['samples_per_ui']; dt = model['dt']; ui = model['ui']
    stride = 1 if model['modulation'] == 'nrz' else 2
    # Shared bit guards keep elapsed time and the retained data identical.
    guard_bits = max(8192, 2 * int(np.ceil(model['impulse_ns'] * model['rate_gbps'])))
    guard_bits += guard_bits % 2
    all_bits = np.random.default_rng(seed).integers(0, 2, bits + 2 * guard_bits, dtype=np.uint8)
    symbols = 2 * all_bits.astype(float) - 1 if stride == 1 else gray_encode(all_bits)
    ffe_symbols = np.convolve(symbols, model['ffe'], mode='full')[:len(symbols)]
    boundary, rj, sj = jittered_boundaries(len(symbols), ui, seed + 100, stride)
    tx = source_waveform(ffe_symbols * SOURCE_VPP / 2, boundary, dt, RISE_S, stop_s=len(symbols) * ui)
    received = fftconvolve(tx, model['h_eq'], mode='full')
    indices = model['clock'] + np.arange(len(symbols)) * p
    decision_input = received[indices]
    decisions, feedback, corrected = dfe_decisions(decision_input, model['dfe'], model['main'], model['modulation'])
    start, count = guard_bits // stride, bits // stride
    metric, eye = measure_eye(received, symbols, model['clock'], p, dt, start, count,
                              feedback, model['main'], histogram=histogram)
    retained = slice(start, start + count)
    decoded = (decisions[retained] > 0).astype(np.uint8) if stride == 1 else gray_decode(decisions[retained])
    bit_errors = int(np.count_nonzero(decoded != all_bits[guard_bits:guard_bits + bits]))
    symbol_errors = int(np.count_nonzero(decisions[retained] != symbols[retained]))
    metric.update(bit_errors=bit_errors, symbol_errors=symbol_errors,
                  empirical_ber=bit_errors / bits, empirical_ser=symbol_errors / count,
                  measured_bits=int(bits))
    shot_center = model['clock'] + (start + 200) * p
    shot_index = np.arange(shot_center - p // 2, shot_center + 12 * p + p // 2)
    feedback_index = np.floor_divide(shot_index - model['clock'] + p // 2, p)
    source_shot = np.arange((start + 200) * p, (start + 213) * p)
    arrays = {
        'bits': all_bits, 'symbols': symbols, 'ffe_symbols': ffe_symbols,
        'boundary_times_s': boundary, 'rj_s': rj, 'sj_s': sj,
        'h_channel': model['h_channel'], 'h_equalized': model['h_eq'],
        'source_pulse_v': model['source_pulse'], 'sbr_channel_v': model['sbr_channel'],
        'sbr_equalized_v': model['sbr_eq'], 'sbr_ffe_v': model['sbr_ffe'],
        'ffe_taps': model['ffe'], 'dfe_taps_v': model['dfe'], 'thresholds_v': model['thresholds'],
        'ctle': model['ctle'], 'decision_input_v': decision_input,
        'decision_corrected_v': corrected, 'decisions': decisions, 'feedback_v': feedback,
        'clock_sample': np.array(model['clock']), 'dt_s': np.array(dt),
        'samples_per_ui': np.array(p), 'retained_start_symbol': np.array(start),
        'retained_symbol_count': np.array(count), 'main_cursor_v': np.array(model['main']),
        'waveform_time_relative_ps': (shot_index - shot_center) * dt * 1e12,
        'waveform_before_dfe_v': received[shot_index],
        'waveform_after_dfe_v': received[shot_index] - feedback[feedback_index],
        'source_time_relative_ps': (source_shot - source_shot[0]) * dt * 1e12,
        'source_waveform_v': tx[source_shot],
        **eye,
    }
    jitter = {'sj_amplitude_ps': SJ_AMPLITUDE_S * 1e12, 'sj_frequency_mhz': SJ_FREQUENCY_HZ / 1e6,
              'rj_sigma_ps': RJ_SIGMA_S * 1e12, 'measured_rj_mean_ps': float(np.mean(rj[retained]) * 1e12),
              'measured_rj_std_ps': float(np.std(rj[retained]) * 1e12),
              'sj_cycles_in_payload': bits / (model['rate_gbps'] * 1e9) * SJ_FREQUENCY_HZ,
              'minimum_boundary_spacing_ps': float(np.min(np.diff(boundary)) * 1e12),
              'maximum_source_magnitude_v': float(np.max(np.abs(tx)))}
    return metric, arrays, jitter


def _case_key(modulation, rate):
    return f'{modulation}_{rate:g}gbps'


def plot_results(output, f, sdd21, cases, all_arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.grid': True,
                         'grid.alpha': .22, 'figure.dpi': 140, 'savefig.dpi': 180})
    def save(fig, filename):
        fig.tight_layout(); fig.savefig(output / filename, metadata={'Software': 'HW14 reproducible analysis'}); plt.close(fig)
    def db(value):
        return 20 * np.log10(np.maximum(np.abs(value), 1e-14))
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(f / 1e9, -db(sdd21), lw=1.4)
    for rate in sorted({c['bit_rate_gbps'] for c in cases.values()}):
        for mod, factor, style in [('NRZ', .5, '--'), ('PAM4', .25, ':')]:
            freq = rate * factor
            ax.axvline(freq, color='.5', lw=.8, ls=style)
            ax.plot(freq, np.interp(freq * 1e9, f, -db(sdd21)), 'o', ms=4)
    ax.set(xlabel='Frequency (GHz)', ylabel='Differential insertion loss (dB)',
           title='Measured 50 cm channel: Sdd21', xlim=(0, 45), ylim=(0, 65))
    save(fig, 'channel.png')
    rows = int(np.ceil(len(cases) / 2))
    fig, axes = plt.subplots(rows, 2, figsize=(10, 3.5 * rows), squeeze=False)
    for ax, (key, case) in zip(axes.ravel(), cases.items()):
        a = all_arrays[key]; metric = case['metrics']
        probability = a['histogram'] / case['measured_symbols']
        ax.set_facecolor(plt.get_cmap('viridis')(0))
        ax.imshow(probability, origin='lower', aspect='auto', interpolation='nearest', cmap='viridis',
                  extent=[a['phase_ps'][0], a['phase_ps'][-1], a['voltage_edges_v'][0], a['voltage_edges_v'][-1]],
                  norm=LogNorm(vmin=1 / case['measured_symbols'], vmax=max(.1, float(probability.max()))))
        ax.set(title=f"{case['modulation'].upper()} at {case['bit_rate_gbps']:g} Gbps\nEH {metric['worst_eye_height_v']*1e3:.2f} mV; EW {metric['worst_eye_width_ps']:.2f} ps",
               xlabel='Time from fixed decision clock (ps)', ylabel='Differential voltage (V)')
    for unused in axes.ravel()[len(cases):]:
        unused.set_visible(False)
    save(fig, 'comparison.png')
    for key, case in cases.items():
        a = all_arrays[key]; metric = case['metrics']
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        ax.set_facecolor(plt.get_cmap('viridis')(0))
        probability = a['histogram'] / case['measured_symbols']
        im = ax.imshow(probability, origin='lower', aspect='auto', interpolation='nearest', cmap='viridis',
                       extent=[a['phase_ps'][0], a['phase_ps'][-1], a['voltage_edges_v'][0], a['voltage_edges_v'][-1]],
                       norm=LogNorm(vmin=1 / case['measured_symbols'], vmax=max(.1, float(probability.max()))))
        for boundary in range(len(a['levels']) - 1):
            ax.plot(a['phase_ps'], a['level_max_v'][boundary], color='white', lw=.8, ls='--')
            ax.plot(a['phase_ps'], a['level_min_v'][boundary + 1], color='white', lw=.8, ls='--')
        for threshold in case['equalizers']['dfe']['thresholds_v']:
            ax.axhline(threshold, color='white', lw=.7, ls=':')
        ax.axvline(0, color='white', lw=.7, ls=':')
        ax.set(xlabel='Time from fixed decision clock (ps)', ylabel='Differential voltage (V)',
               title=f"{case['modulation'].upper()} at {case['bit_rate_gbps']:g} Gbps: FFE + channel + CTLE + DFE\nMinimum EH {metric['worst_eye_height_v']*1e3:.2f} mV; minimum EW {metric['worst_eye_width_ps']:.2f} ps ({metric['worst_eye_width_ui']:.3f} UI)")
        fig.colorbar(im, ax=ax, label='Probability per voltage bin')
        save(fig, f'eye_{key}.png')
    fig, axes = plt.subplots(rows, 2, figsize=(10, 3.5 * rows), squeeze=False)
    for ax, (key, case) in zip(axes.ravel(), cases.items()):
        a = all_arrays[key]
        ax.plot(f / 1e9, db(a['ctle']), label='CTLE')
        ax.plot(f / 1e9, db(sdd21 / 2 * a['ctle']), label='Channel + CTLE')
        ax.axhline(CTLE_PEAK_DB, color='.4', ls=':', label='Common peak gain cap')
        ax.set(title=f"{case['modulation'].upper()} / {case['bit_rate_gbps']:g} Gbps", xlabel='Frequency (GHz)', ylabel='Gain (dB)', xlim=(0, 45), ylim=(-50, 22))
        ax.legend(fontsize=8)
    for unused in axes.ravel()[len(cases):]:
        unused.set_visible(False)
    save(fig, 'frequency_responses.png')
    fig, axes = plt.subplots(rows, 2, figsize=(10, 3.5 * rows), squeeze=False)
    for ax, (key, case) in zip(axes.ravel(), cases.items()):
        a = all_arrays[key]; dt = case['dt_ps']; clock = case['design']['clock_time_ps']
        for name, label in [('sbr_channel_v', 'Channel'), ('sbr_equalized_v', 'Channel + CTLE'), ('sbr_ffe_v', 'FFE + channel + CTLE')]:
            ax.plot(np.arange(len(a[name])) * dt - clock, a[name], label=label)
        sample_idx = int(a['clock_sample']) + np.arange(6) * int(a['samples_per_ui'])
        ax.plot(sample_idx * dt - clock, a['sbr_ffe_v'][sample_idx], 'o', ms=3)
        ax.set(title=f"{case['modulation'].upper()} / {case['bit_rate_gbps']:g} Gbps", xlabel='Time from main cursor (ps)', ylabel='Differential voltage (V)', xlim=(-3 * case['ui_ps'], 8 * case['ui_ps']))
        ax.legend(fontsize=8)
    for unused in axes.ravel()[len(cases):]:
        unused.set_visible(False)
    save(fig, 'pulse_responses.png')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--channel', type=Path, default=root / 'instructions' / 'HW14_diff_channel_50cm.s4p')
    parser.add_argument('--output', type=Path, default=root / 'results')
    parser.add_argument('--bits', type=int, default=65536)
    parser.add_argument('--rates', type=float, nargs='+', default=[20.0, 50.0])
    parser.add_argument('--samples-per-ui', type=int, default=64)
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED)
    parser.add_argument('--impulse-ns', type=float, default=40.0)
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--design-from', type=Path, help='Evaluate saved metrics.json EQ and physical clock without retuning')
    parser.add_argument('--validation', action='store_true', help='Independent seed, doubled time resolution and doubled impulse length, with frozen EQ settings')
    args = parser.parse_args(argv)
    if args.bits < 1024 or args.bits % 2 or args.samples_per_ui < 16 or args.samples_per_ui % 2:
        parser.error('Use an even bit count >=1024 and even samples/UI >=16')
    if any(rate <= 0 for rate in args.rates) or not 0 < args.impulse_ns < 100:
        parser.error('Positive rates and impulse duration below the 100 ns frequency-record period are required')
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    f, s, z0 = read_touchstone(args.channel)
    if not np.isclose(z0, 50):
        raise ValueError('The assignment requires all single-ended ports referenced to 50 ohms')
    sdd21 = differential_forward(s)
    saved = json.loads(args.design_from.read_text(encoding='utf-8')) if args.design_from else None
    if saved:
        p = saved['parameters']
        expected = {'source_vpp_v': SOURCE_VPP, 'rise_time_ps': RISE_S * 1e12,
                    'sj_amplitude_ps': SJ_AMPLITUDE_S * 1e12, 'sj_frequency_mhz': SJ_FREQUENCY_HZ / 1e6,
                    'rj_sigma_ps': RJ_SIGMA_S * 1e12}
        if any(not np.isclose(p.get(name, np.nan), value) for name, value in expected.items()):
            raise ValueError('Saved physical source and jitter parameters do not match this model')
    cases, arrays, models, validation = {}, {}, {}, {}
    for rate in args.rates:
        for modulation in ('nrz', 'pam4'):
            key = _case_key(modulation, rate)
            print(f'Designing {key}', flush=True)
            model = (load_design(f, sdd21, saved['cases'][key], args.samples_per_ui, args.impulse_ns)
                     if saved else design_equalizers(f, sdd21, rate, modulation, args.samples_per_ui, args.impulse_ns))
            print(f'Simulating {key}: predicted noiseless eye {model["design"]["pda_eye_height_v"]:.6g} V', flush=True)
            metric, archive, jitter = simulate_record(model, args.bits, args.seed)
            case = {'modulation': modulation, 'bit_rate_gbps': float(rate),
                    'symbol_rate_gbaud': rate / (1 if modulation == 'nrz' else 2),
                    'ui_ps': model['ui'] * 1e12, 'dt_ps': model['dt'] * 1e12,
                    'measured_bits': args.bits, 'measured_symbols': args.bits // (1 if modulation == 'nrz' else 2),
                    'equalizers': model['equalizers'], 'design': model['design'], 'metrics': metric, 'jitter': jitter,
                    'archive': f'{key}.npz', 'eye_figure': f'eye_{key}.png'}
            cases[key], arrays[key], models[key] = case, archive, model
            np.savez_compressed(args.output / case['archive'], **archive)
            print(json.dumps({'case': key, **metric}), flush=True)
    if args.validation:
        for key, model in models.items():
            print(f'Validating frozen design {key}', flush=True)
            second, _, _ = simulate_record(model, args.bits, args.seed + 1, histogram=False)
            refined = fixed_design(f, sdd21, model, args.samples_per_ui * 2, args.impulse_ns)
            fine, _, _ = simulate_record(refined, args.bits, args.seed, histogram=False)
            duration = min(args.impulse_ns * 2, 90.0)
            extended = fixed_design(f, sdd21, model, args.samples_per_ui, duration)
            long, _, _ = simulate_record(extended, args.bits, args.seed, histogram=False)
            validation[key] = {'independent_seed': {'seed': args.seed + 1, 'metrics': second},
                               'time_refinement': {'samples_per_ui': args.samples_per_ui * 2, 'metrics': fine},
                               'impulse_extension': {'impulse_ns': duration, 'metrics': long}}
    comparison = []
    for rate in args.rates:
        nrz = cases[_case_key('nrz', rate)]['metrics']; pam = cases[_case_key('pam4', rate)]['metrics']
        il_nrz = float(-20 * np.log10(abs(np.interp(rate * .5e9, f, sdd21))))
        il_pam = float(-20 * np.log10(abs(np.interp(rate * .25e9, f, sdd21))))
        comparison.append({'bit_rate_gbps': rate, 'nrz_nyquist_ghz': rate / 2,
                           'pam4_nyquist_ghz': rate / 4, 'nrz_nyquist_loss_db': il_nrz,
                           'pam4_nyquist_loss_db': il_pam, 'nyquist_loss_difference_db': il_nrz - il_pam,
                           'larger_minimum_eye_height': 'nrz' if nrz['worst_eye_height_v'] > pam['worst_eye_height_v'] else 'pam4',
                           'larger_minimum_eye_width_ps': 'nrz' if nrz['worst_eye_width_ps'] > pam['worst_eye_width_ps'] else 'pam4',
                           'larger_minimum_eye_width_ui': 'nrz' if nrz['worst_eye_width_ui'] > pam['worst_eye_width_ui'] else 'pam4'})
    compliance = {
        'same_source_vpp_and_rise': True, 'same_absolute_sj_and_rj_law': True,
        'same_ffe_and_dfe_tap_counts': True, 'same_ctle_pole_count_and_peak_cap': True,
        'all_ffe_l1_within_limit': all(c['equalizers']['ffe']['l1_norm'] <= 1 + 1e-10 for c in cases.values()),
        'all_final_eyes_open': all(c['metrics']['worst_eye_height_v'] > 0 and c['metrics']['worst_eye_width_ps'] > 0 for c in cases.values()),
        'all_payload_decisions_correct': all(c['metrics']['bit_errors'] == 0 for c in cases.values()),
    }
    result = {
        'schema_version': 1,
        'software': {'solver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                     'numpy_version': np.__version__, 'frozen_design_source': args.design_from.name if args.design_from else None},
        'parameters': {'source_vpp_v': SOURCE_VPP, 'source_plane': 'Open-circuit differential Thevenin source; matched transfer Sdd21/2',
                       'rise_time_ps': RISE_S * 1e12, 'rise_definition': '0%-100% linear transition, identical duration for every transition height',
                       'sj_amplitude_ps': SJ_AMPLITUDE_S * 1e12, 'sj_frequency_mhz': SJ_FREQUENCY_HZ / 1e6,
                       'rj_sigma_ps': RJ_SIGMA_S * 1e12, 'seed': args.seed, 'bits': args.bits,
                       'samples_per_ui': args.samples_per_ui, 'impulse_ns': args.impulse_ns,
                       'ctle_maximum_absolute_gain_db': CTLE_PEAK_DB, 'ffe_taps': FFE_COUNT, 'dfe_taps': DFE_COUNT},
        'channel': {'input_filename': args.channel.name, 'input_sha256': hashlib.sha256(args.channel.read_bytes()).hexdigest(),
                    'frequency_points': len(f), 'frequency_min_ghz': f[0] / 1e9, 'frequency_max_ghz': f[-1] / 1e9,
                    'frequency_step_mhz': (f[1] - f[0]) / 1e6, 'idft_period_ns': 1e9 / (f[1] - f[0]),
                    'sdd21_dc_real': float(sdd21[0].real), 'sdd21_dc_imag': float(sdd21[0].imag),
                    'single_ended_reference_ohm': z0, 'differential_reference_ohm': 2 * z0},
        'measurement': {'eye_height': 'Minimum voltage of adjacent upper transmitted level minus maximum voltage of adjacent lower level at one shared fixed clock; all retained symbols, no clipping or percentiles',
                        'eye_width': 'Contiguous interval containing the fixed clock where both conditional boundaries straddle their fixed noiseless midpoint threshold; crossings linearly interpolated',
                        'pam4_summary': 'Worst of the three openings; report each opening as well',
                        'dfe': 'Recovered normalized decisions only; correction held on clock_n-UI/2 <= t < clock_n+UI/2; right endpoint uses the next hold',
                        'band_limit': 'Measured response is used directly through 80 GHz and set to zero above it; real Hermitian IDFT and causal truncation; no DC renormalization',
                        'scope': 'Finite random-record experiment with Tx timing jitter; no added voltage noise, no CDR, no FEC, no low-BER extrapolation or hardware validation'},
        'cases': cases, 'comparison': comparison, 'validation': validation, 'compliance': compliance,
        'runtime_seconds': time.perf_counter() - started,
    }
    np.savez_compressed(args.output / 'channel.npz', frequency_hz=f, sdd21=sdd21)
    (args.output / 'metrics.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    if not args.no_plots:
        plot_results(args.output, f, sdd21, cases, arrays)
    print(json.dumps({'comparison': comparison, 'compliance': compliance, 'runtime_seconds': result['runtime_seconds']}, indent=2), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
