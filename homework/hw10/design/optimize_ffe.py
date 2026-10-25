"""Reproduce the finite-training-pattern FFE search used for the HW10 design.

The linear program assumes correct past decisions during design only. The final
solver validates a recursive DFE on separate patterns and a finer time grid.
This search is not an all-pattern guarantee. Published coefficients are fixed
in src/hw10.py so ordinary reproduction does not need to repeat optimization.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.fft import irfft
from scipy.optimize import linprog
from scipy.signal import fftconvolve


def read_channel(path):
    numbers = []
    for raw in path.read_text().splitlines():
        line = raw.split('!')[0].strip()
        if line and not line.startswith('#'):
            numbers.extend(map(float, line.split()))
    data = np.asarray(numbers).reshape(-1, 33)
    s = (data[:, 1::2] * np.exp(1j * np.deg2rad(data[:, 2::2]))).reshape(-1, 4, 4)
    matched = (s[:, 1, 0] + s[:, 3, 2] - s[:, 1, 2] - s[:, 3, 0]) / 4
    return data[:, 0], matched


def search(channel, seed=70, bits=16384):
    frequency, matched = read_channel(channel)
    dt, samples_per_ui, nfft = .5e-12, 50, 200000
    ui, rise = 25e-12, 5e-12
    window = np.i0(np.pi * np.sqrt(np.maximum(0, 1 - (frequency / frequency[-1]) ** 2))) / np.i0(np.pi)
    pulse = 2 * ui * np.sinc(frequency * ui) * np.sinc(frequency * rise) * np.exp(-1j * np.pi * frequency * (ui + rise))
    symbols = 2 * np.random.default_rng(seed).integers(0, 2, bits) - 1
    selected = np.arange(1024, bits - 1024)
    candidates = []
    for dc_db in (-14, -12, -10, -8, -6):
        pole = 2 * np.pi * 40e9
        adc, aac = 10 ** (dc_db / 20), 10 ** (14 / 20)
        s = 2j * np.pi * frequency
        ctle = aac * pole ** 2 * (s + adc * pole / aac) / (s + pole) ** 3
        # Search on the sampled 6 dB boundary. The released design then lowers
        # the entire CTLE by 0.01 dB and verifies its continuous-frequency peak.
        scale = 10 ** (6 / 20) / max(abs(ctle))
        ctle *= scale
        spectrum = np.zeros(nfft // 2 + 1, complex)
        spectrum[:len(frequency)] = matched * ctle * window * pulse
        base = (irfft(spectrum, n=nfft) / dt)[:40000]
        peak = int(np.argmax(base))
        for offset in (-4, 0, 4):
            clock = peak + 2 * samples_per_ui + offset
            post = base[clock + np.arange(1, 6)[:, None] * samples_per_ui - np.arange(5)[None, :] * samples_per_ui]
            correction = np.stack([
                sum(symbols[selected - i] * post[i - 1, j] for i in range(1, 6))
                for j in range(5)
            ], axis=1)
            rows, bounds = [], []
            for delta, margin in ((-16, .004), (-8, .015), (0, 0), (8, .015), (16, .004)):
                quotient, phase = divmod(clock + delta, samples_per_ui)
                waveform = fftconvolve(symbols, base[phase::samples_per_ui])
                coefficient = np.stack([waveform[selected + quotient - j] for j in range(5)], axis=1)
                coefficient = (coefficient - correction) * symbols[selected, None]
                block = np.zeros((len(selected), 11))
                block[:, :5] = -coefficient
                if delta == 0:
                    block[:, -1] = 1
                rows.append(block)
                bounds.extend(np.full(len(selected), -margin))
            # Variables: five taps, their five absolute-value bounds, and
            # the worst signed voltage at the decision point.
            absolute = np.zeros((11, 11))
            for j in range(5):
                absolute[2 * j, j] = 1
                absolute[2 * j, 5 + j] = -1
                absolute[2 * j + 1, j] = -1
                absolute[2 * j + 1, 5 + j] = -1
            absolute[-1, 5:10] = 1
            rows.append(absolute)
            bounds.extend([0] * 10 + [1])
            objective = np.zeros(11)
            objective[-1] = -1
            solution = linprog(objective, A_ub=np.vstack(rows), b_ub=np.asarray(bounds),
                               bounds=[(-1, 1)] * 5 + [(0, 1)] * 5 + [(0, None)], method='highs')
            if solution.success:
                candidates.append({
                    'ctle_dc_before_headroom_db': float(dc_db + 20 * np.log10(scale)),
                    'ctle_numerator_before_headroom_db': float(14 + 20 * np.log10(scale)),
                    'headroom_db_to_subtract': .01,
                    'poles_ghz': [40, 40, 40],
                    'clock_phase_ps': (clock % samples_per_ui) * .5,
                    'five_tap_vector': solution.x[:5].tolist(),
                    'minimum_signed_training_voltage_v': float(solution.x[-1]),
                })
    return {'training_seed': seed, 'generated_bits': bits, 'retained_training_bits': len(selected),
            'design_assumes_correct_feedback': True, 'candidates': candidates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    home = Path(__file__).resolve().parents[1]
    parser.add_argument('--channel', type=Path, default=home / 'instructions' / 'HW10_complicated_diff_channel.s4p')
    parser.add_argument('--output', type=Path, default=home / 'results' / 'design_search.json')
    parser.add_argument('--seed', type=int, default=70)
    parser.add_argument('--bits', type=int, default=16384)
    args = parser.parse_args()
    if args.bits < 4096:
        parser.error('Use at least 4096 training bits')
    result = search(args.channel, args.seed, args.bits)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(f"Saved {len(result['candidates'])} feasible training candidates to {args.output.name}")


if __name__ == '__main__':
    main()
