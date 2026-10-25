"""Independently audit the saved HW11 numerical artifacts.

Run from any directory. The default audit reads ../results and writes
verification.json beside this script. A reduced-bit CI run can provide
--results and --output explicitly. No external circuit simulator is claimed.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import unittest

import numpy as np
from numpy.testing import assert_allclose

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parent
sys.path.insert(0, str(TESTS))
from reference import direct_mimo_convolution, scalar_nrz, wave_terminated_transfer


def quantile_opening(counts, edges, threshold):
    """Reconstruct observations and use NumPy's empirical inverse CDF.

    This deliberately does not call the production histogram metric helper.
    The finite-record voltage-cluster metric is not conditioned on Tx bits.
    """
    centers = (edges[:-1] + edges[1:]) / 2
    gaps = np.full(counts.shape[1], -np.inf)
    boundaries = {}
    for phase in range(counts.shape[1]):
        observations = np.repeat(centers, counts[:, phase].astype(np.int64))
        lo = observations[observations < threshold]
        hi = observations[observations >= threshold]
        if min(len(lo), len(hi)) < .2 * len(observations):
            continue
        lower = float(np.quantile(lo, .999, method='inverted_cdf'))
        upper = float(np.quantile(hi, .001, method='inverted_cdf'))
        gaps[phase] = upper-lower
        boundaries[phase] = (lower, upper)
    return gaps, boundaries


def verify(results):
    metrics = json.loads((results / 'metrics.json').read_text(encoding='utf-8'))
    config = metrics['configuration']
    analysis_sha = hashlib.sha256((ROOT / 'analyze.py').read_bytes()).hexdigest()
    assert analysis_sha == metrics['reproducibility']['analysis_source_sha256']
    with np.load(results / 'channel.npz', allow_pickle=False) as archive:
        channel = {key: archive[key] for key in archive.files}
    with np.load(results / 'simulation.npz', allow_pickle=False) as archive:
        simulation = {key: archive[key] for key in archive.files}
    frequency, scattering, transfer = (channel[k] for k in ['frequency_hz', 's', 'h'])
    npui, nbits, dt = config['samples_per_ui'], config['bits'], config['dt_s']
    assert frequency.shape == (4001,)
    assert scattering.shape == (4001, 16, 16)
    assert transfer.shape == (4001, 8, 8)
    assert_allclose(np.diff(frequency), 10e6, atol=0, rtol=0)
    assert_allclose(frequency[[0, -1]], [0, 40e9], atol=0, rtol=0)
    source_sha = hashlib.sha256((ROOT / 'instructions' / 'HW11_crosstalk_8lines.s16p').read_bytes()).hexdigest()
    assert source_sha == metrics['input']['sha256']
    wave = wave_terminated_transfer(frequency, scattering, metrics['input']['z0_ohm'],
                                    config['rs_ohm'], config['cs_f'],
                                    config['rl_ohm'], config['cl_f'])
    transfer_error = float(np.max(np.abs(transfer-wave)))
    assert_allclose(transfer, wave, atol=3e-13, rtol=2e-10)

    # A direct trigonometric sum checks selected signed impulse lags across
    # all 64 paths; this is independent of any inverse-FFT implementation.
    taps = channel['impulse_coefficients'].astype(float)
    first_lag = int(channel['first_lag'])
    transform_length = int(round(1 / (dt*(frequency[1]-frequency[0]))))
    alpha = config['kaiser_a']
    window = np.i0(np.pi*alpha*np.sqrt(np.maximum(0, 1-(frequency/frequency[-1])**2))) / np.i0(np.pi*alpha)
    peak_index = int(np.argmax(np.max(np.abs(taps), axis=(1, 2))))
    indices = np.unique(np.clip(np.r_[np.linspace(0, len(taps)-1, 9).astype(int),
                                      -first_lag, peak_index+np.arange(-4, 5)], 0, len(taps)-1))
    trig = []
    for index in indices:
        lag = int(index) + first_lag
        terms = transfer * window[:, None, None]
        weights = np.exp(2j*np.pi*frequency*lag*dt)
        # The supplied measured band is strictly below the sampling Nyquist.
        assert frequency[-1] < .5/dt
        coefficient = (terms[0].real + 2*np.real(np.einsum('fij,f->ij', terms[1:], weights[1:]))) / transform_length
        trig.append(coefficient)
    trig_error = float(np.max(np.abs(taps[indices]-np.asarray(trig))))
    assert_allclose(taps[indices], trig, atol=2e-9, rtol=2e-6)

    bits = simulation['bits']
    assert bits.shape == (nbits + 2*config['guard_ui'] + 2, 8)
    assert np.isin(bits, [0, 1]).all()
    retained_bits = bits[config['guard_ui']:config['guard_ui']+nbits]
    source_correlation = np.corrcoef(retained_bits.astype(float), rowvar=False)
    max_correlation = float(np.max(np.abs(source_correlation-np.eye(8))))
    assert max_correlation < .09
    edges = simulation['voltage_edges_v']
    assert len(edges) == 3001
    assert np.all(np.diff(edges) > 0)
    conservation = {}
    for name, receivers in [('all_on', list(range(8))), ('only_tx0', [0]),
                             ('only_tx4', [1, 2, 3, 4])]:
        counts = simulation[name + '_counts']
        assert counts.shape == (len(receivers), len(edges)-1, npui)
        assert np.issubdtype(counts.dtype, np.integer) and np.all(counts >= 0)
        assert np.all(counts.sum(axis=1) == nbits)
        assert np.all(counts.sum(axis=(1, 2)) == nbits*npui)
        conservation[name] = {'receivers': receivers,
                              'samples_per_receiver': int(nbits*npui),
                              'samples_per_phase': int(nbits)}
        for local, receiver in enumerate(receivers):
            metric = metrics['cases'][name][f'rx{receiver}']
            assert metric['outside_voltage_range'] == 0
            assert metric['samples'] == nbits*npui
            assert metric['minimum_samples_per_phase'] == nbits
            # FEXT-only receiver channels do not have two logic clusters.
            if name == 'only_tx4' and receiver != 4:
                continue
            threshold = transfer[0, receiver, receiver].real / 2
            assert_allclose(metric['threshold_v'], threshold, atol=1e-15)
            gaps, bounds = quantile_opening(counts[local], edges, threshold)
            assert_allclose(metric['opening_v'], np.max(gaps), atol=1e-12)
            phase = int(round(metric['best_phase_ui']*npui))
            assert_allclose(gaps[phase], metric['opening_v'], atol=1e-12)
            assert_allclose(bounds[phase], [metric['lower_cluster_99_9_percentile_v'],
                                            metric['upper_cluster_0_1_percentile_v']], atol=1e-12)

    times = simulation['snapshot_times_s']
    snapshot = simulation['tx4_only_snapshot_v']
    assert snapshot.shape == (round(4e-9/dt), 4)
    assert_allclose(np.diff(times), dt, atol=1e-19, rtol=1e-10)
    assert_allclose(times[0], 145e-9, atol=1e-19)
    assert_allclose(simulation['only_tx4_receivers'], [1, 2, 3, 4], atol=0, rtol=0)
    physical_start = config['guard_ui']*npui + round(times[0]/dt)
    last_lag = first_lag+len(taps)-1
    source_start = physical_start-last_lag
    source_stop = physical_start+len(times)-first_lag
    rise_samples = round(config['Tr_s']/dt)
    delays = np.rint(np.asarray(config['delays_s'])/dt).astype(int)
    source = scalar_nrz(bits[:, [4]], npui, rise_samples, delays[[4]], source_stop)
    source = source[source_start:source_stop, 0]
    direct_snapshot = np.column_stack([np.convolve(source, taps[:, rx, 4], mode='valid')
                                       for rx in [1, 2, 3, 4]])
    snapshot_error = float(np.max(np.abs(snapshot-direct_snapshot)))
    assert_allclose(snapshot, direct_snapshot, atol=2e-6, rtol=2e-4)

    # Also record an absolute double-precision OLA/direct oracle error.
    import test_hw11
    rng = np.random.default_rng(123)
    x = rng.normal(size=(37, 3))
    fixture = rng.normal(size=(11, 2, 3))
    kernel = test_hw11.hw11.ImpulseKernel(fixture, -4, .2, {})
    actual = test_hw11.hw11.overlap_add_mimo(x, kernel, block_size=8)
    convolution_error = float(np.max(np.abs(actual-direct_mimo_convolution(x, fixture))))
    return {
        'status': 'passed',
        'scope': 'Independent mathematical and saved-artifact checks; no ADS execution.',
        'retained_bits_per_transmitter': int(nbits),
        'loaded_transfer_checked_frequencies': int(len(frequency)),
        'loaded_transfer_checked_paths': 64,
        'max_mna_vs_wave_boundary_abs_error': transfer_error,
        'max_saved_impulse_vs_direct_trigonometric_sum_error': trig_error,
        'max_float64_overlap_add_vs_direct_convolution_error': convolution_error,
        'max_tx4_snapshot_vs_direct_convolution_error_v': snapshot_error,
        'maximum_absolute_offdiagonal_source_bit_correlation': max_correlation,
        'histogram_conservation': conservation,
        'percentile_openings_recomputed_from_expanded_observations': 10,
        'source_sha256': source_sha,
        'analysis_source_sha256': analysis_sha,
        'metrics_sha256': hashlib.sha256((results / 'metrics.json').read_bytes()).hexdigest(),
        'audited_artifact_sha256': {
            name: hashlib.sha256((results / name).read_bytes()).hexdigest()
            for name in ['metrics.json', 'channel.npz', 'simulation.npz']
        },
        'limits': [
            'Empirical voltage-cluster percentile openings are neither BER estimates nor guaranteed worst-case openings.',
            'Finite measured bandwidth, Kaiser weighting, sample interval and signed impulse support remain numerical model approximations.',
            'Saved snapshot uses float32 convolution; independent reconstruction uses float64 and the saved float32 impulse coefficients.',
            'The all-eight full-record histogram counts are conserved; the entire full-record waveform is not stored or independently reconstructed.'
        ]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT/'results')
    parser.add_argument('--output', type=Path, default=TESTS/'verification.json')
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.discover(str(TESTS), pattern='test_hw11.py')
    log = io.StringIO()
    outcome = unittest.TextTestRunner(stream=log, verbosity=1).run(suite)
    if not outcome.wasSuccessful():
        print(log.getvalue())
        raise SystemExit(1)
    summary = verify(args.results)
    summary['unit_tests'] = {'run': outcome.testsRun,
                             'failures': len(outcome.failures), 'errors': len(outcome.errors)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
