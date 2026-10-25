"""Independent regression tests for the eight-channel crosstalk solver.

Run: python -m unittest discover -s homework/hw11/tests -v
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference import (direct_mimo_convolution, ideal_through_rc_transfer,
                       ideal_through_scattering, scalar_eye_counts, scalar_nrz,
                       wave_terminated_transfer)

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('hw11_under_test', ROOT / 'analyze.py')
hw11 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = hw11
SPEC.loader.exec_module(hw11)


class TouchstoneTests(unittest.TestCase):
    def load(self, contents, extension):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tests') as folder:
            path = Path(folder) / ('fixture.' + extension)
            path.write_text(contents, encoding='ascii')
            return hw11.load_touchstone(path)

    def test_nonreciprocal_four_port_row_order(self):
        # Four-port ordering differs from the special two-port convention.
        s = np.arange(1, 17).reshape(4, 4) / 100
        pairs = ' '.join(f'{v:.5f} {-v/2:.5f}' for v in s.ravel())
        freq, actual, z0 = self.load(
            '# MHz S RI R 75\n0 ' + pairs + '\n10 ' + pairs + '\n', 's4p')
        assert_allclose(freq, [0, 10e6])
        assert_allclose(actual[0], s * (1 - 0.5j))
        self.assertEqual(z0, 75)
        self.assertNotEqual(actual[0, 1, 2], actual[0, 2, 1])

    def test_two_port_special_order_ma_units_and_degrees(self):
        freq, s, z0 = self.load(
            '! asymmetric example\n# GHz S MA R 50\n'
            '0 0.1 0 0.8 90 0.2 -90 0.3 180\n'
            '1 0.1 0 0.8 90 0.2 -90 0.3 180\n', 's2p')
        assert_allclose(freq, [0, 1e9])
        assert_allclose(s[0], [[0.1, -0.2j], [0.8j, -0.3]], atol=1e-15)
        self.assertEqual(z0, 50)

    def test_db_inline_comments_and_multiline_records(self):
        freq, s, _ = self.load(
            '# kHz S DB R 50\n0 -20 0 -6.020599913279624 180\n'
            ' -40 90 0 0 ! trailing comment\n'
            '100 -20 0 -6.020599913279624 180 -40 90 0 0\n', 's2p')
        assert_allclose(freq, [0, 100e3])
        assert_allclose(s[0], [[0.1, 0.01j], [-0.5, 1]], atol=1e-14)

    def test_incomplete_record_rejected(self):
        with self.assertRaises(ValueError):
            self.load('# Hz S RI R 50\n0 1 0 2 0\n', 's4p')


class LoadedCircuitTests(unittest.TestCase):
    def test_ideal_through_reduces_to_scalar_rc_divider(self):
        frequency = np.array([0, 1e8, 1e9, 10e9, 40e9])
        s = np.repeat(ideal_through_scattering(3)[None], len(frequency), axis=0)
        actual = hw11.solve_transfer(frequency, s, 50)
        expected = ideal_through_rc_transfer(frequency)
        for i in range(3):
            assert_allclose(actual[:, i, i], expected, atol=2e-14, rtol=2e-13)
            for j in range(3):
                if i != j:
                    assert_allclose(actual[:, i, j], 0, atol=2e-14)
        self.assertAlmostEqual(actual[0, 0, 0].real, 65 / 95, places=13)

    def test_matched_terminals_give_half_s_wave_gain(self):
        rng = np.random.default_rng(42)
        frequency = np.array([0, 1e9, 20e9])
        s = (rng.normal(size=(3, 6, 6)) + 1j*rng.normal(size=(3, 6, 6))) / 20
        actual = hw11.solve_transfer(frequency, s, 50, rs=50, cs=0, rl=50, cl=0)
        assert_allclose(actual, s[:, 1::2, 0::2] / 2, atol=1e-14, rtol=1e-13)

    def test_nonreciprocal_network_against_wave_boundary_solution(self):
        rng = np.random.default_rng(118)
        frequency = np.array([0, 250e6, 5e9, 40e9])
        s = (rng.normal(size=(4, 16, 16)) + 1j*rng.normal(size=(4, 16, 16))) / 40
        actual = hw11.solve_transfer(frequency, s, 50)
        expected = wave_terminated_transfer(frequency, s, 50)
        assert_allclose(actual, expected, atol=3e-14, rtol=2e-12)

    def test_supplied_channel_against_independent_wave_solution(self):
        frequency, s, z0 = hw11.load_touchstone(
            ROOT / 'instructions' / 'HW11_crosstalk_8lines.s16p')
        self.assertEqual(s.shape, (4001, 16, 16))
        self.assertEqual(frequency[0], 0)
        self.assertEqual(frequency[-1], 40e9)
        assert_allclose(np.diff(frequency), 10e6, atol=0, rtol=0)
        chosen = np.unique(np.r_[np.arange(0, len(frequency), 199), len(frequency)-1])
        actual = hw11.solve_transfer(frequency[chosen], s[chosen], z0)
        expected = wave_terminated_transfer(frequency[chosen], s[chosen], z0)
        assert_allclose(actual, expected, atol=3e-13, rtol=2e-10)


class FrequencyTransformTests(unittest.TestCase):
    def test_full_band_pure_delay_has_unit_discrete_area(self):
        frequency = np.arange(9) / 16
        transfer = np.exp(-2j*np.pi*frequency*3)[:, None, None]
        kernel = hw11.frequency_to_impulse(frequency, transfer, dt=1,
                                           kaiser_a=0, pre_time=2, post_time=5)
        expected = np.zeros_like(kernel.samples)
        expected[3 - kernel.first_lag, 0, 0] = 1
        assert_allclose(kernel.samples, expected, atol=2e-14)


class ConvolutionTests(unittest.TestCase):
    def test_mimo_overlap_add_matches_direct_convolution_at_every_sample(self):
        rng = np.random.default_rng(123)
        x = rng.normal(size=(37, 3))
        taps = rng.normal(size=(11, 2, 3))
        kernel = hw11.ImpulseKernel(taps, first_lag=-4, dt=0.2, diagnostics={})
        expected = direct_mimo_convolution(x, taps)
        for block_size in [1, 8, 16, 100]:
            with self.subTest(block_size=block_size):
                actual = hw11.overlap_add_mimo(x, kernel, block_size=block_size)
                self.assertEqual(actual.shape, (47, 2))
                assert_allclose(actual, expected, atol=8e-14, rtol=8e-13)

    def test_individual_tx_outputs_sum_to_simultaneous_result(self):
        rng = np.random.default_rng(124)
        x = rng.normal(size=(43, 4))
        taps = rng.normal(size=(9, 4, 4))
        kernel = hw11.ImpulseKernel(taps, first_lag=-3, dt=1, diagnostics={})
        total = hw11.overlap_add_mimo(x, kernel, block_size=13)
        separate = np.zeros_like(total)
        for j in range(4):
            single = np.zeros_like(x)
            single[:, j] = x[:, j]
            separate += hw11.overlap_add_mimo(single, kernel, block_size=13)
        assert_allclose(total, separate, atol=8e-14, rtol=1e-12)


class SourceTests(unittest.TestCase):
    def test_rng_is_reproducible_and_per_channel(self):
        bits = hw11.independent_bits(8192, 8, seed=832)
        self.assertEqual(bits.shape, (8192, 8))
        self.assertTrue(np.all((bits == 0) | (bits == 1)))
        np.testing.assert_array_equal(bits, hw11.independent_bits(8192, 8, seed=832))
        np.testing.assert_array_equal(bits[:, :4], hw11.independent_bits(8192, 4, seed=832))
        self.assertFalse(np.array_equal(bits, hw11.independent_bits(8192, 8, seed=833)))
        correlation = np.corrcoef(bits.astype(float), rowvar=False)
        self.assertLess(np.max(np.abs(correlation - np.eye(8))), 0.06)

    def test_finite_rise_delay_and_no_circular_wrap(self):
        bits = np.array([1, 0, 1], dtype=np.uint8)
        actual = hw11.nrz_waveform(bits, npui=4, rise_samples=2,
                                   delay_samples=2, start=-2, stop=17)
        # Absolute sample positions -2..16. The delay is zero filling,
        # not a circular roll. The finite bit record is zero outside
        # its support; production simulations discard guarded boundaries.
        expected = np.array([0, 0, 0, 0, 0, .5, 1, 1, 1, .5,
                             0, 0, 0, .5, 1, 1, 0, 0, 0])
        assert_allclose(actual, expected, atol=0, rtol=0)

    def test_nrz_chunk_boundaries_preserve_samples_and_all_tail_bits(self):
        bits = np.array([[1, 0], [0, 1], [1, 1], [1, 0], [1, 1]], np.uint8)
        kwargs = dict(npui=10, rise_samples=2, delay_samples=np.array([0, 5]))
        full = hw11.nrz_waveform(bits, **kwargs, start=0, stop=55)
        chunks = [hw11.nrz_waveform(bits, **kwargs, start=a, stop=b)
                  for a, b in [(0, 13), (13, 27), (27, 47), (47, 55)]]
        assert_allclose(np.concatenate(chunks), full, atol=0, rtol=0)
        self.assertEqual(full[49, 0], 1)
        self.assertEqual(full[54, 1], 1)
        assert_allclose(full[:5, 1], 0, atol=0, rtol=0)


class EyeHistogramTests(unittest.TestCase):
    def test_every_sample_bin_and_per_phase_normalization(self):
        y = np.array([[-.1, .9], [.2, .7], [.8, .1], [.9, -.1],
                      [.4, .4], [.6, .6], [.1, .8], [.7, .2], [.3, .5]])
        histogram = hw11.EyeHistogram(3, 2, voltage_min=-.2,
                                     voltage_max=1, voltage_bins=12)
        histogram.update(y[:4], start_sample=0)
        histogram.update(y[4:], start_sample=4)
        edges = np.linspace(-.2, 1, 13)
        for j in range(2):
            expected = scalar_eye_counts(y[:, j], 3, edges)
            np.testing.assert_array_equal(histogram.counts[j], expected)
        probability = histogram.probabilities()
        assert_allclose(probability.sum(axis=1), 1, atol=1e-15)
        self.assertEqual(int(histogram.counts.sum()), y.size)


class StreamingCaseTests(unittest.TestCase):
    def test_guard_alignment_and_last_retained_sample_with_negative_lags(self):
        rng = np.random.default_rng(730)
        bits = rng.integers(0, 2, size=(20, 2), dtype=np.uint8)
        taps = np.zeros((6, 2, 2))
        taps[0] = [[.017, -.007], [-.009, .013]]
        taps[2] = [[.57, .021], [.023, .61]]
        taps[5] = [[.083, -.011], [-.013, .071]]
        kernel = hw11.ImpulseKernel(taps, first_lag=-2, dt=1e-9, diagnostics={})
        npui, nbits, guard_ui = 4, 9, 3
        delays = np.array([0, 2])
        full_x = scalar_nrz(bits, npui, 2, delays, len(bits)*npui)
        for transmitters in [[0, 1], [1]]:
            with self.subTest(transmitters=transmitters):
                # Reorder the receiver list to test explicit receiver selection.
                receivers = [1, 0]
                expected_full = direct_mimo_convolution(
                    full_x[:, transmitters], taps[:, receivers][:, :, transmitters])
                raw_start = guard_ui*npui - kernel.first_lag
                expected = expected_full[raw_start:raw_start + nbits*npui]
                eye, times, snapshot, stats = hw11.simulate_case(
                    bits, kernel, npui, 2, delays, nbits, guard_ui,
                    receivers, transmitters, block_size=7,
                    snapshot_start=0, snapshot_duration=nbits*npui*kernel.dt)
                assert_allclose(snapshot, expected, atol=2e-7, rtol=2e-6)
                assert_allclose(snapshot[-1], expected[-1], atol=2e-7, rtol=2e-6)
                assert_allclose(times, np.arange(nbits*npui)*kernel.dt, atol=1e-20)
                assert_allclose(stats['mean_v'], expected.mean(axis=0), atol=2e-7)
                np.testing.assert_array_equal(eye.counts.sum(axis=1), nbits)
                self.assertEqual(stats['retained_samples_per_receiver'], nbits*npui)


if __name__ == '__main__':
    unittest.main()
