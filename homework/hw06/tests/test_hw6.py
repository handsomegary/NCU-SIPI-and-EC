"""Small independent numerical tests for the HW6 implementation.

Run with ``python -m unittest discover -s homework/hw06/tests -v``.
Only NumPy and the Python standard library are needed.
"""

import importlib.util
import itertools
from pathlib import Path
import tempfile
import unittest

import numpy as np


BASE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hw6_under_test", BASE / "src" / "hw6.py")
hw6 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hw6)


class PDATests(unittest.TestCase):
    def setUp(self):
        self.p = np.array([.05, -.03, .4, .2, 1., .7, -.5, .01, .02, -.4, .05, 0.])
        self.n = 2

    def test_all_symbol_combinations_match_exact_inner_and_outer_bounds(self):
        result = hw6.pda(self.p, self.n, tolerance=.1)
        for q in range(2 * self.n):
            main_index = result["kstart"] + q
            all_indices = range(main_index % self.n, len(self.p), self.n)
            neighbors = [i for i in all_indices if i != main_index]
            contributions = [sum(self.p[i] * b for i, b in zip(neighbors, bits))
                             for bits in itertools.product((0, 1), repeat=len(neighbors))]
            zero, one = np.array(contributions), np.array(contributions) + self.p[main_index]
            self.assertAlmostEqual(result["wc0_exact"][q], zero.max(), places=14)
            self.assertAlmostEqual(result["wc1_exact"][q], one.min(), places=14)
            self.assertAlmostEqual(result["outer_min"][q], min(zero.min(), one.min()), places=14)
            self.assertAlmostEqual(result["outer_max"][q], max(zero.max(), one.max()), places=14)

    def test_tolerance_omits_both_signs_and_accumulates_slack(self):
        result = hw6.pda(self.p, self.n, tolerance=.1)
        q = self.n  # Main tap 1.0; positive neighbors .05, .4, .02, .05, negative -.5.
        self.assertAlmostEqual(result["wc1"][q], .5)
        self.assertAlmostEqual(result["wc0"][q], .4)
        self.assertAlmostEqual(result["omitted_positive"][q], .12)
        self.assertGreater(result["omitted_positive"][q], .1)
        np.testing.assert_allclose(result["wc0_exact"], result["wc0"] + result["omitted_positive"], atol=1e-15)
        np.testing.assert_allclose(result["wc1_exact"], result["wc1"] - result["omitted_negative"], atol=1e-15)
        # An exact negative tolerance boundary is retained; tiny positive taps
        # must never be added to the minimum-one boundary as negative ISI.
        p = np.array([.1, 0, .2, 0, 1., 0, -.1, 0, .01, 0])
        result = hw6.pda(p, 2, tolerance=.1)
        self.assertAlmostEqual(result["wc1"][2], .9)
        self.assertAlmostEqual(result["wc0"][2], .3)

    def test_temporal_origin_is_derived_from_main_cursor(self):
        result = hw6.pda(self.p, self.n)
        self.assertEqual(result["kmain"], 4)
        self.assertEqual(result["kstart"], 2)
        np.testing.assert_array_equal(result["main"], self.p[2:6])
        self.assertEqual(np.argmax(result["main"]), self.n)


class PulseTests(unittest.TestCase):
    @staticmethod
    def synthetic_spectrum():
        # Eight-UI Fourier period, 64 samples/UI, measured below Nyquist.
        return np.arange(256) / 8

    def test_idft_matches_explicit_positive_frequency_sum(self):
        f = self.synthetic_spectrum()
        delay = .5
        s21 = .8 * np.exp(-2j * np.pi * f * delay)
        p, dt, meta = hw6.single_bit_response(f, s21, samples_per_ui=64,
                                             bit_period=1, rise_time=.25, vpp=2,
                                             kaiser_alpha=0)
        spectrum = np.sinc(f) * np.sinc(f * .25) * np.exp(-1j * np.pi * f * 1.25) * s21
        indices = np.array([0, 16, 32, 40, 64, 80, 96, 128, 511])
        expected = np.array([meta["df_hz"] * (spectrum[0].real + 2 * np.real(
            np.sum(spectrum[1:] * np.exp(2j * np.pi * f[1:] * i * dt)))) for i in indices])
        np.testing.assert_allclose(p[indices], expected, rtol=1e-12, atol=2e-14)

    def test_identity_channel_has_causal_ramps_and_matched_half_gain(self):
        f = self.synthetic_spectrum()
        p, dt, _ = hw6.single_bit_response(f, np.ones_like(f), samples_per_ui=64,
                                          bit_period=1, rise_time=.25, vpp=2,
                                          kaiser_alpha=0)
        t = np.arange(len(p)) * dt
        expected = np.minimum(np.clip(t / .25, 0, 1), np.clip((1.25 - t) / .25, 0, 1))
        # Finite bandwidth rounds the ramp corners, so compare with the
        # analytic continuous trapezoid using its bandwidth error allowance.
        np.testing.assert_allclose(p, expected, atol=.015)
        self.assertAlmostEqual(p[40], 1., delta=.002)
        doubled, _, _ = hw6.single_bit_response(f, 2 * np.ones_like(f), samples_per_ui=64,
                                               bit_period=1, rise_time=.25, vpp=2,
                                               kaiser_alpha=0)
        np.testing.assert_allclose(doubled, 2 * p, atol=1e-15)
        self.assertAlmostEqual(p.sum() * dt, 1., places=13)


class ConvolutionAndEyeTests(unittest.TestCase):
    def test_linear_overlap_add_matches_direct_convolution_including_partial_block(self):
        bits = np.array([1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1], dtype=np.uint8)
        kernel = np.array([0, .2, .8, 1, .6, -.1, .03, .01])
        n = 3
        impulses = np.zeros(len(bits) * n)
        impulses[::n] = bits
        expected = np.convolve(impulses, kernel)
        for block in (1, 4, 5, 11, 17):
            actual = hw6.sectional_convolution(bits, kernel, n, block_bits=block)
            np.testing.assert_allclose(actual, expected, atol=2e-15)
            self.assertEqual(len(actual), len(impulses) + len(kernel) - 1)
        # A final nonzero bit proves the last partial block was retained.
        self.assertGreater(expected[-n - 1], 0)

    def test_aligned_eye_conditional_extrema_histogram_and_duplicate(self):
        n = 2
        p = np.array([0., 0., .1, 1., .2, 0., -.1, 0.])
        bits = np.random.default_rng(7).integers(0, 2, 80, dtype=np.uint8)
        impulses = np.zeros(len(bits) * n)
        impulses[::n] = bits
        y = np.convolve(impulses, p)
        boundaries = hw6.pda(p, n, tolerance=.15)
        edges = np.linspace(-1, 2, 25)
        result = hw6.analyze_eye(y, bits, boundaries, n, edges, sbr_length=len(p))
        guard = int(np.ceil(len(p) / n)) + 2
        labels = bits[guard:-guard].astype(bool)
        windows = np.array([[y[boundaries["kstart"] + j * n + q] for q in range(2 * n)]
                            for j in range(guard, len(bits) - guard)])
        np.testing.assert_allclose(result["observed_min_one"], windows[labels].min(axis=0))
        np.testing.assert_allclose(result["observed_max_zero"], windows[~labels].max(axis=0))
        self.assertTrue(np.all(result["observed_min_one"] >= boundaries["wc1_exact"] - 1e-14))
        self.assertTrue(np.all(result["observed_max_zero"] <= boundaries["wc0_exact"] + 1e-14))
        expected_hist = np.column_stack([np.histogram(windows[:, q], bins=edges)[0] for q in range(n)])
        np.testing.assert_array_equal(result["histogram"], expected_hist)
        np.testing.assert_array_equal(result["histogram"].sum(axis=0), len(labels))
        np.testing.assert_allclose(result["probability"].sum(axis=0), 1, atol=1e-14)
        duplicated = np.hstack((result["probability"], result["probability"]))
        np.testing.assert_array_equal(duplicated[:, :n], duplicated[:, n:])
        self.assertEqual(result["summary"]["maximum_exact_bound_violation_v"], 0.)


class ReaderTests(unittest.TestCase):
    def test_touchstone_uses_s21_order_for_nonreciprocal_network(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as folder:
            path = Path(folder) / "example.s2p"
            path.write_text("# MHz S RI R 75\n0 0 0 .7 -.2 .3 .1 0 0\n"
                            "10 0 0 .6 -.3\n .2 .1 0 0 ! continuation\n", encoding="utf-8")
            f, s21, z0 = hw6.read_touchstone(path)
        np.testing.assert_array_equal(f, [0, 10e6])
        np.testing.assert_allclose(s21, [.7 - .2j, .6 - .3j])
        self.assertEqual(z0, 75)

    def test_bit_reader_accepts_whitespace_and_rejects_other_characters(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as folder:
            path = Path(folder) / "bits.txt"
            path.write_text("01 10\n1\t0", encoding="utf-8")
            np.testing.assert_array_equal(hw6.read_bits(path), [0, 1, 1, 0, 1, 0])
            path.write_text("010x", encoding="utf-8")
            with self.assertRaises(ValueError):
                hw6.read_bits(path)


if __name__ == "__main__":
    unittest.main()
