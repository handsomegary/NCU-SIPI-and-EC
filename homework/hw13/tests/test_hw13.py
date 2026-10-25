"""Behavioral tests for the numerical analysis, using independent erfc oracles."""
from __future__ import annotations

import csv
import importlib.util
import io
import json
import math
from pathlib import Path
import platform
import sys
import tempfile
import unittest

import numpy as np

BASE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hw13_under_test", BASE / "src/hw13.py")
hw = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = hw
SPEC.loader.exec_module(hw)
SQRT2 = math.sqrt(2)


def mixture_tail(x, p, left=True):
    m1, s1, m2, s2, w = p
    sign = -1 if left else 1
    return (w * .5 * math.erfc(sign * (x - m1) / (s1 * SQRT2))
            + (1 - w) * .5 * math.erfc(sign * (x - m2) / (s2 * SQRT2)))


def reference_ber(t, p):
    return .5 * (mixture_tail(t - 12.5, p, False) + mixture_tail(t - 37.5, p))


def bisect_root(fn, a, b):
    fa = fn(a)
    for _ in range(100):
        mid = (a + b) / 2
        if fn(mid) * fa > 0:
            a, fa = mid, fn(mid)
        else:
            b = mid
    return (a + b) / 2


class GaussianModelTests(unittest.TestCase):
    def test_extreme_left_and_right_tails_keep_tiny_probability(self):
        p = (0, 1, 0, 1, .37)
        for z in (8., 20., 37., 38.):
            expected = .5 * math.erfc(z / SQRT2)
            for x, left in ((-z, True), (z, False)):
                actual = float(np.exp(hw.log_tail(x, p, left)))
                self.assertGreater(actual, 0)
                self.assertAlmostEqual(actual / expected, 1, delta=2e-11)
        self.assertTrue(math.isfinite(float(hw.log_tail(-100, p))))

    def test_two_gaussian_tails_match_independent_erfc(self):
        p = (-2.1, 1.3, 1.4, 1.7, .32)
        for x in (-20, -8, -2, 0, 3, 10, 20):
            for left in (True, False):
                expected = mixture_tail(x, p, left)
                self.assertAlmostEqual(float(np.exp(hw.log_tail(x, p, left))) / expected, 1, delta=1e-12)

    def test_zero_weight_components_are_valid_single_gaussians(self):
        for p in ((0, 2, 5, 3, 1), (-7, 3, 0, 2, 0)):
            np.testing.assert_allclose(hw.model_ber([15, 25, 35], p),
                                       [reference_ber(t, p) for t in (15, 25, 35)], rtol=2e-13, atol=0)
            self.assertTrue(np.all(np.isfinite(hw.log_pdf(np.array([-10, 0, 10]), p))))

    def test_single_gaussian_minimum_and_eye_width(self):
        p = (0, 2, 0, 2, .3)
        rec = hw.metrics(p, (1e-6, 1e-12))
        expected_min = .5 * math.erfc(12.5 / (2 * SQRT2))
        self.assertAlmostEqual(rec["optimal_time_ps"], 25, delta=2e-5)
        self.assertAlmostEqual(rec["minimum_ber"] / expected_min, 1, delta=1e-10)
        eye = rec["eye_widths"]["1e-06"]
        left = bisect_root(lambda t: reference_ber(t, p) - 1e-6, 12.5, 25)
        self.assertTrue(eye["exists"])
        self.assertAlmostEqual(eye["left_ps"], left, delta=2e-9)
        self.assertAlmostEqual(eye["right_ps"], 50 - left, delta=2e-9)
        self.assertAlmostEqual(eye["width_ps"], 50 - 2 * left, delta=2e-9)
        self.assertAlmostEqual(eye["width_ui"], eye["width_ps"] / 25, delta=1e-14)
        closed = rec["eye_widths"]["1e-12"]
        self.assertFalse(closed["exists"])
        self.assertEqual((closed["width_ps"], closed["width_ui"]), (0., 0.))
        self.assertIsNone(closed["left_ps"])
        self.assertIsNone(closed["right_ps"])

    def test_positive_jitter_mean_shifts_bathtub_right(self):
        p0 = (0, 1.3, 0, 1.3, .4)
        shifted = (1.75, 1.3, 1.75, 1.3, .4)
        t = np.linspace(5, 45, 101)
        np.testing.assert_allclose(hw.model_ber(t, shifted), hw.model_ber(t - 1.75, p0), rtol=5e-14, atol=0)
        self.assertAlmostEqual(hw.metrics(shifted)["optimal_time_ps"], 26.75, delta=2e-5)

    def test_invalid_model_parameters_raise(self):
        invalid = ([0, 1, 0, 1], [0, 0, 0, 1, .5], [0, -1, 0, 1, .5],
                   [0, 1, 0, 0, .5], [0, 1, 0, 1, -.01], [0, 1, 0, 1, 1.01],
                   [math.nan, 1, 0, 1, .5], [0, math.inf, 0, 1, .5])
        for p in invalid:
            for fn in (hw.log_pdf, hw.log_tail, hw.model_ber):
                with self.subTest(parameters=p, function=fn.__name__):
                    with self.assertRaises(ValueError):
                        fn(0, p)

    def test_invalid_ber_controls_and_targets_raise(self):
        p = (0, 1, 0, 1, .5)
        for kwargs in ({"rho": 0}, {"rho": -1}, {"rho": 1.1}, {"rho": math.nan},
                       {"ti": math.nan}, {"tip1": math.inf}, {"ti": 40, "tip1": 30}):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    hw.model_ber(25, p, **kwargs)
        for target in (0, 1, -1, math.nan):
            with self.assertRaises(ValueError):
                hw.metrics(p, (target,))


class EmpiricalHistogramTests(unittest.TestCase):
    def test_cdf_sf_equal_direct_cell_mass_integration(self):
        centers = np.arange(-2., 3.)
        masses = np.array([1., 2., 0., 3., 0.]) / 6
        h = hw.Histogram(centers, masses)
        t = np.linspace(-5, 5, 81)
        left = [math.fsum(m * min(1., max(0., q - (x - .5))) for x, m in zip(centers, masses)) for q in t]
        right = [math.fsum(m * min(1., max(0., (x + .5) - q)) for x, m in zip(centers, masses)) for q in t]
        np.testing.assert_allclose(h.cdf(t), left, rtol=0, atol=4e-16)
        np.testing.assert_allclose(h.sf(t), right, rtol=0, atol=4e-16)
        np.testing.assert_allclose(h.cdf(t) + h.sf(t), 1, rtol=0, atol=4e-16)

    def test_reverse_survival_preserves_tiny_right_mass(self):
        h = hw.Histogram([-1, 0, 1], [.8, .2, 1e-18])
        self.assertEqual(float(h.sf(.5)), 1e-18)
        self.assertEqual(float(h.cdf(.5)), 1.)

    def test_normalization_pdf_units_and_true_zero_bins(self):
        h = hw.Histogram([-.25, -.125, 0, .125], [0, 2, 6, 0])
        self.assertEqual(h.dt, .125)
        np.testing.assert_allclose(h.pmf, [0, .25, .75, 0])
        np.testing.assert_allclose(h.pdf, [0, 2, 6, 0])
        self.assertEqual(float(np.sum(h.pdf) * h.dt), 1)
        self.assertEqual(float(h.ber(25)), 0)

    def test_sampled_bathtub_shift_and_boundary_orientation(self):
        centers = np.arange(-4., 5.)
        m0 = np.zeros(9); m0[4] = 1
        m2 = np.zeros(9); m2[6] = 1
        h0, h2 = hw.Histogram(centers, m0), hw.Histogram(centers, m2)
        t = np.linspace(0, 50, 401)
        np.testing.assert_allclose(h2.ber(t), h0.ber(t - 2), rtol=0, atol=1e-15)
        self.assertEqual(float(h2.ber(14)), .5)
        self.assertEqual(float(h2.ber(15)), 0)
        self.assertEqual(float(h2.ber(39)), 0)
        self.assertEqual(float(h2.ber(40)), .5)

    def test_invalid_histograms_raise(self):
        examples = (([0], [1]), ([0, 1], [1]), ([0, 0], [1, 1]),
                    ([1, 0], [1, 1]), ([0, 1, 3], [1, 1, 1]),
                    ([0, math.nan], [1, 1]), ([0, 1], [0, 0]),
                    ([0, 1], [-1, 2]), ([0, 1], [1, math.inf]),
                    ([[0, 1]], [[1, 1]]), ([0, 1], [1e308, 1e308]))
        for centers, masses in examples:
            with self.subTest(centers=centers, masses=masses):
                with np.errstate(over="ignore"):
                    with self.assertRaises(ValueError):
                        hw.Histogram(centers, masses)

    def test_invalid_or_insufficient_tail_window_raises(self):
        h = hw.Histogram([0, 1, 2], [0, 1, 0])
        for lo, hi in ((0, .1), (.1, .1), (.2, .1), (.1, .5), (math.nan, .1), (.001, .03)):
            with self.assertRaises(ValueError):
                hw.tail_masks(h, lo, hi)

    def test_legacy_display_offset_and_physical_grid_are_corrected(self):
        with tempfile.TemporaryDirectory(dir=BASE / "tests") as tmp:
            self.assertTrue(Path(tmp).resolve().is_relative_to((BASE / "tests").resolve()))
            path = Path(tmp) / "legacy.csv"
            with path.open("w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["jtx1", "jtx2", "jtx3", "jtx4", "jty1", "jty2", "jty3", "jty4", "tjp"])
                for i, x in enumerate(np.linspace(-12.5, 12.5, 200)):
                    mass = (2. if i == 90 else 1. if i == 110 else 0.) + 1e-15
                    writer.writerow([mass] * 8 + [x])
            groups, record = hw.load_legacy(path)
            self.assertEqual(len(groups), 8)
            for hist in groups.values():
                self.assertEqual(hist.dt, .125)
                self.assertEqual(hist.centers[-1], 12.375)
                self.assertEqual(np.count_nonzero(hist.pmf), 2)
                np.testing.assert_allclose(hist.pmf[[90, 110]], [2/3, 1/3], rtol=1e-15)
            self.assertEqual(record["subtracted_display_offset"], 1e-15)
            self.assertGreater(record["original_coordinate_spacing_ps"], .125)

    def test_known_continuous_mixture_is_recovered_from_bin_integrals(self):
        known = (-2., 1.2, 1.4, 1.1, .35)
        centers = np.arange(-12.5, 12.5, .05)
        masses = [mixture_tail(x + .025, known) - mixture_tail(x - .025, known) for x in centers]
        hist = hw.Histogram(centers, masses, 1e8)
        fit = hw.fit_tails(hist)
        np.testing.assert_allclose(fit["parameter_vector"], known, rtol=0, atol=2e-5)
        self.assertLess(fit["tail_log10_rmse"], 1e-8)
        self.assertFalse(fit["at_bound"])


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    print(stream.getvalue())
    record = {"passed": result.wasSuccessful(), "tests_run": result.testsRun,
              "failures": len(result.failures), "errors": len(result.errors),
              "skipped": len(result.skipped), "python": platform.python_version(),
              "runner": "python tests/test_hw13.py", "details": stream.getvalue()}
    (BASE / "results/test_results.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    raise SystemExit(0 if result.wasSuccessful() else 1)
