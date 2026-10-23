"""Boundary, counting, normalization, crossing, and timing-drift regression tests."""

import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hw5 import (
    EyeConfig, analyze, fixed_time_bin_indices, histogram,
    positive_zero_crossings, triggered_time_bin_indices,
    vertical_normalize, voltage_bin_indices,
)


class EyeTests(unittest.TestCase):
    def setUp(self):
        self.config = EyeConfig()

    def test_voltage_offset_and_negative_voltage(self):
        actual = voltage_bin_indices([-0.1, -0.095, -0.09, -0.001, 0, 0.015, 1.099, 1.1], self.config)
        np.testing.assert_array_equal(actual, [0, 0, 1, 9, 10, 11, 119, 119])

    def test_voltage_exact_decimal_edges(self):
        actual = voltage_bin_indices(np.linspace(-0.1, 1.1, 121), self.config)
        np.testing.assert_array_equal(actual, np.minimum(np.arange(121), 119))

    def test_voltage_bounds_reject_real_outside_values(self):
        for values in ([-0.1001], [1.1001], [np.nan]):
            with self.assertRaises(ValueError):
                voltage_bin_indices(values, self.config)

    def test_fixed_time_grid_does_not_lose_decimal_boundaries(self):
        timestamps_ps = (1e-9 + np.arange(250_001) * 1e-12) * 1e12
        actual = fixed_time_bin_indices(timestamps_ps, self.config)
        np.testing.assert_array_equal(actual, (1000 + np.arange(250_001)) % 50)

    def test_fixed_time_edges_and_origin(self):
        actual = fixed_time_bin_indices([-0.1, 0, 0.99999, 1, 49.99999, 50, 137], self.config)
        np.testing.assert_array_equal(actual, [49, 0, 0, 1, 49, 0, 37])
        np.testing.assert_array_equal(fixed_time_bin_indices([10, 11, 60], EyeConfig(origin_ps=10)), [0, 1, 0])

    def test_interpolated_crossing_and_zero_plateau(self):
        actual = positive_zero_crossings([0, 2, 4, 6, 8, 10], [-1, 3, 0, 0, -2, 0])
        np.testing.assert_allclose(actual, [0.5, 10], rtol=0, atol=1e-14)

    def test_clock_requires_strict_time(self):
        with self.assertRaises(ValueError):
            positive_zero_crossings([0, 0], [-1, 1])

    def test_trigger_window_has_no_wrap_or_duplicate_assignment(self):
        time = np.array([9, 10, 59.999, 60, 65, 110, 160])
        indices, valid = triggered_time_bin_indices(time, [10, 65], self.config)
        np.testing.assert_array_equal(valid, [False, True, True, False, True, True, False])
        np.testing.assert_array_equal(indices[valid], [0, 49, 0, 45])

    def test_newest_trigger_wins_in_overlapping_nominal_windows(self):
        indices, valid = triggered_time_bin_indices([49.5, 50, 99], [0, 49, 98], self.config)
        np.testing.assert_array_equal(indices, [0, 1, 1])
        self.assertTrue(np.all(valid))

    def test_empty_triggers_and_columns_are_safe(self):
        indices, valid = triggered_time_bin_indices([0, 1], [], self.config)
        np.testing.assert_array_equal(indices, [-1, -1])
        self.assertFalse(np.any(valid))
        probability, totals = vertical_normalize([[0, 1, 0], [0, 3, 0]])
        np.testing.assert_array_equal(totals, [0, 4, 0])
        np.testing.assert_allclose(probability, [[0, 0.25, 0], [0, 0.75, 0]])
        self.assertTrue(np.all(np.isfinite(probability)))

    def test_histogram_matches_hand_count(self):
        counts = histogram(np.array([0, 1, 1, 119]), np.array([0, 0, 0, 49]), self.config)
        self.assertEqual(int(counts.sum()), 4)
        self.assertEqual(counts[0, 0], 1)
        self.assertEqual(counts[1, 0], 2)
        self.assertEqual(counts[119, 49], 1)

    def test_histogram_rejects_invalid_indices(self):
        with self.assertRaises(ValueError):
            histogram([0], [50], self.config)

    def test_config_rejects_fractional_bins_and_nonfinite_values(self):
        for kwargs in ({"bins_per_ui": 2.5}, {"dv": 0.07}, {"ui_ps": 0}, {"vmin": np.nan}):
            with self.assertRaises(ValueError):
                EyeConfig(**kwargs)

    def test_synthetic_rate_drift_is_removed_by_trigger_alignment(self):
        # True period 51 ps versus a fixed 50 ps window creates deterministic
        # accumulated phase error, while the clock remains perfectly periodic.
        time_ps = np.arange(0, 10_200, 0.1)
        phase = np.sin(2 * np.pi * time_ps / 51)
        output = (phase >= 0).astype(float)
        data = np.column_stack((time_ps * 1e-12, output, phase, np.zeros(time_ps.size)))
        arrays, summary = analyze(data, self.config)
        occupied = arrays["triggered_column_counts"] > 0
        fixed_mixing = np.mean(1 - arrays["fixed_probability"].max(axis=0))
        triggered_mixing = np.mean(1 - arrays["triggered_probability"][:, occupied].max(axis=0))
        self.assertGreater(fixed_mixing, 0.45)
        self.assertLess(triggered_mixing, 0.02)
        self.assertEqual(summary["fixed_rate"]["total_samples"], time_ps.size)
        self.assertLessEqual(summary["triggered"]["total_samples"], time_ps.size)
        self.assertTrue(np.all(np.isfinite(arrays["triggered_probability"])))
        np.testing.assert_allclose(arrays["triggered_probability"][:, occupied].sum(axis=0), 1)


if __name__ == "__main__":
    unittest.main()
