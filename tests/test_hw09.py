"""Independent mathematical tests for HW9's single-ended DFE.

Run at the repository root: python -m unittest discover -s tests -p test_hw09.py -v
The feedback oracle uses a scalar history recurrence and clock-age indexing;
it does not reproduce the production future-waveform slice updates.
"""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
from pathlib import Path
import sys
import unittest

import numpy as np
from numpy.testing import assert_allclose


ROOT = Path(__file__).resolve().parents[1]
HW9 = ROOT / "homework" / "hw09"
sys.path.insert(0, str(HW9))
SPEC = importlib.util.spec_from_file_location("hw09_under_test", HW9 / "analyze.py")
dfe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dfe)
CHANNEL_SPEC = importlib.util.spec_from_file_location("hw09_channel_under_test", HW9 / "channel.py")
channel = importlib.util.module_from_spec(CHANNEL_SPEC)
CHANNEL_SPEC.loader.exec_module(channel)


def sampled_oracle(raw_samples, taps, threshold):
    """Compute each summer voltage from raw input and prior binary decisions."""
    decision = []
    corrected = []
    for value in raw_samples:
        voltage = float(value)
        for age, coefficient in enumerate(taps, start=1):
            if len(decision) >= age:
                voltage += coefficient * decision[-age]
        corrected.append(voltage)
        decision.append(int(voltage >= threshold))
    return np.asarray(corrected), np.asarray(decision, dtype=np.int8)


def held_oracle(raw, samples_per_ui, taps, threshold, phase):
    """Reconstruct held DAC feedback from decision age at each waveform sample."""
    raw = np.asarray(raw, dtype=float)
    phase %= samples_per_ui
    positions = np.arange(phase, len(raw), samples_per_ui, dtype=int)
    corrected, decisions = sampled_oracle(raw[positions], taps, threshold)
    wave = raw.copy()
    for time_index in range(len(raw)):
        latest = (time_index - phase - samples_per_ui // 2) // samples_per_ui
        for delay, coefficient in enumerate(taps):
            decision_index = latest - delay
            if 0 <= decision_index < len(decisions):
                wave[time_index] += coefficient * decisions[decision_index]
    return wave, positions, decisions, corrected


def exhaustive_bounds(sbr, samples_per_ui, latency, ntaps):
    """Enumerate binary states of all cursors, then cancel the selected history."""
    offsets = [j for j in range(-latency // samples_per_ui - 1,
                                len(sbr) // samples_per_ui + 1)
               if 0 <= latency + j * samples_per_ui < len(sbr) and j != 0]
    cursors = [float(sbr[latency + j * samples_per_ui]) for j in offsets]
    low, high = [], []
    for state in itertools.product((0, 1), repeat=len(offsets)):
        raw = sum(value * bit for value, bit in zip(cursors, state))
        correction = -sum(value * bit for j, value, bit in zip(offsets, cursors, state)
                          if 1 <= j <= ntaps)
        low.append(raw + correction)
        high.append(float(sbr[latency]) + raw + correction)
    return max(low), min(high)


class SampledDecisionTests(unittest.TestCase):
    def test_threshold_equality_is_logic_one_and_feedback_uses_past_only(self):
        corrected, decision = dfe.sample_reference([0.5, 0.4, 0.5], [-0.1], 0.5)
        assert_allclose(corrected, [0.5, 0.3, 0.5], atol=1e-14)
        assert_allclose(decision, [1, 0, 1], atol=0)

    def test_single_ended_taps_have_voltage_units_without_main_normalization(self):
        # Main cursor is 0.4 V; binary decisions still take values 0 and 1.
        truth = np.array([1, 1, 0, 1, 0, 0, 1, 1])
        raw = np.convolve(truth, [0.4, 0.13, -0.03])
        corrected, decisions = dfe.sample_reference(raw.tolist(), [-0.13, 0.03], 0.2)
        expected = np.r_[truth * 0.4, 0.0, 0.0]
        assert_allclose(corrected, expected, atol=1e-14)
        assert_allclose(decisions, np.r_[truth, 0, 0], atol=0)

    def test_independent_scalar_history_matches_mixed_sign_multitap_decisions(self):
        raw = np.random.default_rng(507).normal(0.25, 0.5, 109)
        taps = [-0.31, 0.17, -0.08]
        expected = sampled_oracle(raw, taps, 0.22)
        actual = dfe.sample_reference(raw, taps, 0.22)
        assert_allclose(actual[0], expected[0], atol=1e-14)
        assert_allclose(actual[1], expected[1], atol=0)

    def test_zero_taps_are_a_plain_slicer(self):
        corrected, decisions = dfe.sample_reference([-0.2, 0.49, 0.5, 1.2], [], 0.5)
        assert_allclose(corrected, [-0.2, 0.49, 0.5, 1.2])
        assert_allclose(decisions, [0, 0, 1, 1])

    def test_wrong_decision_propagates_to_the_next_symbol(self):
        clean = np.array([1.0, 0.6, 0.0, 0.0])
        _, correct = dfe.sample_reference(clean, [-0.6], 0.5)
        corrupted = clean.copy()
        corrupted[0] = 0.4
        summer, erroneous = dfe.sample_reference(corrupted, [-0.6], 0.5)
        assert_allclose(correct, [1, 0, 0, 0])
        assert_allclose(erroneous, [0, 1, 0, 0])
        assert_allclose(summer, [0.4, 0.6, -0.6, 0.0], atol=1e-14)


class HeldFeedbackTests(unittest.TestCase):
    def test_waveform_matches_independent_clock_age_reference_including_tail(self):
        raw = np.random.default_rng(201).normal(0.2, 0.5, 103)
        original = raw.copy()
        taps = [-0.31, 0.17, -0.08]
        expected, positions, decisions, samples = held_oracle(raw, 8, taps, 0.1, 3)
        actual, actual_positions, actual_decisions = dfe.fixed_rate_dfe(raw, 8, taps, 0.1, 3)
        assert_allclose(raw, original, atol=0)
        assert_allclose(actual_positions, positions, atol=0)
        assert_allclose(actual_decisions, decisions, atol=0)
        assert_allclose(actual, expected, atol=1e-14)
        assert_allclose(actual[positions], samples, atol=1e-14)

    def test_phase_is_a_modulo_clock_phase_not_an_absolute_start_latency(self):
        raw = np.random.default_rng(72).normal(size=69).tolist()
        expected = dfe.fixed_rate_dfe(raw, 8, [-0.2], 0.1, 3)
        for phase in (35, -5):
            with self.subTest(phase=phase):
                actual = dfe.fixed_rate_dfe(raw, 8, [-0.2], 0.1, phase)
                for left, right in zip(actual, expected):
                    assert_allclose(left, right, atol=0)

    def test_update_starts_half_ui_after_decision_and_does_not_change_current_bit(self):
        raw = np.zeros(29)
        raw[3] = 1
        actual, positions, decision = dfe.fixed_rate_dfe(raw, 8, [-0.25], 0.5, 3)
        assert_allclose(positions, [3, 11, 19, 27])
        assert_allclose(decision, [1, 0, 0, 0])
        expected = raw.copy()
        expected[7:15] -= 0.25
        assert_allclose(actual, expected, atol=0)
        self.assertEqual(actual[3], 1)
        self.assertEqual(actual[6], 0)
        self.assertEqual(actual[7], -0.25)
        self.assertEqual(actual[15], 0)

    def test_zero_taps_preserve_full_input_and_slice_the_selected_phase(self):
        raw = np.linspace(-0.3, 1.3, 57)
        actual, positions, decision = dfe.fixed_rate_dfe(raw.tolist(), 8, [], 0.5, 5)
        assert_allclose(actual, raw, atol=0)
        assert_allclose(positions, np.arange(5, len(raw), 8), atol=0)
        assert_allclose(decision, raw[positions] >= 0.5, atol=0)

    def test_delayed_channel_recovers_known_bits_with_full_latency_alignment(self):
        # A delay exceeding two UIs makes phase-only truth alignment incorrect.
        samples_per_ui, delay = 8, 19
        truth = np.random.default_rng(148).integers(0, 2, 97)
        truth[:5] = [1, 0, 1, 1, 0]
        h = np.zeros(delay + 2 * samples_per_ui + 1)
        h[[delay, delay + samples_per_ui, delay + 2 * samples_per_ui]] = [1.0, 0.7, -0.2]
        raw = np.convolve(np.repeat(truth, samples_per_ui), h)
        latency = delay + samples_per_ui // 2
        positions = np.arange(len(truth)) * samples_per_ui + latency
        post, clocks, decisions = dfe.fixed_rate_dfe(raw, samples_per_ui, [-0.7, 0.2], 0.5, latency)
        self.assertGreater(np.count_nonzero((raw[positions] >= 0.5) != truth), 0)
        assert_allclose(post[positions], truth, atol=2e-15)
        clock_indices = (positions - (latency % samples_per_ui)) // samples_per_ui
        assert_allclose(clocks[clock_indices], positions, atol=0)
        assert_allclose(decisions[clock_indices], truth, atol=0)


class CursorDesignTests(unittest.TestCase):
    def setUp(self):
        self.sbr = np.zeros(28)
        self.sbr[::4] = [0.03, -0.06, 1.0, 0.3, -0.12, 0.07, -0.02]

    def test_mixed_sign_residual_bounds_match_exhaustive_symbol_patterns(self):
        design = dfe.design_at(self.sbr, 4, 8, 2)
        low, high = exhaustive_bounds(self.sbr, 4, 8, 2)
        self.assertAlmostEqual(low, 0.10)
        self.assertAlmostEqual(high, 0.92)
        self.assertAlmostEqual(design["low_upper_v"], low)
        self.assertAlmostEqual(design["high_lower_v"], high)
        self.assertAlmostEqual(design["threshold_v"], 0.51)
        self.assertAlmostEqual(design["worst_case_eye_v"], 0.82)
        self.assertAlmostEqual(design["main_cursor_v"], 1.0)
        assert_allclose(design["taps_v"], [-0.3, 0.12], atol=0)
        self.assertEqual(design["latency_samples"], 8)
        self.assertEqual(design["phase_samples"], 0)

    def test_no_dfe_keeps_all_precursor_and_postcursor_isi_in_the_bounds(self):
        design = dfe.design_at(self.sbr, 4, 8, 0)
        low, high = exhaustive_bounds(self.sbr, 4, 8, 0)
        self.assertAlmostEqual(design["low_upper_v"], low)
        self.assertAlmostEqual(design["high_lower_v"], high)
        self.assertAlmostEqual(design["threshold_v"], (low + high) / 2)
        assert_allclose(design["taps_v"], [], atol=0)

    def test_closed_predicted_eye_is_preserved_without_claiming_guaranteed_recovery(self):
        sbr = np.zeros(24)
        sbr[[4, 8, 12, 16]] = [0.8, 1.0, 0.3, -0.6]
        design = dfe.design_at(sbr, 4, 8, 1)
        self.assertAlmostEqual(design["low_upper_v"], 0.8)
        self.assertAlmostEqual(design["high_lower_v"], 0.4)
        self.assertAlmostEqual(design["worst_case_eye_v"], -0.4)
        self.assertAlmostEqual(design["threshold_v"], 0.6)

    def test_phase_search_maximizes_residual_gap_instead_of_main_cursor_peak(self):
        sbr = np.zeros(24)
        sbr[6:10] = [0.9, 0.92, 1.0, 0.98]
        sbr[10:14] = [0.1, 0.08, 0.6, 0.05]
        sbr[14:18] = [0.55, 0.35, 0.7, 0.02]
        best, candidates = dfe.design_dfe(sbr, 4, 1)
        self.assertEqual(len(candidates), 4)
        self.assertEqual([item["latency_samples"] for item in candidates], [6, 7, 8, 9])
        self.assertEqual({item["phase_samples"] for item in candidates}, {0, 1, 2, 3})
        expected_gaps = []
        for candidate in candidates:
            low, high = exhaustive_bounds(sbr, 4, candidate["latency_samples"], 1)
            expected_gaps.append(high - low)
            self.assertAlmostEqual(candidate["worst_case_eye_v"], high - low)
            self.assertAlmostEqual(candidate["threshold_v"], (low + high) / 2)
        self.assertEqual(best["latency_samples"], 9)
        self.assertNotEqual(best["latency_samples"], int(np.argmax(sbr)))
        self.assertAlmostEqual(best["worst_case_eye_v"], max(expected_gaps))
        assert_allclose(best["taps_v"], [-0.05])

    def test_equal_objective_uses_first_candidate(self):
        sbr = np.zeros(28)
        sbr[8:12] = 1
        best, candidates = dfe.design_dfe(sbr, 4, 1)
        scores = [item["worst_case_eye_v"] for item in candidates]
        first = int(np.argmax(scores))
        self.assertEqual(best["latency_samples"], candidates[first]["latency_samples"])
        self.assertEqual(best["latency_samples"], 8)

    def test_postcursors_beyond_finite_response_are_zero(self):
        sbr = np.array([0, 0, 0, 0, 0.7, 0.1])
        design = dfe.design_at(sbr, 4, 4, 3)
        assert_allclose(design["taps_v"], [0, 0, 0], atol=0)
        self.assertAlmostEqual(design["threshold_v"], 0.35)
        self.assertAlmostEqual(design["worst_case_eye_v"], 0.7)


class ChannelReferenceTests(unittest.TestCase):
    def test_sectional_convolution_matches_direct_linear_convolution_across_boundaries(self):
        rng = np.random.default_rng(591)
        x, h = rng.normal(size=77), rng.normal(size=19)
        x[[0, 8, 16, -1]] += 2
        expected = np.convolve(x, h)
        for section in (1, 8, 33, 79):
            with self.subTest(section=section):
                actual = channel.sectional_convolve(x, h, section_samples=section)
                assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
                self.assertEqual(len(actual), len(x) + len(h) - 1)

    def test_matched_source_divider_is_applied_exactly_once(self):
        n, dt, delay = 64, 0.5e-12, 9
        frequency = np.fft.rfftfreq(n, dt)
        transfer = 0.5 * np.exp(-2j * np.pi * frequency * delay * dt)
        _, h = channel.frequency_impulse(frequency, transfer, dt, alpha=0)
        expected = np.zeros(n)
        expected[delay] = 0.5
        assert_allclose(h, expected, atol=2e-15)
        output = channel.sectional_convolve(np.full(91, 2.0), h, section_samples=13)
        assert_allclose(output[delay:delay + 91], 1.0, atol=2e-14)

    def test_causal_finite_ramp_pulse_matches_explicit_rise_and_fall(self):
        expected = [0, 1, 2, 2, 2, 2, 2, 2, 2, 1, 0]
        pulse = channel.trapezoid_pulse(8.0, 2.0, 1.0, vpp=2)
        assert_allclose(pulse, expected, atol=0)
        explicit = channel.levels_to_waveform([2, 0], 8, 2)
        assert_allclose(explicit[:len(expected)], expected, atol=0)

    def test_actual_touchstone_grid_dc_scaling_and_input_hash(self):
        candidates = [HW9 / "data" / "HW7_stripline_50cm.s2p",
                      HW9 / "instructions" / "HW7_stripline_50cm.s2p"]
        data = next((path for path in candidates if path.is_file()), candidates[0])
        self.assertEqual(hashlib.sha256(data.read_bytes()).hexdigest(),
                         "564560062e750405e89cadbb8be73cc7a40ad3188d19383ebbfe1e22723f5efc")
        frequency, scattering, z0 = channel.read_touchstone(data)
        self.assertEqual(len(frequency), 8001)
        self.assertEqual(frequency[0], 0)
        self.assertEqual(frequency[-1], 80e9)
        assert_allclose(np.diff(frequency), 10e6, rtol=0, atol=0)
        self.assertEqual(z0, 50)
        self.assertAlmostEqual(scattering[0, 1, 0].real, 0.969466574, places=10)
        self.assertAlmostEqual((scattering[0, 1, 0] / 2).real, 0.484733287, places=10)


if __name__ == "__main__":
    unittest.main()
