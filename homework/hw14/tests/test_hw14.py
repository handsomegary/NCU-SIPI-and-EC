"""Independent physical oracles for the HW14 NRZ/PAM4 solver.

Run with ``python -m unittest discover -s homework/hw14/tests -v``.
All fixtures are anonymous and synthetic.  Oracles use excited odd-mode
waves, analytic piecewise-linear ramps, and a serial feedback recurrence.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hw14_under_test", ROOT / "src" / "hw14.py")
hw14 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = hw14
SPEC.loader.exec_module(hw14)


def reference_ramps(levels, boundaries, dt, rise, count):
    """Exact superposition of each voltage increment's continuous ramp."""
    time = np.arange(count, dtype=float) * dt
    voltage = np.zeros(count, dtype=float)
    prior = 0.0
    for level, boundary in zip(levels, boundaries):
        voltage += (float(level) - prior) * np.clip((time - boundary) / rise, 0, 1)
        prior = float(level)
    return voltage


class TouchstoneAndSourcePlaneTests(unittest.TestCase):
    def test_four_port_reader_preserves_nonreciprocal_row_order(self):
        matrix = np.arange(1, 17, dtype=float).reshape(4, 4) / 20
        pairs = " ".join(f"{v:.8g} 0" for v in matrix.ravel())
        fixture = "! anonymous nonreciprocal fixture\n# MHz S RI R 50\n"
        fixture += f"0 {pairs}\n10 {pairs} ! inline comment\n"
        with tempfile.TemporaryDirectory(dir=ROOT / "tests") as directory:
            path = Path(directory) / "anonymous.s4p"
            path.write_text(fixture, encoding="ascii")
            frequency, loaded, z0 = hw14.read_touchstone(path)
        assert_array_equal(frequency, [0, 10e6])
        assert_allclose(loaded[0], matrix)
        self.assertEqual(z0, 50)
        self.assertEqual(loaded[0, 1, 0], 0.25)

    def test_mixed_mode_matches_direct_odd_mode_excitation(self):
        rng = np.random.default_rng(1401)
        scattering = rng.normal(size=(9, 4, 4)) + 1j * rng.normal(size=(9, 4, 4))
        incident = np.array([1, 0, -1, 0], dtype=float) / np.sqrt(2)
        outgoing = scattering @ incident
        expected = (outgoing[:, 1] - outgoing[:, 3]) / np.sqrt(2)
        assert_allclose(hw14.differential_forward(scattering), expected, atol=2e-15)

    def test_ideal_through_odd_mode_and_matched_source_divider(self):
        scattering = np.zeros((3, 4, 4), complex)
        scattering[:, 1, 0] = scattering[:, 0, 1] = 1
        scattering[:, 3, 2] = scattering[:, 2, 3] = 1
        odd_gain = hw14.differential_forward(scattering)
        assert_allclose(odd_gain, 1, atol=1e-15)
        # A 1 Vpp open-circuit source has +/-0.5 V source levels.  A matched
        # ideal-through channel receives +/-0.25 V (0.5 Vpp), not +/-0.5 V.
        load_levels = odd_gain[0] / 2 * np.array([-0.5, 0.5])
        assert_allclose(load_levels, [-0.25, 0.25], atol=1e-15)

    def test_actual_model_applies_the_matched_divider_once(self):
        # Complete band for p=16 at 20 Gbps: dt=3.125 ps, Fs=320 GHz.
        # Stub only the numerical optimizer; source/channel construction
        # remains the actual production implementation.
        frequency = np.arange(16001) * 10e6
        through = np.ones(len(frequency), complex)
        def fixed_ffe(pulse, clock, samples_per_ui, modulation, *unused):
            main = float(pulse[clock - samples_per_ui])
            return np.array([0., 1., 0.]), main, main, 0.
        with patch.object(hw14, "design_ffe", side_effect=fixed_ffe):
            model = hw14.design_equalizers(frequency, through, 20, "nrz", 16, 40)
        expected = np.zeros_like(model["h_channel"])
        expected[0] = 0.5
        assert_allclose(model["h_channel"], expected, atol=1e-15)
        source = hw14.source_waveform(np.array([0.5, 0.5]), np.array([0., 50e-12]),
                                     model["dt"], 10e-12, 100e-12)
        actual_load = np.convolve(source, model["h_channel"])[:len(source)]
        assert_allclose(actual_load[4:], 0.25, atol=1e-15)


class TransferAndImpulseTests(unittest.TestCase):
    def test_full_band_integer_delay_is_a_discrete_half_gain_impulse(self):
        dt = 1e-12
        frequency = np.arange(17) / (32 * dt)
        transfer = 0.5 * np.exp(-2j * np.pi * frequency * 7 * dt)
        kernel, full = hw14.impulse_response(frequency, transfer, dt, 24 * dt)
        expected = np.zeros(24)
        expected[7] = 0.5
        assert_allclose(kernel, expected, atol=1e-15)
        assert_allclose(np.fft.rfft(full), transfer, atol=1e-15)
        self.assertAlmostEqual(float(np.sum(full)), 0.5, places=14)
        # No extra dt or second factor 1/2 belongs in sample convolution.
        pulse = hw14.source_pulse(8, dt, rise_s=2 * dt, vpp=1)
        output = np.convolve(pulse, kernel)
        assert_allclose(output[7:7 + len(pulse)], pulse / 2, atol=1e-15)
        self.assertAlmostEqual(float(np.max(output)), 0.25, places=14)

    def test_input_band_is_preserved_and_unmeasured_band_is_zero(self):
        dt = 1e-12
        frequency = np.arange(9) / (64 * dt)
        transfer = np.exp(-2j * np.pi * frequency * 11 * dt) / (1 + 1j * frequency / 30e9)
        _, full = hw14.impulse_response(frequency, transfer, dt, 40 * dt)
        reconstructed = np.fft.rfft(full)
        assert_allclose(reconstructed[:len(frequency)], transfer, atol=3e-16)
        assert_allclose(reconstructed[len(frequency):], 0, atol=3e-16)

    def test_idft_period_and_bandwidth_constraints_are_enforced(self):
        frequency = np.arange(11) * 1e9
        with self.assertRaises(ValueError):
            hw14.impulse_response(frequency, np.ones(11), 1e-12, 1e-9)
        with self.assertRaises(ValueError):
            hw14.impulse_response(frequency, np.ones(11), 100e-12, 0.5e-9)


class GrayCodingTests(unittest.TestCase):
    def test_all_four_bit_pairs_have_correct_order_and_roundtrip(self):
        bits = np.array([0, 0, 0, 1, 1, 1, 1, 0], np.uint8)
        levels = hw14.gray_encode(bits)
        assert_allclose(levels, [-1, -1 / 3, 1 / 3, 1], atol=1e-15)
        assert_array_equal(hw14.gray_decode(levels), bits)
        pairs = bits.reshape(-1, 2)
        assert_array_equal(np.count_nonzero(np.diff(pairs.astype(int), axis=0), axis=1), [1, 1, 1])

    def test_random_payload_roundtrip_is_bit_exact(self):
        bits = np.random.default_rng(1402).integers(0, 2, size=2048, dtype=np.uint8)
        assert_array_equal(hw14.gray_decode(hw14.gray_encode(bits)), bits)

    def test_adjacent_symbol_error_changes_one_bit(self):
        truth = hw14.gray_decode(np.array([-1, -1 / 3, 1 / 3]))
        mistaken = hw14.gray_decode(np.array([-1 / 3, 1 / 3, 1]))
        self.assertEqual(np.count_nonzero(truth != mistaken), 3)


class TransmitterClockTests(unittest.TestCase):
    def test_continuous_ramps_match_independent_superposition(self):
        levels = np.array([-0.5, -1 / 6, 0.5, 1 / 6, -0.5])
        boundaries = np.array([-1.3, 21.1, 39.4, 60.8, 80.2]) * 1e-12
        dt, rise, stop = 0.5e-12, 5e-12, 100e-12
        actual = hw14.source_waveform(levels, boundaries, dt, rise, stop)
        assert_allclose(actual, reference_ramps(levels, boundaries, dt, rise, len(actual)), atol=2e-13)
        self.assertEqual(len(actual), 200)

    def test_positive_jitter_delays_transition(self):
        actual = hw14.source_waveform(np.array([0, 1]), np.array([0, 12e-12]),
                                     1e-12, 4e-12, 20e-12)
        self.assertAlmostEqual(actual[12], 0, places=13)
        self.assertAlmostEqual(actual[14], 0.5, places=13)
        self.assertAlmostEqual(actual[16], 1, places=13)

    def test_ffe_symbol_levels_are_shaped_on_current_jittered_clock(self):
        # This is a transmitter-level convention check.  The source receives
        # already filtered DAC levels, so all FFE tap contributions share
        # the current boundary timing instead of delayed independent clocks.
        symbols = np.array([-1, 1, -1, -1, 1], float)
        coefficients = np.array([-0.1, 0.8, 0.1])
        levels = np.convolve(symbols, coefficients)[:len(symbols)] / 2
        boundaries = np.arange(len(symbols)) * 25e-12 + np.array([0, 1.2, -0.7, 0.3, -0.4]) * 1e-12
        actual = hw14.source_waveform(levels, boundaries, 0.5e-12, 10e-12, 125e-12)
        expected = reference_ramps(levels, boundaries, 0.5e-12, 10e-12, len(actual))
        assert_allclose(actual, expected, atol=2e-13)
        self.assertLessEqual(np.max(np.abs(actual)), 0.5 + 1e-14)

    def test_overlapping_true_transitions_are_rejected(self):
        with self.assertRaises(ValueError):
            hw14.source_waveform(np.array([-0.5, 0.5, -0.5]),
                                 np.array([0, 4, 25]) * 1e-12, 1e-12, 10e-12, 40e-12)

    def test_jitter_uses_absolute_ps_and_shared_underlying_bit_clock(self):
        times, rj, sj = hw14.jittered_boundaries(2048, 25e-12, seed=1403, bit_stride=1)
        p_times, p_rj, p_sj = hw14.jittered_boundaries(1024, 50e-12, seed=1403, bit_stride=2)
        assert_allclose(times, np.arange(2048) * 25e-12 + rj + sj, atol=1e-23)
        assert_allclose(sj, 1.5e-12 * np.sin(2 * np.pi * 10e6 * np.arange(2048) * 25e-12), atol=1e-23)
        assert_allclose(p_rj, rj[::2], atol=1e-25)
        assert_allclose(p_sj, sj[::2], atol=1e-25)
        assert_allclose(p_times, times[::2], atol=1e-23)


class ContinuousCTLETests(unittest.TestCase):
    def test_actual_continuous_peak_is_the_declared_cap(self):
        zero, pole1, pole2, peak_db = 1e9, 15e9, 40e9, 18
        objective = lambda log_f: -abs(hw14.ctle_response(np.array([10 ** log_f]), zero, pole1, pole2, peak_db)[0])
        result = minimize_scalar(objective, bounds=(4, 13), method="bounded", options={"xatol": 1e-12})
        self.assertAlmostEqual(20 * np.log10(-result.fun), peak_db, places=8)
        peak_frequency = hw14.ctle_peak_frequency(zero, pole1, pole2)
        self.assertLess(abs(10 ** result.x - peak_frequency) / peak_frequency, 2e-6)

    def test_ctle_has_real_dc_and_real_impulse_symmetry(self):
        frequency = np.array([0, 1e9, 10e9, 30e9])
        positive = hw14.ctle_response(frequency, 2e9, 20e9, 40e9)
        negative = hw14.ctle_response(-frequency, 2e9, 20e9, 40e9)
        assert_allclose(negative, positive.conj(), atol=2e-15)
        self.assertAlmostEqual(positive[0].imag, 0, places=14)
        self.assertGreater(positive[0].real, 0)


class RecursiveDFETests(unittest.TestCase):
    def check_recovery(self, modulation, truth):
        truth = np.asarray(truth, dtype=float)
        main = 0.27
        taps = np.array([0.041, -0.023, 0.011])
        received = np.convolve(truth, np.r_[main, taps])[:len(truth)]
        decisions, feedback, corrected = hw14.dfe_decisions(received, taps, main, modulation)
        assert_allclose(decisions, truth, atol=1e-15)
        assert_allclose(corrected, main * truth, atol=1e-15)
        assert_allclose(received - feedback, corrected, atol=1e-15)
        # The second feedback must follow the first recovered decision when
        # an observation is deliberately forced across the wrong threshold.
        changed = received.copy()
        changed[0] = main if truth[0] < 0 else -main
        changed_decisions, changed_feedback, _ = hw14.dfe_decisions(changed, taps, main, modulation)
        self.assertNotEqual(changed_decisions[0], truth[0])
        self.assertAlmostEqual(changed_feedback[1], taps[0] * changed_decisions[0], places=14)

    def test_nrz_feedback_is_recursive(self):
        self.check_recovery("nrz", [-1, 1, 1, -1, -1, 1, -1, 1, 1])

    def test_pam4_feedback_is_recursive_and_preserves_inner_levels(self):
        self.check_recovery("pam4", [-1, -1 / 3, 1 / 3, 1, -1, 1 / 3, -1 / 3, 1, 1])

    def test_pam4_thresholds_are_midpoints_at_the_main_cursor_scale(self):
        main = 0.3
        samples = main * np.array([-0.9, -0.5, -0.01, 0.01, 0.5, 0.9])
        decisions, feedback, corrected = hw14.dfe_decisions(samples, [], main, "pam4")
        assert_allclose(decisions, [-1, -1 / 3, -1 / 3, 1 / 3, 1 / 3, 1])
        assert_array_equal(feedback, np.zeros(len(samples)))
        assert_allclose(corrected, samples)


class EyeAndPayloadAlignmentTests(unittest.TestCase):
    def test_vertical_gap_above_fixed_threshold_has_zero_width(self):
        truth = np.tile(np.array([-1., 1.]), 20)
        p, clock = 8, 16
        position = np.arange(clock + len(truth) * p)
        current = np.clip((position - clock + p // 2) // p, 0, len(truth) - 1)
        received = 2. + .3 * truth[current]
        metrics, _ = hw14.measure_eye(received, truth, clock, p, 1e-12, 4, 30,
                                     main_cursor_v=.3)
        self.assertAlmostEqual(metrics["worst_eye_height_v"], .6, places=14)
        self.assertEqual(metrics["worst_eye_width_ps"], 0)

    def test_interpolated_width_has_known_triangle_crossings(self):
        width, left, right = hw14.contiguous_width(
            np.array([-2., -1., 0., 1., 2.]), np.array([-1., 1., 2., 1., -1.]), 2)
        assert_allclose([width, left, right], [3., -1.5, 1.5], atol=1e-15)
        self.assertEqual(hw14.contiguous_width(np.arange(3.), np.array([1., -0.1, 1.]), 1)[0], 0)

    def test_eye_feedback_hold_uses_next_feedback_at_right_endpoint(self):
        rng = np.random.default_rng(1404)
        truth = rng.choice([-1., -1 / 3, 1 / 3, 1.], size=120)
        feedback = rng.normal(0, 0.04, size=len(truth))
        p, clock, main, dt = 8, 16, 0.3, 1e-12
        positions = np.arange(clock + (len(truth) - 1) * p + p)
        held_index = np.clip((positions - clock + p // 2) // p, 0, len(truth) - 1)
        received = main * truth[held_index] + feedback[held_index]
        start, count = 10, 100
        metrics, arrays = hw14.measure_eye(received, truth, clock, p, dt, start, count,
                                          feedback_v=feedback, main_cursor_v=main)
        assert_allclose(metrics["eye_heights_v"], [0.2, 0.2, 0.2], atol=3e-16)
        assert_array_equal(arrays["histogram"].sum(axis=0), np.full(p + 1, count))
        selected_truth = truth[start:start + count]
        endpoint = main * truth[start + 1:start + count + 1]
        for row, level in enumerate(arrays["levels"]):
            assert_allclose(arrays["level_min_v"][row, -1], np.min(endpoint[selected_truth == level]), atol=1e-15)
            assert_allclose(arrays["level_max_v"][row, -1], np.max(endpoint[selected_truth == level]), atol=1e-15)

    def check_synthetic_pipeline(self, modulation):
        # Known three-cursor link at the fixed clock.  With FFE [.1,.8,.1]
        # and ideal matched gain .5, cursors are [.025,.2,.025] volts.
        # Jittered transitions remain far from this plateau sampling time.
        p, rate = 16, 20.
        stride = 1 if modulation == "nrz" else 2
        ui = stride / (rate * 1e9)
        dt = ui / p
        delay = 3 * p
        # Trailing zeros retain sufficient convolution output for the final
        # guard decisions after the source's finite observation interval.
        kernel = np.zeros(delay + 4 * p + 1)
        kernel[delay] = 0.5
        pulse = hw14.source_pulse(p, dt)
        sbr = np.convolve(kernel, pulse)
        ffe = np.array([.1, .8, .1])
        main = .2
        model = {"modulation": modulation, "rate_gbps": rate, "samples_per_ui": p,
                 "ui": ui, "dt": dt, "impulse_ns": 1., "h_channel": kernel,
                 "h_eq": kernel, "source_pulse": pulse, "sbr_channel": sbr,
                 "sbr_eq": sbr, "sbr_ffe": hw14.apply_ffe_pulse(sbr, ffe, p),
                 "ffe": ffe, "dfe": np.array([.025, 0., 0., 0., 0.]),
                 "main": main, "clock": delay + p + p // 2,
                 "thresholds": np.array([0.]) if stride == 1 else main * np.array([-2/3, 0., 2/3]),
                 "ctle": np.ones(2)}
        metrics, archive, _ = hw14.simulate_record(model, bits=1024, seed=1405)
        start, count = int(archive["retained_start_symbol"]), int(archive["retained_symbol_count"])
        idx = np.arange(start, start + count)
        symbols = archive["symbols"]
        expected_input = .025 * symbols[idx + 1] + .2 * symbols[idx] + .025 * symbols[idx - 1]
        assert_allclose(archive["decision_input_v"][idx], expected_input, atol=8e-15)
        assert_allclose(archive["feedback_v"][idx], .025 * symbols[idx - 1], atol=1e-15)
        assert_allclose(archive["decision_corrected_v"][idx], .2 * symbols[idx] + .025 * symbols[idx + 1], atol=8e-15)
        assert_allclose(archive["decisions"][idx], symbols[idx], atol=1e-15)
        # The payload extraction must discard exactly the bit guards and
        # retain the original bits, including Gray pair boundaries.
        decoded = (archive["decisions"][idx] > 0).astype(np.uint8) if stride == 1 else hw14.gray_decode(archive["decisions"][idx])
        assert_array_equal(decoded, archive["bits"][start * stride:(start + count) * stride])
        self.assertEqual(metrics["bit_errors"], 0)
        self.assertEqual(metrics["symbol_errors"], 0)
        self.assertEqual(metrics["measured_bits"], 1024)
        self.assertEqual(metrics["measured_symbols"], 1024 // stride)
        assert_array_equal(archive["histogram"].sum(axis=0), np.full(p + 1, 1024 // stride))
        expected_height = .35 if stride == 1 else 2 * .2 / 3 - .05
        assert_allclose(metrics["eye_heights_v"], np.full(stride + (stride == 2), expected_height), atol=2e-14)

    def test_nrz_actual_pipeline_payload_clock_and_dfe_alignment(self):
        self.check_synthetic_pipeline("nrz")

    def test_pam4_actual_pipeline_payload_clock_and_dfe_alignment(self):
        self.check_synthetic_pipeline("pam4")


if __name__ == "__main__":
    unittest.main()
