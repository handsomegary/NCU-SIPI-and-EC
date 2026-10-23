"""Independent mathematical regression tests for the HW8 FFE simulation.

Run from the HW8 directory with ``python -m unittest discover -s tests -v``.
The oracles below use closed-form signals, handout numbers, or direct linear
convolution rather than a second copy of the production algorithms.
"""

from __future__ import annotations

import itertools
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
from numpy.testing import assert_allclose


ROOT = Path(__file__).resolve().parents[1]
if (ROOT / "src" / "hw8_ffe.py").is_file():
    MODULE = ROOT / "src" / "hw8_ffe.py"
    DATA = ROOT / "data" / "HW7_stripline_50cm.s2p"
else:
    MODULE = ROOT / "homework" / "hw08" / "analyze.py"
    DATA = ROOT / "homework" / "hw08" / "instructions" / "HW7_stripline_50cm.s2p"
SPEC = importlib.util.spec_from_file_location("hw08", MODULE)
ffe = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ffe)


class TouchstoneTests(unittest.TestCase):
    def read_text(self, text):
        with tempfile.TemporaryDirectory(dir=ROOT / "tests") as folder:
            path = Path(folder) / "synthetic.s2p"
            path.write_text(text, encoding="ascii")
            return ffe.read_touchstone(path)

    def test_ma_degrees_units_and_asymmetric_port_order(self):
        freq, s, z0 = self.read_text(
            "! Deliberately nonreciprocal fixture\n# GHz S MA R 75\n"
            "0 0.1 0 0.8 90 0.2 -90 0.3 180\n"
            "0.25 0.2 0 0.7 45 0.4 -45 0.5 0 ! inline comment\n"
        )
        assert_allclose(freq, [0.0, 250e6])
        self.assertEqual(s.shape, (2, 2, 2))
        self.assertEqual(z0, 75)
        assert_allclose(s[0], [[0.1, -0.2j], [0.8j, -0.3]], atol=1e-15)
        assert_allclose(s[1, 1, 0], 0.7 * np.exp(1j * np.pi / 4))

    def test_ri_multiline_records_and_scientific_notation(self):
        freq, s, _ = self.read_text(
            "# MHz S RI R 50\n"
            "0 1 0 0.5 -0.25\n 0.125 0.2 0.25 -0.5\n"
            "1e1 2 0 0.7 -0.1 0.3 0.4 0.6 -0.2\n"
        )
        assert_allclose(freq, [0.0, 10e6])
        assert_allclose(s[0], [[1.0, 0.125 + 0.2j], [0.5 - 0.25j, 0.25 - 0.5j]])

    def test_db_magnitude_and_khz_units(self):
        freq, s, _ = self.read_text(
            "# kHz S DB R 50\n0 -20 0 -6.020599913279624 180 -40 90 0 0\n"
            "100 -20 0 -6.020599913279624 180 -40 90 0 0\n"
        )
        assert_allclose(freq, [0.0, 100e3])
        assert_allclose(s[0], [[0.1, 0.01j], [-0.5, 1.0]], atol=1e-14)

    def test_actual_channel_preserves_known_dc_s21(self):
        freq, s, z0 = ffe.read_touchstone(DATA)
        self.assertEqual(freq[0], 0.0)
        assert_allclose(np.diff(freq), 10e6)
        self.assertEqual(z0, 50)
        self.assertAlmostEqual(s[0, 1, 0].real, 0.969466574, places=10)
        self.assertAlmostEqual((s[0, 1, 0] / 2).real, 0.484733287, places=10)

    def test_malformed_record_rejected(self):
        with self.assertRaises(ValueError):
            self.read_text("# Hz S RI R 50\n0 1 0 2 0 3 0\n")


class TransformTests(unittest.TestCase):
    def setUp(self):
        self.n = 32
        self.dt = 0.125
        self.freq = np.fft.rfftfreq(self.n, self.dt)

    def impulse(self, spectrum):
        return ffe.frequency_impulse(self.freq, spectrum, self.dt, alpha=0)

    def test_identity_channel_has_unit_discrete_impulse(self):
        time, impulse = self.impulse(np.ones(self.freq.size, dtype=complex))
        expected = np.zeros(self.n)
        expected[0] = 1
        assert_allclose(time, np.arange(self.n) * self.dt)
        assert_allclose(impulse, expected, atol=1e-14)

    def test_dc_bin_is_not_doubled(self):
        spectrum = np.zeros(self.freq.size, dtype=complex)
        spectrum[0] = 1.5
        _, impulse = self.impulse(spectrum)
        assert_allclose(impulse, np.full(self.n, 1.5 / self.n), atol=1e-14)
        self.assertAlmostEqual(np.sum(impulse), 1.5)

    def test_nyquist_bin_is_not_doubled(self):
        spectrum = np.zeros(self.freq.size, dtype=complex)
        spectrum[-1] = 0.75
        _, impulse = self.impulse(spectrum)
        assert_allclose(impulse, 0.75 * (-1.0) ** np.arange(self.n) / self.n, atol=1e-14)

    def test_positive_interior_bin_has_hermitian_partner(self):
        spectrum = np.zeros(self.freq.size, dtype=complex)
        spectrum[3] = 0.6 * np.exp(0.3j)
        _, impulse = self.impulse(spectrum)
        expected = 1.2 * np.cos(2 * np.pi * 3 * np.arange(self.n) / self.n + 0.3) / self.n
        assert_allclose(impulse, expected, atol=1e-14)

    def test_bandlimited_input_is_zero_padded_without_extra_gain(self):
        spectrum = np.array([0.5, 0.25 + 0.1j, 0.2 - 0.3j])
        _, impulse = ffe.frequency_impulse(self.freq[:3], spectrum, self.dt, alpha=0)
        time_index = np.arange(self.n)
        expected = np.full(self.n, 0.5 / self.n)
        expected += 2 * np.real(spectrum[1] * np.exp(2j * np.pi * time_index / self.n)) / self.n
        expected += 2 * np.real(spectrum[2] * np.exp(4j * np.pi * time_index / self.n)) / self.n
        assert_allclose(impulse, expected, atol=1e-14)

    def test_integer_delay_and_matched_source_divider_once(self):
        delay = 7
        h = 0.5 * np.exp(-2j * np.pi * self.freq * delay * self.dt)
        _, impulse = self.impulse(h)
        expected = np.zeros(self.n)
        expected[delay] = 0.5
        assert_allclose(impulse, expected, atol=1e-14)
        # A 0..2 V source through an ideal matched S21=1 channel reaches 1 V.
        x = np.full(70, 2.0)
        y = ffe.sectional_convolve(x, impulse, section_samples=11)
        assert_allclose(y[delay : delay + 70], 1.0, atol=1e-13)

    def test_transform_requires_dc_and_uniform_frequency_grid(self):
        with self.assertRaises(ValueError):
            ffe.frequency_impulse(self.freq[1:], np.ones(self.freq.size - 1), self.dt)
        bad_frequency = self.freq.copy()
        bad_frequency[5] += 0.03
        with self.assertRaises(ValueError):
            ffe.frequency_impulse(bad_frequency, np.ones(self.freq.size), self.dt)


class PulseAndWaveformTests(unittest.TestCase):
    def test_analytic_sbr_has_causal_trapezoid_phase_and_area(self):
        dt = 0.025e-12
        n = 4096
        tb, tr, vpp = 25e-12, 5e-12, 2.0
        freq = np.fft.rfftfreq(n, dt)
        time, pulse = ffe.get_sbr(freq, np.ones(freq.size), tb, tr, dt, vpp=vpp, alpha=0)
        expected = np.maximum(0.0, np.minimum.reduce([time / tr, np.ones(time.size), (tb + tr - time) / tr])) * vpp
        # Fourier reconstruction of a piecewise-linear pulse has a small error
        # at slope discontinuities. 0.5% is below one 200-sample ramp step.
        assert_allclose(pulse, expected, atol=0.005, rtol=0)
        self.assertAlmostEqual(np.sum(pulse) * dt / (vpp * tb), 1.0, places=10)
        self.assertLess(np.max(np.abs(pulse[time > tb + tr + 4 * dt])), 0.005)

    def test_sampled_single_bit_matches_explicit_waveform(self):
        tb, tr, dt = 25e-12, 5e-12, 0.5e-12
        pulse = ffe.trapezoid_pulse(tb, tr, dt, vpp=2)
        waveform = ffe.levels_to_waveform([2.0, 0.0], 50, 10)
        common = min(len(pulse), len(waveform))
        assert_allclose(pulse[:common], waveform[:common], atol=1e-14)
        self.assertAlmostEqual(np.sum(pulse) * dt, 2 * tb, delta=1e-24)

    def test_pulse_rejects_nonintegral_ui_and_risetime(self):
        with self.assertRaises(ValueError):
            ffe.trapezoid_pulse(25.2e-12, 5e-12, 0.5e-12)
        with self.assertRaises(ValueError):
            ffe.trapezoid_pulse(25e-12, 5.2e-12, 0.5e-12)

    def test_fixed_risetime_for_small_and_large_transitions(self):
        n_ui, n_rise = 20, 5
        levels = [2.0, -0.2, 0.3, 1.8]
        waveform = ffe.levels_to_waveform(levels, n_ui, n_rise, previous=0)
        for k, value in enumerate(levels):
            previous = 0.0 if k == 0 else levels[k - 1]
            segment = waveform[k * n_ui : (k + 1) * n_ui]
            with self.subTest(transition=k):
                assert_allclose(segment[:n_rise], previous + (value - previous) * np.arange(n_rise) / n_rise)
                assert_allclose(segment[n_rise:], value)

    def test_zero_rise_is_ideal_nrz_and_source_is_unipolar(self):
        levels = ffe.ffe_levels([0, 1, 1, 0], [1])
        assert_allclose(levels, [0, 2, 2, 0])
        waveform = ffe.levels_to_waveform(levels, 3, 0)
        assert_allclose(waveform, np.repeat([0, 2, 2, 0], 3))

    def test_handout_two_tap_sequence_and_fir_tail(self):
        bits = [0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0]
        levels = ffe.ffe_levels(bits, [0.8, -0.2], vpp=1)
        assert_allclose(levels[:5], [0, 0.8, 0.6, -0.2, 0])
        self.assertEqual(len(levels), len(bits) + 1)
        assert_allclose(ffe.ffe_levels([0, 1], [0.8, -0.2], vpp=1), [0, 0.8, -0.2])

    def test_tap_scale_independence_and_maximum_swing(self):
        taps = np.array([-0.13, 0.61, -0.21, 0.05])
        states = np.array(list(itertools.product([0, 1], repeat=len(taps))))
        levels = []
        for bits in states:
            response = ffe.ffe_levels(bits, taps * 17, vpp=2)
            levels.append(response[len(taps) - 1])
        self.assertAlmostEqual(max(levels) - min(levels), 2.0, places=13)
        assert_allclose(ffe.ffe_levels([0, 1, 1, 0], taps), ffe.ffe_levels([0, 1, 1, 0], 17 * taps))


class ZeroForcingTests(unittest.TestCase):
    def test_no_equalizer_is_exactly_one(self):
        assert_allclose(ffe.zf_taps([0.33], 0), [1.0])

    def test_lecture_page_21_numerical_example(self):
        cursors = np.array([36, 230, 97, 37, 18], dtype=float)
        taps = ffe.zf_taps(cursors, 1)
        # Solve the numbers in the displayed handout matrix. The lecture's
        # printed rounded coefficients differ slightly from that matrix.
        assert_allclose(taps, [-0.0975361791612765, 0.6231478113081557,
                               -0.2494519337422591, 0.0149290470221032,
                               -0.0149350287662054], atol=1e-14, rtol=0)
        self.assertAlmostEqual(np.sum(np.abs(taps)), 1.0, places=13)
        output = np.convolve(cursors, taps)
        self.assertLess(np.max(np.abs(output[[1, 3, 4, 5]])), 1e-12)
        self.assertGreater(output[2], 100)

    def test_general_two_precursor_indexing_by_full_convolution(self):
        # Unequal precursor samples and nonzero postcursors prevent accidental
        # agreement with the broken special-case indexing for npre <= 1.
        cursors = np.array([0.03, 0.09, 1.0, 0.31, 0.07, 0.025])
        npre = 2
        taps = ffe.zf_taps(cursors, npre)
        selected = np.convolve(cursors, taps)[npre : npre + cursors.size]
        self.assertAlmostEqual(np.sum(abs(taps)), 1.0, places=13)
        self.assertGreater(selected[npre], 0)
        assert_allclose(np.delete(selected, npre), 0.0, atol=1e-14)

    def test_cursor_sampling_zero_pads_both_ends(self):
        sbr = np.array([0.1, 0.2, 0.8, 0.7, 0.3, 0.2])
        assert_allclose(ffe.get_cursor(sbr, 2, 3, 4), [0, 0, 0.1, 0.8, 0.3, 0, 0, 0])

    def test_zero_finite_design_residual_does_not_mean_zero_physical_isi(self):
        # The selected finite model ignores a[-2] and a[+2]. Full SBR
        # convolution includes them, so claiming all physical cursors vanish
        # from the finite matrix residual would be incorrect.
        full_channel = np.array([0.06, 0.12, 1.0, 0.30, 0.09])
        design_channel = full_channel[1:4]
        taps = ffe.zf_taps(design_channel, 1)
        finite = np.convolve(design_channel, taps)[1:4]
        assert_allclose(finite[[0, 2]], 0.0, atol=1e-14)
        physical = np.convolve(full_channel, taps)
        main = 3
        self.assertGreater(np.max(np.abs(physical[[main - 1, main + 1]])), 1e-3)


class SectionalConvolutionTests(unittest.TestCase):
    def test_random_overlap_add_matches_direct_convolution(self):
        rng = np.random.default_rng(728)
        x, h = rng.normal(size=137), rng.normal(size=43)
        expected = np.convolve(x, h)
        for section_samples in [1, 7, 32, 42, 43, 64, 137, 200]:
            with self.subTest(section_samples=section_samples):
                actual = ffe.sectional_convolve(x, h, section_samples=section_samples)
                self.assertEqual(actual.size, expected.size)
                assert_allclose(actual, expected, atol=2e-13, rtol=2e-13)

    def test_impulses_on_section_boundaries_and_no_circular_wrap(self):
        x = np.zeros(35)
        x[[0, 7, 8, 15, 16, 34]] = [0.2, 0.7, -0.4, 1.3, 0.1, 2]
        h = np.zeros(21)
        h[[0, 3, 20]] = [0.3, -0.1, 0.9]
        actual = ffe.sectional_convolve(x, h, section_samples=8)
        assert_allclose(actual, np.convolve(x, h), atol=1e-14)
        self.assertAlmostEqual(actual[-1], 1.8, places=13)
        self.assertAlmostEqual(actual[0], 0.06, places=13)

    def test_delayed_last_sample_does_not_wrap_to_start(self):
        x = np.zeros(17)
        x[-1] = 2
        h = np.zeros(12)
        h[-1] = 0.5
        actual = ffe.sectional_convolve(x, h, section_samples=6)
        assert_allclose(actual[:-1], 0, atol=1e-14)
        self.assertAlmostEqual(actual[-1], 1.0, places=13)

    def test_finite_rise_waveform_preserves_previous_level_across_sections(self):
        levels = ffe.ffe_levels([0, 1, 1, 0, 1, 0, 0, 1], [-0.1, 0.7, -0.2])
        waveform = ffe.levels_to_waveform(levels, 20, 5)
        split_bit = 4
        left = ffe.levels_to_waveform(levels[:split_bit], 20, 5)
        right = ffe.levels_to_waveform(levels[split_bit:], 20, 5, previous=levels[split_bit - 1])
        assert_allclose(np.concatenate([left, right]), waveform)
        h = np.array([0, 0.3, 0.2, 0.1])
        assert_allclose(ffe.sectional_convolve(waveform, h, section_samples=80), np.convolve(waveform, h), atol=1e-14)


class EyeHistogramTests(unittest.TestCase):
    def test_voltage_rows_and_time_columns_have_correct_orientation(self):
        waveform = np.tile([0.1, 0.8, -0.8, -0.1], 7)
        probability, extent = ffe.eye_histogram(waveform, 4, voltage_bins=4, voltage_range=(-1, 1))
        self.assertEqual(probability.shape, (4, 4))
        self.assertEqual(list(extent), [0, 1, -1, 1])
        assert_allclose(np.sum(probability, axis=0), 1.0)
        self.assertEqual(np.argmax(probability[:, 0]), 2)
        self.assertEqual(np.argmax(probability[:, 1]), 3)
        self.assertEqual(np.argmax(probability[:, 2]), 0)
        self.assertEqual(np.argmax(probability[:, 3]), 1)

    def test_histogram_is_conditional_probability_even_for_partial_ui(self):
        probability, _ = ffe.eye_histogram([0, 1, 0, 1, 0], 2, voltage_bins=2, voltage_range=(0, 1))
        assert_allclose(np.sum(probability, axis=0), 1.0)
        assert_allclose(probability, [[1.0, 0.0], [0.0, 1.0]])

    def test_eye_phase_shift_corresponds_to_column_rotation(self):
        waveform = np.tile([0.1, 0.8, -0.8, -0.1], 9)
        probability, _ = ffe.eye_histogram(waveform, 4, voltage_bins=10, voltage_range=(-1, 1))
        shifted, _ = ffe.eye_histogram(np.roll(waveform, 1), 4, voltage_bins=10, voltage_range=(-1, 1))
        assert_allclose(shifted, np.roll(probability, 1, axis=1))


class FullSimulationTests(unittest.TestCase):
    """Small end-to-end run checks deliverables and shared-input comparison."""

    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory(dir=ROOT / "tests")
        cls.output = Path(cls.folder.name)
        cls.inputs = []
        apply_ffe = ffe.ffe_levels

        def record_input(bits, taps, vpp=2):
            cls.inputs.append((np.array(bits, copy=True), np.array(taps, copy=True), vpp))
            return apply_ffe(bits, taps, vpp)

        with mock.patch.object(ffe, "ffe_levels", side_effect=record_input):
            with contextlib.redirect_stdout(io.StringIO()):
                cls.summary = ffe.run_simulation(DATA, cls.output, bits=128, seed=19)

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def test_all_required_cases_and_four_case4_snapshots_exist(self):
        expected_cases = [(0, 0), (1, 0), (0, 1), (1, 1), (1, 2)]
        cases = self.summary["cases"]
        self.assertEqual([(c["npre"], c["npost"]) for c in cases], expected_cases)
        for case in cases:
            path = self.output / case["eye_file"]
            self.assertGreater(path.stat().st_size, 1000)
            self.assertAlmostEqual(case["tap_l1_norm"], 1.0, places=13)
            self.assertLessEqual(case["tx_swing_v"], 2.0 + 1e-12)
        for name in ["case4_tx_snapshot1.png", "case4_tx_snapshot2.png",
                     "case4_rx_snapshot1.png", "case4_rx_snapshot2.png",
                     "metrics.csv", "summary.json",
                     "eyes_overview.png", "ffe_taps.png", "equalized_sbr.png"]:
            with self.subTest(file=name):
                self.assertGreater((self.output / name).stat().st_size, 100)
        for signal in ["tx", "rx"]:
            for snapshot in [1, 2]:
                values = np.loadtxt(self.output / f"case4_{signal}_snapshot{snapshot}.csv",
                                    delimiter=",", skiprows=1)
                self.assertEqual(values.shape[1], 2)
                self.assertGreater(values.shape[0], 100)
                assert_allclose(np.diff(values[:, 0]), 0.5, atol=1e-12)
        snapshots = self.summary["cases"][3]["snapshots"]
        self.assertEqual(len(snapshots), 4)
        self.assertEqual(snapshots[0]["source_payload_start_bit"], snapshots[1]["source_payload_start_bit"])
        self.assertEqual(snapshots[2]["source_payload_start_bit"], snapshots[3]["source_payload_start_bit"])
        self.assertNotEqual(snapshots[0]["source_payload_start_bit"], snapshots[2]["source_payload_start_bit"])

    def test_all_cases_use_identical_unipolar_source_sequence(self):
        self.assertEqual(len(self.inputs), 5)
        for bits, taps, vpp in self.inputs:
            assert_allclose(bits, self.inputs[0][0], atol=0)
            self.assertTrue(np.all(np.isin(bits, [0, 1])))
            self.assertEqual(vpp, 2.0)
            self.assertAlmostEqual(np.sum(np.abs(taps)), 1.0, places=13)
        self.assertEqual(self.summary["parameters"]["source_levels_v"], [0, 2])
        self.assertAlmostEqual(self.summary["channel"]["dc_gain"], 0.484733287, places=10)

    def test_full_channel_metrics_match_independent_reference(self):
        self.assertEqual(self.summary["channel"]["main_sbr_index"], 5804)
        expected_taps = [
            [1.0], [-0.09958441, 0.90041559], [0.68714306, -0.31285694],
            [-0.07062923, 0.63861067, -0.29076009],
            [-0.07102079, 0.64215105, -0.27609499, -0.01073317],
        ]
        for case, taps in zip(self.summary["cases"], expected_taps):
            assert_allclose(case["taps"], taps, atol=1e-8, rtol=0)
            self.assertLess(case["zf_matrix_residual_v"], 1e-13)
        case4 = self.summary["cases"][3]
        assert_allclose(case4["physical_selected_cursors_at_design_phase_v"],
                        [-0.00081067, 0.22466316, -0.00587062], atol=1e-8, rtol=0)
        self.assertGreater(case4["physical_selected_residual_v"], 0.005)
        # Saving/reloading the summary must retain actual finite JSON values.
        saved = json.loads((self.output / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(saved["parameters"]["payload_bits"], 128)
        self.assertNotIn("NaN", (self.output / "summary.json").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
