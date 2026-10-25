"""Independent physical and numerical checks for HW12.

Run with ``python -m unittest discover -s homework/hw12/tests -v``.
The oracles use loaded-network admittance, analytic ramps, direct convolution,
and integrated theoretical distributions, rather than reproducing the streaming
simulation. All fixtures are anonymous and generated locally.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.special import ndtr


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hw12_under_test", ROOT / "src" / "hw12.py")
hw12 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = hw12
SPEC.loader.exec_module(hw12)


def reference_ramps(bits, jitter, dt, ui, rise, start, stop, bit_start=0, vpp=1.0):
    """Superpose the exact response of each level change to a linear ramp."""
    times = np.arange(start, stop, dtype=float) * dt
    result = np.zeros(times.size)
    previous = 0.0
    for k, bit in enumerate(bits):
        level = float(bit) * vpp
        onset = (bit_start + k) * ui + jitter[k]
        result += (level - previous) * np.clip((times - onset) / rise, 0.0, 1.0)
        previous = level
    return result


def eye_from_row(row):
    """Give each time phase the same exposure while isolating one voltage bin."""
    row = np.asarray(row, dtype=float)
    exposure = max(1.0, float(np.max(row)) * 2)
    return np.vstack((exposure - row, row, np.zeros_like(row)))


class TouchstoneAndCircuitTests(unittest.TestCase):
    def read_text(self, text):
        with tempfile.TemporaryDirectory(dir=ROOT / "tests") as folder:
            path = Path(folder) / "anonymous_fixture.s2p"
            path.write_text(text, encoding="ascii")
            return hw12.read_touchstone_2port(path)

    def test_ma_degrees_units_and_nonreciprocal_orientation(self):
        frequency, scattering, z0 = self.read_text(
            "! S2P record order is S11,S21,S12,S22\n# GHz S MA R 75\n"
            "0 0.1 0 0.8 90 0.2 -90 0.3 180\n"
            "0.25 0.2 0 0.7 45 0.4 -45 0.5 0 ! inline comment\n"
        )
        assert_allclose(frequency, [0, 250e6])
        self.assertEqual(z0, 75)
        self.assertEqual(scattering.shape, (2, 2, 2))
        assert_allclose(scattering[0], [[0.1, -0.2j], [0.8j, -0.3]], atol=1e-15)
        assert_allclose(scattering[1, 1, 0], 0.7 * np.exp(1j * np.pi / 4))

    def test_ri_wrapped_record(self):
        frequency, scattering, _ = self.read_text(
            "# MHz S RI R 50\n"
            "0 0.1 0 0.5 -0.25\n 0.125 0.2 0.25 -0.5\n"
            "10 0.2 0 0.7 -0.1 0.3 0.4 0.6 -0.2\n"
        )
        assert_allclose(frequency, [0, 10e6])
        assert_allclose(scattering[0], [[0.1, 0.125 + 0.2j], [0.5 - 0.25j, 0.25 - 0.5j]])

    def test_db_decibels_and_khz(self):
        frequency, scattering, _ = self.read_text(
            "# kHz S DB R 50\n"
            "0 -20 0 -6.020599913279624 180 -40 90 -20 0\n"
            "100 -20 0 -6.020599913279624 180 -40 90 -20 0\n"
        )
        assert_allclose(frequency, [0, 100e3])
        assert_allclose(scattering[0], [[0.1, 0.01j], [-0.5, 0.1]], atol=1e-14)

    def test_loaded_circuit_agrees_with_independent_admittance_solution(self):
        frequency = np.array([0, 100e6, 3e9, 20e9])
        scattering = np.array([
            [[0.07, 0.11], [0.76, -0.04]],
            [[0.08 - 0.02j, 0.13 + 0.01j], [0.73 - 0.12j, -0.03 + 0.01j]],
            [[0.02 + 0.04j, 0.06 - 0.02j], [0.42 + 0.18j, -0.08 + 0.04j]],
            [[0.13 - 0.01j, 0.03 + 0.05j], [-0.1 + 0.25j, 0.06 - 0.03j]],
        ], dtype=complex)
        z0, rs, rl, cs, cl = 50.0, 20.0, 85.0, 0.3e-12, 0.3e-12
        expected = []
        identity = np.eye(2)
        for f, s in zip(frequency, scattering):
            # V=(I+S)a, Iport=(I-S)a/Z0. The circuit is solved independently
            # as two nodal voltages with Norton source 1/Rs at node 1.
            admittance = (identity - s) @ np.linalg.inv(identity + s) / z0
            loading = np.diag([1 / rs + 2j * np.pi * f * cs,
                               1 / rl + 2j * np.pi * f * cl])
            expected.append(np.linalg.solve(admittance + loading, [1 / rs, 0])[1])
        result = hw12.channel_transfer(frequency, scattering, z0, rs=rs, cs=cs, rl=rl, cl=cl)
        assert_allclose(result, expected, rtol=2e-13, atol=1e-14)
        reversed_result = hw12.channel_transfer(frequency, scattering.transpose(0, 2, 1), z0)
        self.assertGreater(np.max(np.abs(result - reversed_result)), 0.1)

    def test_zero_length_through_network_has_resistor_divider_dc_gain(self):
        frequency = np.array([0, 1e9, 10e9])
        scattering = np.tile([[0.0, 1.0], [1.0, 0.0]], (3, 1, 1)).astype(complex)
        actual = hw12.channel_transfer(frequency, scattering, 50.0)
        expected = 1 / (1 + 20 / 85 + 2j * np.pi * frequency * 20 * 0.6e-12)
        assert_allclose(actual, expected, rtol=1e-12, atol=1e-13)


class ImpulseAndStreamingTests(unittest.TestCase):
    def test_complete_band_integer_delay_is_discrete_unit_impulse(self):
        dt = 1e-12
        frequency = np.arange(501) * 1e9
        transfer = np.exp(-2j * np.pi * frequency * (23 * dt))
        actual = hw12.build_fir(frequency, transfer, dt, truncation=60 * dt, kaiser_a=0)
        expected = np.zeros_like(actual)
        expected[23] = 1.0
        assert_allclose(actual, expected, atol=2e-13)

    def test_streaming_matches_full_direct_convolution_for_irregular_blocks(self):
        rng = np.random.default_rng(1212)
        signal = rng.normal(size=263)
        kernel = rng.normal(size=31) * np.exp(-np.arange(31) / 7)
        boundaries = [0, 1, 4, 41, 42, 105, 263]
        blocks = [signal[a:b] for a, b in zip(boundaries[:-1], boundaries[1:])]
        output = list(hw12.stream_fir(iter(blocks), kernel))
        self.assertEqual([len(block) for block in output], [len(block) for block in blocks])
        expected = np.convolve(signal, kernel, mode="full")[:len(signal)]
        assert_allclose(np.concatenate(output), expected, atol=2e-12, rtol=1e-12)

    def test_streaming_delta_kernel_and_single_sample_chunks(self):
        signal = np.array([1, -2, 0, 5, 8], dtype=float)
        for kernel in [np.array([1.0]), np.array([0, 0, 0.25])]:
            with self.subTest(kernel=kernel.tolist()):
                output = np.concatenate(list(hw12.stream_fir((signal[i:i + 1] for i in range(len(signal))), kernel)))
                assert_allclose(output, np.convolve(signal, kernel)[:len(signal)], atol=1e-13)


class WaveformBoundaryTests(unittest.TestCase):
    def test_negative_first_jitter_evaluates_existing_ramp_at_time_zero(self):
        bits = np.array([1, 1, 0], dtype=np.int8)
        jitter = np.array([-2.0, 0, 0]) * 1e-12
        actual = hw12.generate_waveform_segment(bits, jitter, 1e-12, 10e-12, 4e-12, 0, 8)
        assert_allclose(actual, [0.5, 0.75, 1, 1, 1, 1, 1, 1], atol=1e-14)

    def test_jitter_sign_and_repeated_bits(self):
        bits = np.array([0, 1, 1, 0, 1, 0, 0], dtype=np.int8)
        jitter = np.array([0, 2, -2, -1, 1, -2, 2], dtype=float) * 1e-12
        actual = hw12.generate_waveform_segment(bits, jitter, 0.5e-12, 10e-12, 4e-12, 0, 140)
        expected = reference_ramps(bits, jitter, 0.5e-12, 10e-12, 4e-12, 0, 140)
        assert_allclose(actual, expected, atol=1e-13)
        # Positive jitter delays the rise from 10 ps to 12 ps, whose midpoint is 14 ps.
        self.assertAlmostEqual(actual[28], 0.5, places=12)
        # A jitter value attached to a repeated bit cannot introduce an extra edge.
        self.assertAlmostEqual(actual[39], 1.0, places=12)

    def test_arbitrary_segment_boundaries_equal_one_contiguous_generation(self):
        rng = np.random.default_rng(1203)
        bits = rng.integers(0, 2, size=50, dtype=np.int8)
        jitter = rng.uniform(-2e-12, 2e-12, size=len(bits))
        dt, ui, rise = 0.5e-12, 25e-12, 5e-12
        boundaries = [0, 3, 47, 50, 99, 201, 503, 1001, 2400]
        actual = np.concatenate([hw12.generate_waveform_segment(bits, jitter, dt, ui, rise, a, b)
                                 for a, b in zip(boundaries[:-1], boundaries[1:])])
        expected = reference_ramps(bits, jitter, dt, ui, rise, 0, boundaries[-1])
        assert_allclose(actual, expected, atol=5e-13)
        self.assertEqual(len(actual), boundaries[-1])

    def test_negative_guard_bit_origin_preserves_absolute_time_and_exact_count(self):
        rng = np.random.default_rng(1204)
        bits = rng.integers(0, 2, size=29, dtype=np.int8)
        jitter = rng.uniform(-1e-12, 1e-12, size=29)
        dt, ui, rise, origin = 1e-12, 25e-12, 5e-12, -8
        sample_count = 17 * 25
        actual = hw12.generate_waveform_segment(bits, jitter, dt, ui, rise, 0, sample_count, bit_start=origin, vpp=0.8)
        expected = reference_ramps(bits, jitter, dt, ui, rise, 0, sample_count, bit_start=origin, vpp=0.8)
        assert_allclose(actual, expected, atol=5e-13)
        self.assertEqual(len(actual), sample_count)

    def test_overlapping_true_transition_ramps_are_rejected(self):
        bits = np.array([1, 0, 1], dtype=np.int8)
        jitter = np.array([0, -8, 0], dtype=float) * 1e-12
        with self.assertRaises(ValueError):
            hw12.generate_waveform_segment(bits, jitter, 1e-12, 10e-12, 4e-12, 0, 30)

    def test_repeated_bit_jitter_does_not_create_false_overlap(self):
        bits = np.array([1, 1, 0], dtype=np.int8)
        jitter = np.array([0, -9, 0], dtype=float) * 1e-12
        actual = hw12.generate_waveform_segment(bits, jitter, 1e-12, 10e-12, 4e-12, 0, 30)
        expected = reference_ramps(bits, jitter, 1e-12, 10e-12, 4e-12, 0, 30)
        assert_allclose(actual, expected, atol=1e-13)


class EyeAndJitterStatisticsTests(unittest.TestCase):
    def test_eye_accumulation_exact_counts_phase_and_no_pseudocount(self):
        nphase, start = 10, 17
        waveform = np.tile([0, 0.5, 1.0], 35)[:103]
        actual = np.zeros((3, nphase), dtype=np.int64)
        expected = np.zeros_like(actual)
        for i, voltage in enumerate(waveform):
            expected[int(np.floor((voltage + 0.25) / 0.5)), (start + i) % nphase] += 1
        for a, b in [(0, 7), (7, 19), (19, 81), (81, 103)]:
            hw12.accumulate_eye(actual, waveform[a:b], start + a, nphase, -0.25, 0.5)
        np.testing.assert_array_equal(actual, expected)
        self.assertEqual(actual.sum(), len(waveform))
        np.testing.assert_array_equal(actual.sum(axis=0), np.bincount((np.arange(len(waveform)) + start) % nphase, minlength=nphase))

    def test_one_ui_wrap_keeps_boundary_probability_and_endpoint_excluded_axis(self):
        row = np.zeros(8)
        row[[6, 7, 0]] = [25, 50, 25]
        dt = 0.125e-12
        result = hw12.analyze_eye_slice(eye_from_row(row), 0, 1, 1.5, dt, reference_index=7)
        expected = np.zeros(8)
        expected[[3, 4, 5]] = [0.25, 0.5, 0.25]
        assert_allclose(result["pmf"], expected, atol=1e-15)
        assert_allclose(result["time_seconds"], np.arange(-4, 4) * dt, atol=1e-27)
        self.assertAlmostEqual(float(np.sum(result["pdf"]) * dt), 1.0, places=13)
        self.assertAlmostEqual(result["mean_seconds"] / dt, 0.0, places=13)
        self.assertAlmostEqual(result["variance_seconds2"] / dt**2, 0.5, places=13)
        self.assertEqual(np.count_nonzero(result["pmf"]), 3)

    def test_moments_use_centered_variance_and_preserve_nonzero_mean(self):
        row = np.zeros(10)
        row[[1, 3, 6]] = [1, 2, 1]
        dt = 0.1e-12
        result = hw12.analyze_eye_slice(eye_from_row(row), 0, 1, 1.5, dt, reference_index=3)
        self.assertAlmostEqual(result["mean_seconds"] / dt, 0.25, places=13)
        self.assertAlmostEqual(result["variance_seconds2"] / dt**2, 3.1875, places=12)
        self.assertAlmostEqual(result["rms_seconds"] / dt, np.sqrt(3.1875), places=12)

    def test_zero_jitter_is_exact_delta_not_smoothed(self):
        row = np.zeros(200)
        row[199] = 120
        result = hw12.analyze_eye_slice(eye_from_row(row), 0, 1, 1.5, 0.125e-12)
        self.assertEqual(result["reference_index"], 199)
        self.assertEqual(np.count_nonzero(result["pmf"]), 1)
        self.assertEqual(result["pmf"][100], 1)
        self.assertEqual(result["variance_seconds2"], 0)
        self.assertEqual(result["rms_seconds"], 0)

    def test_empty_threshold_slice_raises_instead_of_fabricating_distribution(self):
        eye = np.zeros((3, 20), dtype=np.int64)
        eye[0, :] = 100
        with self.assertRaises(ValueError):
            hw12.analyze_eye_slice(eye, 0, 1, 1.5, 1e-12)

    def test_analytic_sj_rj_and_independent_combination_variance(self):
        # Integrate continuous probabilities over each exact time bin. The
        # expected values are A^2/2, sigma^2, and their sum, in ps^2.
        nphase, dt_ps = 1000, 0.025
        centers = (np.arange(nphase) - nphase // 2) * dt_ps
        edges = np.r_[centers - dt_ps / 2, centers[-1] + dt_ps / 2]
        amplitude, sigma = 3.0, 1.0
        sj_cdf = 0.5 + np.arcsin(np.clip(edges / amplitude, -1, 1)) / np.pi
        rj_cdf = ndtr(edges / sigma)
        phase = (np.arange(2048) + 0.5) * (2 * np.pi / 2048)
        combined_cdf = np.mean(ndtr((edges[:, None] - amplitude * np.sin(phase)) / sigma), axis=1)
        for name, cdf, expected_variance in [
            ("SJ", sj_cdf, 4.5), ("RJ", rj_cdf, 1.0), ("SJ+RJ", combined_cdf, 5.5)
        ]:
            with self.subTest(case=name):
                row = np.diff(cdf)
                result = hw12.analyze_eye_slice(eye_from_row(row), 0, 1, 1.5, dt_ps * 1e-12, reference_index=nphase // 2)
                self.assertAlmostEqual(float(np.sum(result["pmf"])), 1, places=13)
                self.assertAlmostEqual(result["mean_seconds"] * 1e12, 0, delta=1e-10)
                self.assertAlmostEqual(result["variance_seconds2"] * 1e24, expected_variance, delta=0.002)


if __name__ == "__main__":
    unittest.main(verbosity=2)
