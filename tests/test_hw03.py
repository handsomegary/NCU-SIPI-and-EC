"""Numerical and boundary invariants for the HW3 rewrite (NumPy + unittest)."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


SOURCE = Path(__file__).resolve().parents[1] / "homework" / "hw03" / "analyze.py"
SPEC = importlib.util.spec_from_file_location("hw03", SOURCE)
hw = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hw)


class HW3Tests(unittest.TestCase):
    def test_sectional_matches_direct_including_partial_block(self):
        rng = np.random.default_rng(23)
        bits = rng.integers(0, 2, 103).astype(bool)
        h = rng.normal(size=47) * 0.02
        expected = np.convolve(hw.input_section(bits, 10, 2).astype(float), h)
        for block_bits in (1, 4, 16, 100, 200):
            actual, info = hw.simulate(bits, h, 10, 2, block_bits)
            np.testing.assert_allclose(actual, expected, atol=2e-7, rtol=2e-6)
            self.assertEqual(info["processed_bit_count"], 103)
            self.assertEqual(info["final_block_bits"], 103 - (info["block_count"] - 1) * block_bits)
            self.assertGreaterEqual(info["fft_points"], block_bits * 10 + len(h) - 1)

    def test_boundary_preserves_true_previous_bit(self):
        for values in ([1, 1, 1, 1], [0, 0, 1, 1], [1, 1, 0, 0], [0, 1, 0, 1]):
            bits = np.array(values, dtype=bool)
            whole = hw.input_section(bits, 10, 2)
            joined = np.concatenate([hw.input_section(bits[:2], 10, 2),
                                     hw.input_section(bits[2:], 10, 2, bool(bits[1]))])
            np.testing.assert_array_equal(whole, joined)
        np.testing.assert_allclose(hw.input_section(np.array([1, 0]), 5, 2),
                                   [0, .5, 1, 1, 1, 1, .5, 0, 0, 0])

    def test_impulse_scaling_and_delay(self):
        dt = 1e-12
        # Constant real transfer on all but Nyquist bin: exact known sampled FIR.
        f = np.arange(0, 500e9, 1e9)
        delay = 20 * dt
        transfer = 0.75 * np.exp(-2j * np.pi * f * delay)
        h, periodic, info = hw.sampled_impulse(f, transfer, dt, 0.5e-9, alpha=0)
        self.assertEqual(np.argmax(np.abs(periodic)), 20)
        self.assertAlmostEqual(periodic.sum(), .75, places=12)
        self.assertAlmostEqual(periodic[20], .75 * 999 / 1000, places=12)
        self.assertEqual(info["impulse_period_points"], 1000)
        self.assertAlmostEqual(info["peak_delay_s"], delay, places=20)
        # Independent lecture half-DC twice-real IFFT agrees with irfft scaling.
        spectrum = np.zeros(1000, complex)
        spectrum[:len(f)] = transfer
        spectrum[0] *= .5
        np.testing.assert_allclose(periodic, 2 * np.fft.ifft(spectrum).real, atol=1e-14)

    def test_integer_grid_and_rejections(self):
        self.assertEqual(hw.integer_points(20e-12, 2e-12, "rise"), 10)
        with self.assertRaises(ValueError):
            hw.integer_points(20e-12, 3e-12, "rise")
        with self.assertRaises(ValueError):
            hw.input_section(np.array([1]), 10, 11)
        with self.assertRaises(ValueError):
            hw.sampled_impulse(np.array([0, 1e9]), np.ones(2), 1e-9, .1e-9)

    def test_memmap_matches_memory_and_ads_interpolation(self):
        bits = np.array([1, 1, 0, 1, 1], dtype=bool)
        h = np.array([.1, .7, .2])
        expected, _ = hw.simulate(bits, h, 10, 2, block_bits=2)
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            actual, _ = hw.simulate(bits, h, 10, 2, 2, base / "waveform")
            np.testing.assert_array_equal(actual, expected)
            times = np.array([0, .3, 1.2, 3.8])
            reference = np.interp(times, np.arange(len(expected)) * .1, expected)
            ads = base / "ads.txt"
            np.savetxt(ads, np.column_stack([times, reference]), header="time vout", comments="")
            stats, windows = hw.compare_ads(ads, actual, .1, [(0, 4)], chunk_rows=2)
            self.assertLess(stats["max_absolute_error_v"], 1e-14)
            self.assertEqual(stats["sample_count"], 4)
            self.assertEqual(windows[0]["sample_count"], 4)
            cropped, cropped_windows = hw.compare_ads(ads, actual, .1, [(0, 1)],
                                                      chunk_rows=2, valid_end_s=1)
            self.assertEqual(cropped["sample_count"], 2)
            self.assertEqual(cropped["excluded_sample_count"], 2)
            self.assertEqual(cropped["file_sample_count"], 4)
            self.assertAlmostEqual(cropped["file_end_s"], 3.8)
            self.assertEqual(cropped_windows[0]["sample_count"], 2)
            del actual

    def test_touchstone_and_bits_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            channel = base / "channel.s2p"
            channel.write_text("# GHz S MA R 50\n0 0 0 .5 0 .5 0 0 0\n1 0 0 .5 -90 .5 -90 0 0\n")
            f, s = hw.read_s21(channel)
            np.testing.assert_allclose(f, [0, 1e9])
            np.testing.assert_allclose(s, [.5, -.5j], atol=1e-15)
            pattern = base / "bits.txt"; pattern.write_text("01\n10 ")
            np.testing.assert_array_equal(hw.read_bits(pattern), [0, 1, 1, 0])
            pattern.write_text("01x")
            with self.assertRaises(ValueError):
                hw.read_bits(pattern)


if __name__ == "__main__":
    unittest.main()
