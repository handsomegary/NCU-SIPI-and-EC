"""Behavioral checks for portable HW1 analysis (no assignment answers required)."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "homework" / "hw01" / "analyze.py"
SPEC = importlib.util.spec_from_file_location("hw01_analyze", SCRIPT)
hw01 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hw01)


class SpectrumTests(unittest.TestCase):
    def test_one_sided_dc_and_nyquist_have_physical_amplitudes(self):
        n = np.arange(256)
        signal = 0.4 + 1.7 * np.cos(2 * np.pi * 17 * n / n.size) + 0.6 * (-1.0) ** n
        freq, amplitude, gain = hw01.amplitude_spectrum(signal, 1 / 1024, 0)
        self.assertAlmostEqual(amplitude[0], 0.4, places=12)
        self.assertAlmostEqual(amplitude[17], 1.7, places=12)
        self.assertAlmostEqual(amplitude[-1], 0.6, places=12)
        self.assertAlmostEqual(freq[-1], 512, places=12)
        self.assertEqual(gain, 1)

    def test_coherent_gain_correction_preserves_an_isolated_tone(self):
        n = np.arange(1024)
        signal = 2.3 * np.cos(2 * np.pi * 173 * n / n.size + 0.3)
        for alpha in (1, 4, 5):
            with self.subTest(alpha=alpha):
                _, amplitude, gain = hw01.amplitude_spectrum(signal, 1 / 4096, alpha)
                self.assertLess(abs(amplitude[173] - 2.3), 2e-5)
                self.assertLess(gain, 1)

    def test_odd_fft_last_interior_bin_is_doubled(self):
        n = np.arange(255)
        signal = 0.7 * np.cos(2 * np.pi * 127 * n / n.size)
        _, amplitude, _ = hw01.amplitude_spectrum(signal, 1, 0)
        self.assertAlmostEqual(amplitude[-1], 0.7, places=12)

    def test_sampling_errors_are_rejected(self):
        t = np.arange(32) / 1000
        x = np.sin(2 * np.pi * 53 * t)
        for bad_t in (t[::-1], np.r_[t[:17], t[17:] + 0.0001]):
            with self.subTest(bad_t=bad_t):
                with self.assertRaises(ValueError):
                    hw01.validate_samples(bad_t, x)
        x[4] = np.nan
        with self.assertRaises(ValueError):
            hw01.validate_samples(t, x)


class EstimatorTests(unittest.TestCase):
    def test_recovers_unseen_frequencies_amplitudes_phases_and_weak_tone(self):
        # Deliberately different units/frequencies/phases from the assignment.
        t = np.arange(768) / 4096
        frequencies = np.array([183.4, 207.8, 633.25])
        amplitudes = np.array([1.2, 0.8, 0.006])
        phases = np.array([0.23, -0.71, 1.1])
        signal = 0.15 + np.sin(2 * np.pi * t[:, None] * frequencies + phases) @ amplitudes
        result, reconstructed = hw01.fit_tones(t, signal, max_components=5)
        self.assertEqual(result["q"], 3)
        np.testing.assert_allclose([c["frequency_hz"] for c in result["components"]], frequencies, atol=1e-7, rtol=0)
        np.testing.assert_allclose([c["amplitude_v"] for c in result["components"]], amplitudes, atol=1e-8, rtol=0)
        np.testing.assert_allclose([c["phase_rad"] for c in result["components"]], phases, atol=1e-7, rtol=0)
        self.assertAlmostEqual(result["offset_v"], 0.15, places=9)
        self.assertLess(np.sqrt(np.mean((signal - reconstructed) ** 2)), 1e-9)

    def test_bic_avoids_adding_tones_to_seeded_noise(self):
        rng = np.random.default_rng(23)
        t = np.arange(512) / 2048
        signal = 1.1 * np.sin(2 * np.pi * 137.3 * t + 0.2) + rng.normal(0, 0.025, t.size)
        result, _ = hw01.fit_tones(t, signal, max_components=4)
        self.assertEqual(result["q"], 1)
        self.assertLess(abs(result["components"][0]["frequency_hz"] - 137.3), 0.03)
        self.assertGreater(len(result["order_selection"]), 2)


if __name__ == "__main__":
    unittest.main()
