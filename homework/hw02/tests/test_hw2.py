"""Independent numerical checks for the HW2 reconstruction."""

import sys
from pathlib import Path
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hw2 import (ROOT, TARGET_RISE, alpha_for_rise, circuit_ramp_reflection, circuit_kaiser_reflection,
                 kaiser_impulse, kaiser_spectrum, kaiser_step, normalized_rise,
                 read_touchstone, sparameters_from_circuit, windowed_reflection)


class KaiserTests(unittest.TestCase):
    def test_analytic_impulse_matches_independent_fourier_quadrature(self):
        frequency = np.linspace(-0.5, 0.5, 40001)
        for alpha in (0, 1, 4):
            for q in (-4, -1.5, 0, 1.5, 4):
                integral = np.trapezoid(kaiser_spectrum(frequency, 1, alpha)
                                        * np.exp(2j * np.pi * frequency * q), frequency)
                self.assertAlmostEqual(integral.real, float(kaiser_impulse(np.array(q), 1, alpha)), delta=2e-9)
                self.assertAlmostEqual(integral.imag, 0, delta=1e-14)

    def test_step_derivative_and_symmetry(self):
        q = np.linspace(-4, 4, 81)
        for alpha in (0, 1, 4):
            np.testing.assert_allclose(kaiser_step(q, alpha) + kaiser_step(-q, alpha), 1, atol=1e-13)
            derivative = (kaiser_step(q + 1e-5, alpha) - kaiser_step(q - 1e-5, alpha)) / 2e-5
            np.testing.assert_allclose(derivative, kaiser_impulse(q, 1, alpha), atol=1e-9)

    def test_rise_values_and_40ps_solution(self):
        expected = [0.8916388524, 1.4207550337, 2.0193203531, 2.4834136373, 2.8736927099]
        np.testing.assert_allclose([normalized_rise(a) for a in range(5)], expected, atol=1e-10)
        alpha = alpha_for_rise(40e9, TARGET_RISE)
        self.assertAlmostEqual(alpha, 1.26911023801, places=10)
        self.assertAlmostEqual(normalized_rise(alpha) / 40e9, TARGET_RISE, delta=1e-22)

    def test_window_has_compact_support(self):
        np.testing.assert_array_equal(kaiser_spectrum(np.array([-1, 0, 1]), 1, 2), [0, 1, 0])


class TouchstoneTests(unittest.TestCase):
    def parse(self, content):
        with tempfile.TemporaryDirectory(dir=ROOT) as directory:
            path = Path(directory) / "example.s2p"
            path.write_text(content, encoding="utf-8")
            return read_touchstone(path)

    def test_ri_two_port_order_wrapping_and_comment(self):
        data = self.parse("# MHZ S RI R 75\n1 1 2 3 4 ! comment\n5 6 7 8\n")
        np.testing.assert_array_equal(data.frequency_hz, [1e6])
        np.testing.assert_array_equal(data.s[0], [[1 + 2j, 5 + 6j], [3 + 4j, 7 + 8j]])
        self.assertEqual(data.reference_ohm, 75)

    def test_ma_and_db_phase_conversion(self):
        ma = self.parse("# HZ S MA R 50\n0 1 90 2 180 3 -90 4 0\n")
        db = self.parse("# HZ S DB R 50\n0 0 90 6.020599913279624 180 9.542425094393248 -90 12.041199826559248 0\n")
        np.testing.assert_allclose(ma.s, db.s, atol=1e-14)
        np.testing.assert_allclose(ma.s[0], [[1j, -3j], [-2, 4]], atol=1e-14)

    def test_rejects_incomplete_records(self):
        with self.assertRaises(ValueError):
            self.parse("# GHZ S RI R 50\n0 1 0\n")


class ReconstructionTests(unittest.TestCase):
    def test_full_hermitian_ifft_matches_direct_frequency_sum_and_dc_once(self):
        frequency = np.array([0., 1., 2.])
        response = np.array([.4, .2 + .1j, -.3 + .4j])
        time, impulse, step, diagnostic = windowed_reflection(frequency, response, 0, delay_s=0, n_fft=64)
        expected = response[0].real + 2 * np.real(response[1] * np.exp(2j * np.pi * time))
        expected += np.real(response[2] * np.exp(4j * np.pi * time))
        np.testing.assert_allclose(impulse, expected, atol=2e-14)
        self.assertAlmostEqual(diagnostic["integrated_impulse"], .4, places=14)
        self.assertLess(diagnostic["max_ifft_imaginary_relative"], 1e-14)

    def test_rejects_invalid_frequency_grid_or_dc(self):
        for frequency, response in [(np.array([0, 1, 3]), np.ones(3)),
                                    (np.array([1, 2, 3]), np.ones(3)),
                                    (np.array([0, 1, 2]), np.array([1j, 1, 1]))]:
            with self.assertRaises(ValueError):
                windowed_reflection(frequency, response, 1, n_fft=64)

    def test_supplied_s11_matches_independent_circuit_and_is_passive(self):
        data = read_touchstone(ROOT / "data" / "chebyshev_lpf.s2p")
        circuit = sparameters_from_circuit(data.frequency_hz)
        self.assertLess(np.max(np.abs(circuit[:, 0, 0] - data.s[:, 0, 0])), 1e-6)
        power = np.abs(data.s[:, 0, 0]) ** 2 + np.abs(data.s[:, 1, 0]) ** 2
        np.testing.assert_allclose(power, 1, atol=2e-6)
        self.assertEqual(data.frequency_hz.size, 2001)
        self.assertEqual(data.frequency_hz[-1], 20e9)

    def test_ngspice_matches_independent_state_space_response(self):
        reference = np.loadtxt(ROOT / "simulation" / "spice_reference.csv", delimiter=",", skiprows=1)
        time, source, port, reflected = reference.T
        np.testing.assert_allclose(port - source / 2, reflected, atol=2e-12)
        state = circuit_ramp_reflection(time)
        self.assertLess(np.max(np.abs(reflected - state)), 1e-4)
        # A matched source sees no reflection before the line round trip.
        self.assertLess(np.max(np.abs(state[time < 500e-12])), 1e-12)

    def test_same_kaiser_input_agrees_between_circuit_ode_and_frequency_ifft(self):
        data = read_touchstone(ROOT / "data" / "chebyshev_lpf.s2p")
        F = 2 * data.frequency_hz[-1]
        alpha = alpha_for_rise(F, TARGET_RISE)
        time, _, step, _ = windowed_reflection(data.frequency_hz, data.s[:, 0, 0], alpha)
        comparison = np.linspace(0, 3e-9, 6001)
        ode = circuit_kaiser_reflection(comparison, F, alpha)
        idft = np.interp(comparison, time, step)
        self.assertLess(np.max(np.abs(ode - idft)), 1e-4)


if __name__ == "__main__":
    unittest.main()
