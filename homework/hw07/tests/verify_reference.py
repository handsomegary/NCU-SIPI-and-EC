"""Independent numerical checks for the HW7 implementation.

Run: python homework/hw07/tests/verify_reference.py
References use direct convolution, a full Hermitian IFFT, time-knot source
interpolation, and the analytic continuous trapezoidal-pulse spectrum.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

BASE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("hw7_checked", BASE / "src" / "hw7.py")
hw7 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = hw7
SPEC.loader.exec_module(hw7)

PARAMETERS = dict(fp1_ghz=20.6, fp2_ghz=22.8, adc_db=-4.91, aac_db=11.48)
CHECKS = []


def check(name, condition, **measurements):
    CHECKS.append({"name": name, "passed": bool(condition), **measurements})
    if not condition:
        raise AssertionError(name)


def source_reference(bits, samples_per_ui=50):
    """Piecewise-linear interpolation through transition endpoints, in volts."""
    bits = np.asarray(bits)
    times, voltages = [], []
    previous = 0
    for index, bit in enumerate(bits):
        times.extend([index * samples_per_ui, index * samples_per_ui + samples_per_ui // 5])
        voltages.extend([2 * previous, 2 * bit])
        previous = bit
    return np.interp(np.arange(len(bits) * samples_per_ui), times, voltages)


def independent_peak(parameters):
    p1, p2 = parameters["fp1_ghz"], parameters["fp2_ghz"]
    dc, ac = 10 ** (parameters["adc_db"] / 20), 10 ** (parameters["aac_db"] / 20)
    z2 = (dc / ac * p1) ** 2
    rising = (p1 * p2) ** 2 - (p1 * p1 + p2 * p2) * z2
    squared_frequency = rising / (np.sqrt(z2 * z2 + rising) + z2) if rising > 0 else 0.0
    frequency = np.sqrt(squared_frequency)
    magnitude_squared = ac ** 2 * p2 ** 2 * (squared_frequency + z2) / ((squared_frequency + p1 ** 2) * (squared_frequency + p2 ** 2))
    return float(frequency * 1e9), float(10 * np.log10(magnitude_squared))


def independent_ctle(frequency_hz, parameters):
    """Use angular frequencies rather than the implementation's GHz form."""
    s = 2j * np.pi * frequency_hz
    p1, p2 = 2 * np.pi * parameters["fp1_ghz"] * 1e9, 2 * np.pi * parameters["fp2_ghz"] * 1e9
    dc, ac = 10 ** (parameters["adc_db"] / 20), 10 ** (parameters["aac_db"] / 20)
    return ac * p2 * (s + dc * p1 / ac) / ((s + p1) * (s + p2))


def hermitian_impulse(frequency, transfer, samples_per_ui=50):
    dt = 25e-12 / samples_per_ui
    n = round(1 / ((frequency[1] - frequency[0]) * dt))
    window = np.i0(np.pi * np.sqrt(1 - (frequency / frequency[-1]) ** 2)) / np.i0(np.pi)
    spectrum = np.zeros(n, complex)
    positive = transfer * window
    spectrum[:len(positive)] = positive
    spectrum[-len(positive) + 1:] = positive[1:][::-1].conj()
    return np.fft.ifft(spectrum), window


def independent_pda(pulse, samples_per_ui=50):
    extended = np.pad(pulse, (0, (-len(pulse)) % samples_per_ui)).reshape(-1, samples_per_ui)
    main_positions = extended.argmax(axis=0)
    main = extended.max(axis=0)
    margin = 2 * main - np.abs(extended).sum(axis=0)
    phase = int(margin.argmax())
    return int(main_positions[phase] * samples_per_ui + phase), float(margin[phase])


def parser_checks():
    matrix = np.array([[.15 + .07j, .42 - .11j], [.71 + .19j, -.12 + .03j]])
    ordered = np.array([matrix[0, 0], matrix[1, 0], matrix[0, 1], matrix[1, 1]])
    maximum_error = 0.0
    with tempfile.TemporaryDirectory(prefix="hw7-reference-") as directory:
        for unit, scale in [("HZ", 1), ("KHZ", 1e3), ("MHZ", 1e6), ("GHZ", 1e9)]:
            for representation in ["RI", "MA", "DB"]:
                if representation == "RI":
                    pairs = np.column_stack((ordered.real, ordered.imag))
                else:
                    magnitude = np.abs(ordered)
                    pairs = np.column_stack((20 * np.log10(magnitude) if representation == "DB" else magnitude, np.angle(ordered, deg=True)))
                lines = ["! synthetic two-port reference", f"# {unit.lower()} S {representation} R 75"]
                values = " ".join(format(v, ".17g") for v in pairs.ravel())
                for frequency in [0, 1e9, 2e9]:
                    row = (format(frequency / scale, ".17g") + " " + values).split()
                    lines.extend([" ".join(row[:5]) + " ! continuation", " ".join(row[5:])])
                fixture = Path(directory) / "reference.s2p"
                fixture.write_text("\n".join(lines), encoding="ascii")
                parsed = hw7.read_touchstone(fixture)
                maximum_error = max(maximum_error, float(np.abs(parsed.s - matrix).max()))
                assert np.array_equal(parsed.frequency_hz, [0, 1e9, 2e9]) and parsed.reference_ohm == 75
    check("touchstone_all_units_formats_ordering_and_continuations", maximum_error < 1e-14, maximum_complex_error=maximum_error, combinations=12)
    raw = np.loadtxt(BASE / "instructions" / "HW7_stripline_50cm.s2p", comments=("!", "#"))
    network = hw7.read_touchstone(BASE / "instructions" / "HW7_stripline_50cm.s2p")
    expected = raw[:, 3] * (np.cos(np.deg2rad(raw[:, 4])) + 1j * np.sin(np.deg2rad(raw[:, 4])))
    check("assignment_frequency_and_s21_preserved", len(raw) == 8001 and np.array_equal(network.frequency_hz, raw[:, 0]) and np.all(np.diff(raw[:, 0]) == 10e6) and raw[-1, 0] == 80e9 and np.allclose(network.s[:, 1, 0], expected, atol=1e-15, rtol=0), points=len(raw), df_hz=float(raw[1, 0] - raw[0, 0]), fmax_hz=float(raw[-1, 0]))
    return network


def ctle_checks(frequency):
    test_sets = [PARAMETERS, dict(fp1_ghz=20, fp2_ghz=40, adc_db=-6, aac_db=9), dict(fp1_ghz=1, fp2_ghz=1, adc_db=6, aac_db=-5)]
    rng = np.random.default_rng(61)
    for _ in range(12):
        test_sets.append(dict(fp1_ghz=rng.uniform(1, 40), fp2_ghz=rng.uniform(1, 40), adc_db=rng.uniform(-20, 6), aac_db=rng.uniform(-5, 15)))
    worst_gain_error = worst_frequency_error = worst_dense_gap = 0.0
    dense_frequency = np.r_[0, np.geomspace(1e3, 1e15, 180000)]
    for parameters in test_sets:
        reported = hw7.ctle_peak(**parameters)
        peak_f, peak_db = independent_peak(parameters)
        worst_gain_error = max(worst_gain_error, abs(reported["peak_gain_db"] - peak_db))
        worst_frequency_error = max(worst_frequency_error, abs(reported["peak_frequency_hz"] - peak_f))
        dense_max = float(20 * np.log10(np.abs(independent_ctle(dense_frequency, parameters))).max())
        worst_dense_gap = max(worst_dense_gap, abs(peak_db - dense_max))
    check("global_peak_dc_case_and_dense_frequency_sweep", worst_gain_error < 1e-10 and worst_frequency_error < .02 and worst_dense_gap < 1e-6, parameter_sets=len(test_sets), maximum_gain_error_db=worst_gain_error, maximum_frequency_error_hz=worst_frequency_error, maximum_dense_gap_db=worst_dense_gap)
    actual = hw7.ctle_peak(**PARAMETERS)
    response_error = float(np.abs(hw7.ctle_response(frequency, **PARAMETERS) - independent_ctle(frequency, PARAMETERS)).max())
    check("chosen_ctle_constraints_and_angular_frequency_reference", actual["peak_gain_db"] < 6 and response_error < 1e-14, **actual, response_complex_error=response_error)


def impulse_checks(network):
    frequency = network.frequency_hz
    response = network.s[:, 1, 0] * independent_ctle(frequency, PARAMETERS)
    time, implementation = hw7.discrete_impulse(frequency, response)
    reference, window = hermitian_impulse(frequency, response)
    error = float(np.abs(implementation - reference.real).max())
    imaginary = float(np.abs(reference.imag).max())
    dc_error = float(abs(implementation.sum() - response[0].real))
    check("irfft_matches_independent_full_hermitian_ifft_and_dc", error < 1e-14 and imaginary < 1e-14 and dc_error < 1e-12 and len(time) == 200000, maximum_kernel_error=error, maximum_reference_imaginary=imaginary, dc_sum_error=dc_error, dt_s=float(time[1]), fft_period_s=float(len(time) * time[1]), kernel_sum=float(implementation.sum()), expected_dc_gain=float(response[0].real), maximum_window_error=float(np.abs(hw7.kaiser_weights(frequency) - window).max()))
    return response, implementation


def source_and_convolution_checks():
    rng = np.random.default_rng(372)
    pattern = np.r_[1, 0, 0, 1, 1, 0, rng.integers(0, 2, 131)].astype(np.uint8)
    untouched = pattern.copy()
    source = source_reference(pattern)
    source_error = float(np.abs(hw7.make_nrz(pattern) - source).max())
    kernel = rng.normal(size=129) * np.exp(-np.arange(129) / 30) / 13
    reference_source = source_reference(np.r_[pattern, 0])
    direct = np.convolve(reference_source, kernel) / 2
    output = hw7.simulate(pattern, kernel)
    error = float(np.abs(output - direct).max())
    identity = hw7.simulate(pattern, [1.0])
    identity_error = float(np.abs(identity - reference_source / 2).max())
    check("source_endpoints_all_bits_and_partial_overlap_add_sections", source_error < 1e-14 and error < 1e-13 and np.array_equal(pattern, untouched) and len(output) == len(reference_source) + len(kernel) - 1, bits=len(pattern), source_error_v=source_error, maximum_direct_convolution_error_v=error, output_samples=len(output), expected_output_samples=len(direct))
    check("physical_matching_divider_exactly_once", identity_error < 1e-14 and source.max() == 2 and identity.max() == 1, identity_error_v=identity_error, source_high_v=float(source.max()), matched_identity_high_v=float(identity.max()))


def eye_and_sensitivity_checks(network, response, full_h):
    frequency = network.frequency_hz
    dt = .5e-12
    cutoff = round(20e-9 / dt)
    h = full_h[:cutoff]
    pulse = hw7.pulse_response(h)
    cursor, pda = independent_pda(pulse)
    bits = np.random.default_rng(0).integers(0, 2, 65539, dtype=np.uint8)
    copy = bits.copy()
    y = hw7.simulate(bits, h)
    info = hw7.analyze_eye(y, bits, pulse, skip_bits=800)
    ids = np.arange(800, len(bits))
    levels = y[ids * 50 + cursor]
    labels = bits[ids]
    random_eh = float(levels[labels == 1].min() - levels[labels == 0].max())
    before_h = hw7.discrete_impulse(frequency, network.s[:, 1, 0])[1][:cutoff]
    before_y = hw7.simulate(bits, before_h)
    finite = bool(np.isfinite(y).all() and np.isfinite(before_y).all())
    check("same_pattern_and_bit_conditioned_random_center_eh", finite and np.array_equal(bits, copy) and len(y) == len(before_y) and cursor == info["cursor_index"] and abs(info["center_eye_height_v"] - random_eh) < 1e-14 and abs(info["pda_eye_height_v"] - pda) < 1e-14 and random_eh >= .3 and pda >= .3, bits=len(bits), analyzed_bits=len(ids), cursor_index=cursor, center_time_s=cursor * dt, random_center_eh_v=random_eh, discrete_all_tap_pda_eh_v=pda, excluded_initial_bits=800, all_waveform_samples_finite=finite)

    continuous_results = []
    for samples_per_ui in [25, 50, 100]:
        step = 25e-12 / samples_per_ui
        size = round(1 / ((frequency[1] - frequency[0]) * step))
        pulse_spectrum = 2 * 25e-12 * np.sinc(frequency * 25e-12) * np.sinc(frequency * 5e-12) * np.exp(-1j * np.pi * frequency * 30e-12)
        window = np.i0(np.pi * np.sqrt(1 - (frequency / frequency[-1]) ** 2)) / np.i0(np.pi)
        hermitian = np.zeros(size, complex)
        positive = response * window * pulse_spectrum / 2
        hermitian[:len(positive)] = positive
        hermitian[-len(positive) + 1:] = positive[1:][::-1].conj()
        analytic_pulse = np.fft.ifft(hermitian).real / step
        analytic_cursor, analytic_eh = independent_pda(analytic_pulse, samples_per_ui)
        continuous_results.append({"dt_s": step, "cursor_time_s": analytic_cursor * step, "full_period_analytic_pda_eh_v": analytic_eh})
    half, quarter = continuous_results[1], continuous_results[2]
    check("independent_continuous_trapezoid_spectrum_and_dt_convergence", abs(half["full_period_analytic_pda_eh_v"] - quarter["full_period_analytic_pda_eh_v"]) < .0002 and min(item["full_period_analytic_pda_eh_v"] for item in continuous_results) > .3, reference="2*Tb*sinc(f*Tb)*sinc(f*Tr)*exp(-j*pi*f*(Tb+Tr)); H=channel*CTLE/2", sensitivity=continuous_results)
    tails = []
    for duration in [20, 40, 80, 100]:
        count = min(round(duration * 1e-9 / dt), len(full_h))
        _, margin = independent_pda(hw7.pulse_response(full_h[:count]))
        tails.append({"impulse_length_ns": duration, "discrete_all_tap_pda_eh_v": margin})
    check("tail_cutoff_sensitivity_retains_300mv_margin", min(item["discrete_all_tap_pda_eh_v"] for item in tails) > .3 and abs(tails[0]["discrete_all_tap_pda_eh_v"] - tails[-1]["discrete_all_tap_pda_eh_v"]) < .003, sensitivity=tails, maximum_cutoff_difference_v=abs(tails[0]["discrete_all_tap_pda_eh_v"] - tails[-1]["discrete_all_tap_pda_eh_v"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE / "results" / "independent_verification.json")
    args = parser.parse_args()
    network = parser_checks()
    ctle_checks(network.frequency_hz)
    response, full_h = impulse_checks(network)
    source_and_convolution_checks()
    eye_and_sensitivity_checks(network, response, full_h)
    record = {"passed": all(item["passed"] for item in CHECKS), "all_passed": all(item["passed"] for item in CHECKS), "checks": CHECKS, "parameters": PARAMETERS, "scope": "Finite seeded random eye and conservative all-tap PDA are separately reported; neither is substituted for the other.", "report_checks": [
        {"label": "Independent reader and gain checks", "detail": "All 12 unit/format combinations reproduce a defined two-port matrix. A dense 12-decade frequency sweep and a separately derived peak calculation confirm the global CTLE gain limit."},
        {"label": "Impulse and physical normalization", "detail": "The real IFFT matches a full Hermitian IFFT and preserves the DC kernel sum. An identity transfer maps the 2 V source high level to 1 V with exactly one matched-source division by two."},
        {"label": "Direct versus sectional convolution", "detail": "Overlap-add matches direct linear convolution for a 137-bit pattern, including its partial final block and final falling ramp; transition endpoints independently match time-knot interpolation."},
        {"label": "Eye alignment and sensitivity", "detail": "Bit-conditioned random EH and all-tap PDA use y[k*samples_per_ui+cursor] against bit k and remain distinct measures. An independent continuous-pulse spectrum, 0.25 ps grid and full 100 ns impulse period retain more than 300 mV worst-case margin."}
    ]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": record["passed"], "checks": len(CHECKS)}, indent=2))


if __name__ == "__main__":
    main()
