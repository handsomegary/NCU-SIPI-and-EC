"""Independently verify HW14 saved records, decisions and physical scaling.

The verifier does not import the production solver to recompute results.
It uses complete saved records, a direct full complex IDFT, analytic ramps,
explicit FFE shifts, recovered-decision convolution, and Gray bit labels.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import unittest

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import oaconvolve


ROOT = Path(__file__).resolve().parents[1]
GRAY_LEVELS = np.array([-1., -1 / 3, 1 / 3, 1.])
GRAY_LABELS = np.array([[0, 0], [0, 1], [1, 1], [1, 0]], dtype=np.uint8)


def width_from_boundaries(phase, lower, upper, threshold):
    """Find the strict-positive threshold interval containing phase zero."""
    center = int(np.argmin(np.abs(phase)))
    clearance = np.minimum(lower - threshold, threshold - upper)
    if clearance[center] <= 0:
        return 0.
    good = np.flatnonzero(clearance > 0)
    component = np.split(good, np.flatnonzero(np.diff(good) != 1) + 1)
    selected = next(group for group in component if center in group)
    first, last = int(selected[0]), int(selected[-1])
    left, right = float(phase[first]), float(phase[last])
    if first:
        left = float(np.interp(0., clearance[first - 1:first + 1], phase[first - 1:first + 1]))
    if last + 1 < len(phase):
        right = float(np.interp(0., clearance[last:last + 2][::-1], phase[last:last + 2][::-1]))
    return right - left


def direct_idft(frequency, positive, dt):
    """Explicitly populate positive/negative frequencies, then a full IFFT."""
    n = int(round(1 / ((frequency[1] - frequency[0]) * dt)))
    spectrum = np.zeros(n, dtype=complex)
    spectrum[:len(positive)] = positive
    spectrum[0] = positive[0].real
    for index in range(1, len(positive)):
        if index == n - index:
            spectrum[index] = positive[index].real
        else:
            spectrum[n - index] = np.conj(positive[index])
    full = np.fft.ifft(spectrum)
    return full.real, float(np.max(np.abs(full.imag)))


def verify(results, output, require_validation=True):
    results, output = Path(results), Path(output)
    result = json.loads((results / "metrics.json").read_text(encoding="utf-8"))
    checks, evidence, case_checks, failures = {}, {}, {}, []
    archive_hashes = {}

    def check(name, condition, detail=None):
        checks[name] = bool(condition)
        if detail is not None:
            evidence[name] = detail
        if not condition:
            failures.append(name)

    def close(name, actual, expected, tolerance=1e-11):
        actual, expected = np.asarray(actual), np.asarray(expected)
        error = float(np.max(np.abs(actual - expected))) if actual.size else 0.
        check(name, actual.shape == expected.shape and np.isfinite(actual).all() and error <= tolerance,
              {"maximum_absolute_error": error, "tolerance": tolerance})
        return error

    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"), pattern="test_*.py")
    unit_stream = io.StringIO()
    unit_result = unittest.TextTestRunner(stream=unit_stream, verbosity=1).run(suite)
    check("independent_unit_tests", unit_result.wasSuccessful(), {"tests_run": unit_result.testsRun,
          "failures": len(unit_result.failures), "errors": len(unit_result.errors)})
    if not unit_result.wasSuccessful():
        evidence["unit_test_output"] = unit_stream.getvalue()

    parameters = result["parameters"]
    for key, value in {"source_vpp_v": 1., "rise_time_ps": 10., "sj_amplitude_ps": 1.5,
                       "sj_frequency_mhz": 10., "rj_sigma_ps": 1., "ffe_taps": 3, "dfe_taps": 5,
                       "ctle_maximum_absolute_gain_db": 18.}.items():
        check("parameter_" + key, parameters[key] == value)
    channel_path = ROOT / "instructions" / result["channel"]["input_filename"]
    check("input_sha256", hashlib.sha256(channel_path.read_bytes()).hexdigest() == result["channel"]["input_sha256"])
    # The supplied file is GHz/MA; read its numeric rows independently.
    text = channel_path.read_text(encoding="utf-8-sig")
    numbers = []
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].strip()
        if line and not line.startswith("#"):
            numbers.extend(float(token) for token in line.split())
    source = np.asarray(numbers).reshape(-1, 33)
    source_frequency = source[:, 0] * 1e9
    source_pairs = source[:, 1:].reshape(-1, 4, 4, 2)
    scattering = source_pairs[..., 0] * np.exp(1j * np.pi / 180 * source_pairs[..., 1])
    incident = np.array([1., 0., -1., 0.]) / np.sqrt(2)
    outgoing = scattering @ incident
    source_odd = (outgoing[:, 1] - outgoing[:, 3]) / np.sqrt(2)
    with np.load(results / "channel.npz") as saved_channel:
        frequency = saved_channel["frequency_hz"].copy()
        odd = saved_channel["sdd21"].copy()
    close("source_frequency_hz", frequency, source_frequency, 1e-4)
    close("source_odd_mode", odd, source_odd, 2e-15)
    check("channel_grid", len(frequency) == 8001 and frequency[0] == 0 and frequency[-1] == 80e9
          and np.allclose(np.diff(frequency), 10e6, rtol=0, atol=1e-3))

    cache = {}
    for key, case in result["cases"].items():
        archive_hashes[case["archive"]] = hashlib.sha256((results / case["archive"]).read_bytes()).hexdigest()
        with np.load(results / case["archive"]) as archive:
            data = {name: archive[name].copy() for name in archive.files}
        initial_failures = len(failures)
        prefix = key + "/"
        modulation, rate = case["modulation"], case["bit_rate_gbps"]
        stride = 1 if modulation == "nrz" else 2
        p, dt, clock = int(data["samples_per_ui"]), float(data["dt_s"]), int(data["clock_sample"])
        start, count = int(data["retained_start_symbol"]), int(data["retained_symbol_count"])
        retained = slice(start, start + count)
        truth, bits = data["symbols"], data["bits"]
        decisions, feedback, corrected = data["decisions"], data["feedback_v"], data["decision_corrected_v"]
        main = float(data["main_cursor_v"])
        levels = np.array([-1., 1.]) if stride == 1 else GRAY_LEVELS
        thresholds = main * (levels[:-1] + levels[1:]) / 2
        metric = case["metrics"]
        check(prefix + "record_metadata", count * stride == case["measured_bits"] == metric["measured_bits"] == parameters["bits"]
              and count == case["measured_symbols"] == metric["measured_symbols"]
              and p == parameters["samples_per_ui"] and len(bits) == len(truth) * stride)
        close(prefix + "physical_ui_ps", p * dt * 1e12, stride * 1000 / rate, 1e-10)
        close(prefix + "clock_ps", clock * dt * 1e12, case["design"]["clock_time_ps"], 1e-10)
        expected_truth = 2 * bits.astype(float) - 1 if stride == 1 else np.array([-1., -1/3, 1., 1/3])[2 * bits[::2] + bits[1::2]]
        close(prefix + "payload_encoding", truth, expected_truth, 1e-15)
        taps = data["ffe_taps"]
        ffe_expected = np.zeros(len(truth))
        for offset, coefficient in enumerate(taps):
            ffe_expected[offset:] += coefficient * truth[:len(truth) - offset]
        close(prefix + "ffe_on_symbols", data["ffe_symbols"], ffe_expected, 2e-15)
        check(prefix + "ffe_l1_and_source_range", len(taps) == 3 and abs(np.abs(taps).sum() - 1) < 1e-12
              and np.max(np.abs(ffe_expected)) <= 1 + 1e-12)
        nominal = np.arange(len(truth)) * p * dt
        expected_sj = 1.5e-12 * np.sin(2 * np.pi * 10e6 * nominal)
        close(prefix + "sinusoidal_jitter_seconds", data["sj_s"], expected_sj, 1e-21)
        close(prefix + "jittered_boundaries_seconds", data["boundary_times_s"], nominal + data["rj_s"] + data["sj_s"], 1e-20)
        expected_rj = np.random.default_rng(parameters["seed"] + 100).normal(0, 1e-12, len(truth) * stride)[::stride]
        close(prefix + "reproducible_absolute_rj_seconds", data["rj_s"], expected_rj, 1e-25)
        close(prefix + "retained_rj_std_ps", np.std(data["rj_s"][retained]) * 1e12, case["jitter"]["measured_rj_std_ps"], 1e-12)
        check(prefix + "ramp_spacing", np.min(np.diff(data["boundary_times_s"])) > 10e-12)

        settings = case["equalizers"]["ctle"]
        zero = settings["zero_ghz"] * 1e9
        pole1, pole2 = np.asarray(settings["poles_ghz"]) * 1e9
        raw_ctle = lambda f: (1 + 1j * np.asarray(f) / zero) / ((1 + 1j * np.asarray(f) / pole1) * (1 + 1j * np.asarray(f) / pole2))
        maximum = minimize_scalar(lambda log_f: -abs(raw_ctle(10 ** log_f)), bounds=(3, 14), method="bounded", options={"xatol": 1e-12})
        raw_peak = max(1., -maximum.fun)
        expected_ctle = raw_ctle(frequency) * 10 ** (18 / 20) / raw_peak
        close(prefix + "ctle_reported_corners_and_peak", data["ctle"], expected_ctle, 8e-12)
        check(prefix + "equalizer_complexity", len(settings["poles_ghz"]) == 2 and settings["peak_gain_db"] == 18
              and case["equalizers"]["ffe"]["tap_count"] == 3 and case["equalizers"]["dfe"]["tap_count"] == 5)

        # Independent full complex IFFT verifies gain, DC, phase, and units.
        kernel_errors = {}
        for archive_name, response in (("h_channel", odd / 2), ("h_equalized", odd * data["ctle"] / 2)):
            full, imaginary_error = direct_idft(frequency, response, dt)
            kernel_errors[archive_name] = close(prefix + archive_name, data[archive_name], full[:len(data[archive_name])], 3e-13)
            check(prefix + archive_name + "_idft_real", imaginary_error < 2e-14)
            close(prefix + archive_name + "_full_dc", np.sum(full), response[0].real, 5e-13)
        # The final full IDFT here is the equalized transfer.  Reconstruct
        # the entire transmitted record independently, in bounded chunks,
        # then use overlap-add convolution rather than the solver's FFT
        # convolution.  Startup, guards and every payload sample are kept.
        oracle_kernel = full[:len(data["h_equalized"])].copy()
        ui = stride / (rate * 1e9)
        source_count = int(np.ceil(len(truth) * ui / dt))
        entire_source = np.empty(source_count)
        boundaries = data["boundary_times_s"]
        voltage_levels = ffe_expected / 2
        previous_levels = np.r_[0., voltage_levels[:-1]]
        for block_start in range(0, source_count, 2 ** 18):
            block_stop = min(source_count, block_start + 2 ** 18)
            times = np.arange(block_start, block_stop, dtype=float) * dt
            event = np.searchsorted(boundaries, times, side="right") - 1
            selected = np.maximum(event, 0)
            fraction = np.clip((times - boundaries[selected]) / 10e-12, 0, 1)
            block = previous_levels[selected] + (voltage_levels[selected] - previous_levels[selected]) * fraction
            block[event < 0] = 0.
            entire_source[block_start:block_stop] = block
        oracle_received = oaconvolve(entire_source, oracle_kernel, mode="full")
        all_decision_positions = clock + np.arange(len(truth)) * p
        end_to_end_error = close(prefix + "end_to_end_every_decision_input_v", data["decision_input_v"],
                                 oracle_received[all_decision_positions], 8e-11)
        oracle_minima = np.empty_like(data["level_min_v"])
        oracle_maxima = np.empty_like(data["level_max_v"])
        oracle_histogram = np.empty_like(data["histogram"])
        retained_index = np.arange(start, start + count)
        groups = [truth[retained] == level for level in levels]
        edges = data["voltage_edges_v"]
        for column, phase in enumerate(range(-p // 2, p // 2 + 1)):
            voltage = oracle_received[clock + retained_index * p + phase] - feedback[retained_index + int(phase == p // 2)]
            for row, group in enumerate(groups):
                oracle_minima[row, column] = np.min(voltage[group])
                oracle_maxima[row, column] = np.max(voltage[group])
            bins = np.searchsorted(edges, voltage, side="right") - 1
            bins[voltage == edges[-1]] = len(edges) - 2
            in_range = (bins >= 0) & (bins < len(edges) - 1)
            oracle_histogram[:, column] = np.bincount(bins[in_range], minlength=len(edges) - 1)
        full_envelope_error = max(close(prefix + "end_to_end_all_phase_minima_v", data["level_min_v"], oracle_minima, 8e-11),
                                  close(prefix + "end_to_end_all_phase_maxima_v", data["level_max_v"], oracle_maxima, 8e-11))
        histogram_differences = int(np.count_nonzero(oracle_histogram != data["histogram"]))
        histogram_l1_difference = int(np.abs(oracle_histogram.astype(np.int64) - data["histogram"].astype(np.int64)).sum())
        check(prefix + "end_to_end_every_histogram_cell", histogram_differences == 0,
              {"differing_cells": histogram_differences, "l1_count_difference": histogram_l1_difference,
               "tested_cells": int(oracle_histogram.size)})
        time = np.arange(len(data["source_pulse_v"])) * dt
        pulse = .5 * (np.clip(time / 10e-12, 0, 1) - np.clip((time - p * dt) / 10e-12, 0, 1))
        close(prefix + "source_pulse_v", data["source_pulse_v"], pulse, 2e-15)
        for kernel_name, pulse_name in (("h_channel", "sbr_channel_v"), ("h_equalized", "sbr_equalized_v")):
            close(prefix + pulse_name, data[pulse_name], np.convolve(data[kernel_name], pulse), 3e-13)
        sbr_ffe = np.zeros(len(data["sbr_equalized_v"]) + (len(taps) - 1) * p)
        for offset, coefficient in enumerate(taps):
            sbr_ffe[offset * p:offset * p + len(data["sbr_equalized_v"])] += coefficient * data["sbr_equalized_v"]
        close(prefix + "ffe_single_symbol_response", data["sbr_ffe_v"], sbr_ffe, 3e-13)
        close(prefix + "main_cursor_v", main, sbr_ffe[clock], 3e-13)
        close(prefix + "dfe_postcursors_v", data["dfe_taps_v"], sbr_ffe[clock + np.arange(1, 6) * p], 3e-13)
        close(prefix + "threshold_scale_v", data["thresholds_v"], thresholds, 3e-13)

        expected_feedback = np.convolve(decisions, np.r_[0., data["dfe_taps_v"]])[:len(decisions)]
        close(prefix + "recovered_decision_feedback_v", feedback, expected_feedback, 3e-13)
        close(prefix + "corrected_decision_voltage_v", corrected, data["decision_input_v"] - feedback, 3e-13)
        expected_decisions = levels[np.searchsorted(thresholds, corrected, side="right")]
        close(prefix + "slicer_decisions", decisions, expected_decisions, 1e-15)
        decoded = (decisions[retained] > 0).astype(np.uint8) if stride == 1 else GRAY_LABELS[np.argmin(np.abs(decisions[retained, None] - GRAY_LEVELS), axis=1)].ravel()
        bit_errors = int(np.count_nonzero(decoded != bits[start * stride:(start + count) * stride]))
        symbol_errors = int(np.count_nonzero(decisions[retained] != truth[retained]))
        check(prefix + "bit_and_symbol_error_counts", bit_errors == metric["bit_errors"] and symbol_errors == metric["symbol_errors"])
        close(prefix + "empirical_ber", bit_errors / (count * stride), metric["empirical_ber"], 1e-18)
        close(prefix + "empirical_ser", symbol_errors / count, metric["empirical_ser"], 1e-18)
        minima = np.array([np.min(corrected[retained][truth[retained] == level]) for level in levels])
        maxima = np.array([np.max(corrected[retained][truth[retained] == level]) for level in levels])
        close(prefix + "center_conditional_min_v", data["level_min_v"][:, p // 2], minima, 3e-13)
        close(prefix + "center_conditional_max_v", data["level_max_v"][:, p // 2], maxima, 3e-13)
        heights = minima[1:] - maxima[:-1]
        close(prefix + "finite_record_eye_heights_v", metric["eye_heights_v"], heights, 3e-13)
        close(prefix + "worst_eye_height_v", metric["worst_eye_height_v"], np.min(heights), 3e-13)
        expected_phase = np.arange(-p // 2, p // 2 + 1) * dt * 1e12
        close(prefix + "phase_axis_ps", data["phase_ps"], expected_phase, 1e-10)
        widths = [width_from_boundaries(data["phase_ps"], data["level_min_v"][i + 1], data["level_max_v"][i], threshold)
                  for i, threshold in enumerate(thresholds)]
        close(prefix + "fixed_threshold_eye_widths_ps", metric["eye_widths_ps"], widths, 1e-9)
        close(prefix + "worst_eye_width_ps", metric["worst_eye_width_ps"], np.min(widths), 1e-9)
        close(prefix + "normalized_eye_width_ui", metric["eye_widths_ui"], np.asarray(widths) / (p * dt * 1e12), 1e-11)
        check(prefix + "exact_histogram_exposure", np.all(data["histogram"].sum(axis=0) == count))
        check(prefix + "all_three_or_one_eyes_at_same_clock", len(heights) == len(thresholds) == (1 if stride == 1 else 3))

        # Reconstruct the saved source excerpt from FFE voltages and current
        # jittered boundaries.  This catches an unintended analog FFE delay.
        source_positions = (start + 200) * p + np.rint(data["source_time_relative_ps"] / (dt * 1e12)).astype(int)
        source_times = source_positions * dt
        event = np.searchsorted(data["boundary_times_s"], source_times, side="right") - 1
        old = np.where(event > 0, ffe_expected[np.maximum(0, event - 1)], 0.) / 2
        fraction = np.clip((source_times - data["boundary_times_s"][event]) / 10e-12, 0, 1)
        source_excerpt = old + (ffe_expected[event] / 2 - old) * fraction
        source_excerpt_error = close(prefix + "jittered_ffe_source_excerpt_v", data["source_waveform_v"], source_excerpt, 8e-11)
        shot_center = clock + (start + 200) * p
        shot_positions = shot_center + np.rint(data["waveform_time_relative_ps"] / (dt * 1e12)).astype(int)
        hold = (shot_positions - clock + p // 2) // p
        close(prefix + "coherent_dfe_waveform_hold_v", data["waveform_after_dfe_v"], data["waveform_before_dfe_v"] - feedback[hold], 3e-13)
        at_clock = (shot_positions - clock) % p == 0
        close(prefix + "plotted_decision_voltage_v", data["waveform_after_dfe_v"][at_clock], corrected[(shot_positions[at_clock] - clock) // p], 3e-13)
        cache[key] = {"bits": bits, "rj": data["rj_s"], "sj": data["sj_s"], "stride": stride}
        case_checks[key] = {"successful": len(failures) == initial_failures, "bit_rate_gbps": rate,
                            "modulation": modulation, "measured_bits": count * stride, "measured_symbols": count,
                            "samples_per_ui": p, "dt_ps": dt * 1e12, "bit_errors": bit_errors,
                            "symbol_errors": symbol_errors, "eye_heights_v": heights.tolist(),
                            "eye_widths_ps": widths, "kernel_max_error_v_per_v": kernel_errors,
                            "source_excerpt_max_error_v": source_excerpt_error,
                            "end_to_end_source_samples": source_count,
                            "end_to_end_decisions_checked": int(len(truth)),
                            "end_to_end_max_decision_voltage_error_v": end_to_end_error,
                            "end_to_end_max_conditional_envelope_error_v": full_envelope_error,
                            "end_to_end_histogram_cells_checked": int(oracle_histogram.size)}

    rates = sorted({case["bit_rate_gbps"] for case in result["cases"].values()})
    for rate in rates:
        pair = [(key, case) for key, case in result["cases"].items() if case["bit_rate_gbps"] == rate]
        by_modulation = {case["modulation"]: cache[key] for key, case in pair}
        nrz, pam4 = by_modulation["nrz"], by_modulation["pam4"]
        check(f"{rate:g}gbps/same_bit_record", np.array_equal(nrz["bits"], pam4["bits"]))
        close(f"{rate:g}gbps/shared_rj_clock_seconds", nrz["rj"][::2], pam4["rj"], 1e-25)
        close(f"{rate:g}gbps/shared_sj_clock_seconds", nrz["sj"][::2], pam4["sj"], 1e-25)

    validation = result.get("validation", {})
    check("validation_records_available", bool(validation) or not require_validation)
    validation_summary = {}
    for variant in ("independent_seed", "time_refinement", "impulse_extension"):
        differences = []
        for key, case in result["cases"].items():
            if key not in validation:
                continue
            item = validation[key][variant]
            refined = item["metrics"]
            baseline = case["metrics"]
            check(key + "/" + variant + "_finite_counts", refined["measured_bits"] == case["measured_bits"]
                  and refined["measured_symbols"] == case["measured_symbols"]
                  and np.isfinite(refined["eye_heights_v"]).all() and np.isfinite(refined["eye_widths_ps"]).all())
            if variant == "time_refinement":
                check(key + "/time_refinement_grid", item["samples_per_ui"] == parameters["samples_per_ui"] * 2)
            if variant == "impulse_extension":
                check(key + "/impulse_extension_duration", item["impulse_ns"] > parameters["impulse_ns"] and item["impulse_ns"] < 100)
            differences.append({"case": key,
                                "eye_height_change_mv": (refined["worst_eye_height_v"] - baseline["worst_eye_height_v"]) * 1000,
                                "eye_width_change_ps": refined["worst_eye_width_ps"] - baseline["worst_eye_width_ps"],
                                "bit_errors": refined["bit_errors"]})
        validation_summary[variant] = differences
    comparison_sensitivity = []
    for rate in rates:
        matching = {case["modulation"]: key for key, case in result["cases"].items() if case["bit_rate_gbps"] == rate}
        for variant in ("baseline", "independent_seed", "time_refinement", "impulse_extension"):
            if variant != "baseline" and not all(key in validation for key in matching.values()):
                continue
            metrics = {mod: (result["cases"][key]["metrics"] if variant == "baseline" else validation[key][variant]["metrics"])
                       for mod, key in matching.items()}
            difference = metrics["pam4"]["worst_eye_height_v"] - metrics["nrz"]["worst_eye_height_v"]
            width_difference = metrics["pam4"]["worst_eye_width_ps"] - metrics["nrz"]["worst_eye_width_ps"]
            comparison_sensitivity.append({"bit_rate_gbps": rate, "variant": variant,
                                           "pam4_minus_nrz_minimum_eye_height_mv": difference * 1000,
                                           "pam4_minus_nrz_minimum_eye_width_ps": width_difference,
                                           "larger_minimum_eye_height": "pam4" if difference > 0 else "nrz",
                                           "larger_minimum_eye_width_ps": "pam4" if width_difference > 0 else "nrz"})
    successful = all(checks.values())
    report = {"successful": successful, "passed": successful, "unit_test_count": unit_result.testsRun,
              "check_count": len(checks), "failed_checks": failures, "checks": checks, "evidence": evidence,
              "cases": case_checks, "validation": validation_summary, "comparison_sensitivity": comparison_sensitivity,
              "verified_file_sha256": {"metrics.json": hashlib.sha256((results / "metrics.json").read_bytes()).hexdigest(),
                                       "src/hw14.py": hashlib.sha256((ROOT / "src" / "hw14.py").read_bytes()).hexdigest(),
                                       **archive_hashes},
              "limitations": ["Finite-record bit-error counts are not a low-BER guarantee.",
                              "Saved validation metrics are consistency-checked; their full waveforms are not persisted.",
                              "Eye widths use the documented fixed-threshold DFE hold convention."]}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"successful": successful, "unit_test_count": unit_result.testsRun,
                      "check_count": len(checks), "failed_checks": failures,
                      "validation": validation_summary}, indent=2))
    return 0 if successful else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "verification.json")
    parser.add_argument("--allow-no-validation", action="store_true", help="For smoke records only; full report requires validation")
    args = parser.parse_args()
    return verify(args.results, args.output, not args.allow_no_validation)


if __name__ == "__main__":
    raise SystemExit(main())
