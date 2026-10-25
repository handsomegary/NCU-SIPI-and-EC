"""Independently verify every persisted HW9 decision and waveform snapshot.

The oracle imports no production solver/channel module. It reads the supplied
MA Touchstone directly, uses NumPy's inverse FFT and direct pulse convolution,
then reconstructs source-symbol responses with polyphase convolution. Scalar
history feedback starts at the first clock, including all pre-payload history.
Local waveform samples are rebuilt from pulse superposition and held DAC age.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
VOLTAGE_TOL = 2e-10
TIME_PS_TOL = 1e-9
SOURCE_SHA256 = "564560062e750405e89cadbb8be73cc7a40ad3188d19383ebbfe1e22723f5efc"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def close(actual, expected, label, tolerance=VOLTAGE_TOL):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape, f"{label}: shape mismatch")
    require(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)),
            f"{label}: nonfinite value")
    difference = float(np.max(np.abs(actual - expected))) if actual.size else 0.0
    require(difference <= tolerance, f"{label}: maximum difference {difference:.6g} exceeds {tolerance:.6g}")
    return difference


def exact(actual, expected, label):
    require(np.array_equal(actual, expected), f"{label}: exact values differ")


def csv_values(folder, filename, header):
    require(Path(filename).name == filename, "Artifact filename must be relative to results")
    path = folder / filename
    with path.open(encoding="utf-8-sig") as stream:
        actual_header = stream.readline().strip().split(",")
    exact(actual_header, header, filename + " header")
    result = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    require(result.shape[1] == len(header) and np.all(np.isfinite(result)),
            filename + ": invalid numeric CSV")
    return result


def independent_channel(data, model):
    content = data.read_bytes()
    require(hashlib.sha256(content).hexdigest() == SOURCE_SHA256,
            "Channel data does not match the audited original")
    options = next(line.strip().lower().split() for line in content.decode("ascii").splitlines()
                   if line.lstrip().startswith("#"))
    require(options == ["#", "hz", "s", "ma", "r", "50"],
            "Independent reader expects the supplied Hz/MA/50-ohm two-port input")
    records = np.loadtxt(data, comments=("!", "#"))
    require(records.shape == (8001, 9), "Unexpected channel record count")
    frequency = records[:, 0]
    exact(np.diff(frequency), np.full(8000, 10e6), "Input frequency grid")
    require(frequency[0] == 0 and frequency[-1] == 80e9, "Unexpected measured frequency range")
    dt = model["sample_step_ps"] * 1e-12
    p = int(model["samples_per_ui"])
    rise = round(model["rise_time_ps"] * 1e-12 / dt)
    require(p == 50 and rise == 10, "Assignment requires 50 samples/UI and a 10-sample rise")
    require(model["channel_transfer"] == "S21/2" and model["decision_alphabet"] == [0, 1],
            "Incorrect channel or slicer convention")
    close(model["data_rate_gbps"], 40, "Data rate")
    close(model["ui_ps"], 25, "Unit interval")
    close(model["source_vpp_v"], 2, "Source swing")
    close(model["reference_impedance_ohm"], 50, "Reference impedance")
    transfer = records[:, 3] * np.exp(1j * np.deg2rad(records[:, 4])) / 2
    n = round(1 / (10e6 * dt))
    require(n == 200000, "Unexpected inverse-transform period")
    window = np.i0(np.pi * model["kaiser_alpha"] * np.sqrt(1 - (frequency / frequency[-1]) ** 2)) / np.i0(np.pi * model["kaiser_alpha"])
    spectrum = np.zeros(n // 2 + 1, dtype=complex)
    spectrum[:len(frequency)] = transfer * window
    spectrum[0] = spectrum[0].real
    spectrum[-1] = spectrum[-1].real
    period_impulse = np.fft.irfft(spectrum, n=n)
    count = round(model["impulse_truncation_ns"] * 1e-9 / dt)
    impulse = period_impulse[:count]
    sample = np.arange(p + rise + 1)
    pulse = model["source_vpp_v"] * np.maximum(0, np.minimum.reduce([
        sample / rise, np.ones(len(sample)), (p + rise - sample) / rise]))
    response = np.convolve(pulse, impulse)
    return response, impulse, period_impulse, transfer, pulse, dt


def bounds(response, p, latency, ntaps):
    # A polyphase branch contains all cursors in chronological source latency.
    branch = response[latency % p::p].copy()
    main_index = latency // p
    main = float(branch[main_index])
    selected = branch[main_index + 1:main_index + ntaps + 1]
    taps = -np.pad(selected, (0, max(0, ntaps - len(selected))))
    branch[main_index:main_index + ntaps + 1] = 0
    low = float(sum(value for value in branch if value > 0))
    high = main + float(sum(value for value in branch if value < 0))
    return {"tap_count": ntaps, "latency_samples": latency, "phase_samples": latency % p,
            "main_cursor_v": main, "taps_v": taps,
            "threshold_v": (low + high) / 2, "low_upper_v": low,
            "high_lower_v": high, "worst_case_eye_v": high - low,
            "residual_abs_sum_v": float(np.abs(branch).sum()),
            "residual_sum_v": float(branch.sum()),
            "cursor_offset_min": -main_index,
            "cursor_offset_max": len(branch) - 1 - main_index}


def history_feedback(samples, taps, threshold, force_flip=None):
    decisions = np.empty(len(samples), dtype=np.uint8)
    summer = np.empty(len(samples), dtype=float)
    for clock, raw in enumerate(samples):
        voltage = float(raw)
        for age, coefficient in enumerate(taps, 1):
            if clock >= age:
                voltage += coefficient * decisions[clock - age]
        summer[clock] = voltage
        bit = int(voltage >= threshold)
        decisions[clock] = 1 - bit if clock == force_flip else bit
    return summer, decisions


def metric_reference(samples, truth, threshold, decisions):
    wrong = np.flatnonzero(decisions != truth)
    low = float(max(samples[truth == 0]))
    high = float(min(samples[truth == 1]))
    return {"evaluated_bits": len(truth), "errors": len(wrong),
            "observed_error_fraction": len(wrong) / len(truth),
            "observed_low_max_v": low, "observed_high_min_v": high,
            "observed_eye_v": high - low,
            "minimum_signed_margin_v": float(min(np.where(truth == 1, samples - threshold, threshold - samples))),
            "first_error_payload_bit": int(wrong[0]) if len(wrong) else None}


def verify_metrics(saved, reference, label):
    maximum = 0.0
    for key, expected in reference.items():
        if key in {"evaluated_bits", "errors", "first_error_payload_bit"}:
            require(saved[key] == expected, f"{label}: incorrect {key}")
        else:
            maximum = max(maximum, close(saved[key], expected, label + " " + key))
    return maximum


def local_channel_waveform(bits, response, p, indices):
    # Direct local pulse superposition, independent of production overlap-add.
    first = max(0, int(np.ceil((indices[0] - len(response) + 1) / p)))
    last = min(len(bits) - 1, int(indices[-1] // p))
    source = np.arange(first, last + 1)
    lag = indices[:, None] - source[None, :] * p
    valid = (lag >= 0) & (lag < len(response))
    coefficient = response[np.clip(lag, 0, len(response) - 1)] * valid
    return coefficient @ bits[source].astype(float)


def snapshot_reference(folder, info, bits, response, p, dt, phase, taps,
                       decisions, scored, truth, payload_decisions):
    values = csv_values(folder, info["csv"], ["time_ps", "post_channel_v", "post_dfe_v"])
    require(values.shape == (2001, 3) and info["samples"] == 2001,
            info["csv"] + ": a 1 ns snapshot must have 2001 endpoint-inclusive samples")
    close(values[:, 0], np.arange(2001) * dt * 1e12,
          info["csv"] + " relative time", TIME_PS_TOL)
    close(info["duration_ns"], 1.0, info["csv"] + " duration")
    start = int(info["start_sample"])
    require(start == info["start_sample"], "Snapshot start must be an integer sample")
    close(info["start_time_s"], start * dt, info["csv"] + " start time", 1e-18)
    indices = start + np.arange(2001)
    raw = local_channel_waveform(bits, response, p, indices)
    held = raw.copy()
    latest = (indices - phase - p // 2) // p
    for age, coefficient in enumerate(taps):
        clock = latest - age
        valid = (clock >= 0) & (clock < len(decisions))
        held[valid] += coefficient * decisions[clock[valid]]
    raw_error = close(values[:, 1], raw, info["csv"] + " direct post-channel waveform")
    corrected_error = close(values[:, 2], held, info["csv"] + " held-feedback waveform")
    in_snapshot = (scored >= indices[0]) & (scored <= indices[-1])
    require(int(in_snapshot.sum()) == info["scored_decisions_in_snapshot"],
            info["csv"] + ": incorrect scored clock count")
    errors = int(np.count_nonzero(payload_decisions[in_snapshot] != truth[in_snapshot]))
    require(errors == info["errors_in_snapshot"], info["csv"] + ": incorrect snapshot error count")
    image = folder / info["file"]
    require(image.is_file() and image.stat().st_size > 1000, "Missing waveform figure")
    return {"file": info["csv"], "samples_verified": 2001,
            "post_channel_max_abs_v": raw_error, "post_dfe_max_abs_v": corrected_error,
            "errors_verified": errors}


def verify(results=None):
    folder = Path(results) if results else HERE / "results"
    summary = json.loads((folder / "summary.json").read_text(encoding="utf-8"))
    model = summary["model"]
    response, impulse, period_impulse, transfer, pulse, dt = independent_channel(
        HERE / "data" / "HW7_stripline_50cm.s2p", model)
    require(summary["verification"]["input_data_sha256"] == SOURCE_SHA256, "Summary input hash mismatch")
    close(summary["verification"]["finite_fir_dc_gain"], impulse.sum(), "Finite FIR DC gain")
    close(summary["verification"]["frequency_dc_gain"], transfer[0].real, "Frequency DC gain")
    close(summary["verification"]["finite_fir_dc_difference"], impulse.sum() - transfer[0].real,
          "Finite FIR DC difference")
    p = int(model["samples_per_ui"])
    count, guard = int(model["payload_bits"]), int(model["guard_bits_each_side"])
    require(count >= 1024 and guard >= int(np.ceil(len(response) / p)), "Insufficient payload or guards")
    bits = np.random.default_rng(int(model["seed"])).integers(0, 2, count + 2 * guard, dtype=np.uint8)
    truth = bits[guard:guard + count]
    require(np.any(truth == 0) and np.any(truth == 1), "Payload must include both logic levels")
    cases = summary["cases"]
    exact([case["design"]["tap_count"] for case in cases], [1, 2, 3], "Required DFE cases")
    report = {"passed": False, "source_sha256": SOURCE_SHA256,
              "oracle": "Independent NumPy inverse FFT, direct pulse convolution, source-bit polyphase convolution, scalar history decisions, and local held-feedback reconstruction; no production solver imports",
              "tolerances": {"voltage_abs_v": VOLTAGE_TOL, "time_abs_ps": TIME_PS_TOL},
              "payload_bits_per_case": count, "pre_payload_history_verified": True,
              "cases": []}
    references = {}
    search_start = min(max(int(np.argmax(response)) - p // 2, 0), len(response) - p)
    response40 = np.convolve(pulse, period_impulse[:round(40e-9 / dt)])
    for case in cases:
        design = case["design"]
        ntaps = int(design["tap_count"])
        candidates = [bounds(response, p, latency, ntaps) for latency in range(search_start, search_start + p)]
        expected = candidates[int(np.argmax([row["worst_case_eye_v"] for row in candidates]))]
        for key, value in expected.items():
            close(design[key], value, f"{ntaps}-tap design {key}")
        phase, latency = expected["phase_samples"], expected["latency_samples"]
        taps, threshold = expected["taps_v"], expected["threshold_v"]
        raw = np.convolve(bits.astype(float), response[phase::p])
        summer, decisions = history_feedback(raw, taps, threshold)
        locations = guard + np.arange(count) + latency // p
        scored = (guard + np.arange(count)) * p + latency
        raw_payload, corrected = raw[locations], summer[locations]
        payload_decisions = decisions[locations]
        values = csv_values(folder, f"decisions_{ntaps}tap.csv",
                            ["payload_bit", "truth", "decision", "post_channel_v", "post_dfe_v", "threshold_v"])
        require(values.shape == (count, 6), "Decision CSV must cover the complete payload")
        exact(values[:, 0], np.arange(count), "Decision CSV payload order")
        exact(values[:, 1], truth, "Original source-bit truth")
        exact(values[:, 2], payload_decisions, "Independent decision recurrence")
        raw_error = close(values[:, 3], raw_payload, "Every persisted post-channel sample")
        summer_error = close(values[:, 4], corrected, "Every persisted post-DFE sample")
        close(values[:, 5], np.full(count, threshold), "Every persisted threshold")
        metrics = metric_reference(corrected, truth, threshold, payload_decisions)
        verify_metrics(case["metrics"], metrics, f"{ntaps}-tap corrected metrics")
        raw_threshold = bounds(response, p, latency, 0)["threshold_v"]
        close(case["raw_channel_threshold_v"], raw_threshold, "Raw-channel threshold")
        raw_metrics = metric_reference(raw_payload, truth, raw_threshold, raw_payload >= raw_threshold)
        verify_metrics(case["raw_channel"], raw_metrics, f"{ntaps}-tap raw-channel metrics")
        searches = csv_values(folder, f"phase_search_{ntaps}tap.csv",
                              ["latency_samples", "phase_samples", "threshold_v", "worst_case_eye_v", "main_cursor_v"])
        expected_search = np.array([[candidate[key] for key in
                                    ["latency_samples", "phase_samples", "threshold_v", "worst_case_eye_v", "main_cursor_v"]]
                                   for candidate in candidates])
        close(searches, expected_search, "Every persisted phase-search candidate")
        sensitivity = summary["verification"]["truncation_20_to_40_ns"][str(ntaps)]
        longer = bounds(response40, p, latency, ntaps)
        close(sensitivity["fixed_latency_threshold_change_v"], longer["threshold_v"] - threshold, "40 ns threshold sensitivity")
        close(sensitivity["fixed_latency_eye_bound_change_v"], longer["worst_case_eye_v"] - expected["worst_case_eye_v"], "40 ns eye sensitivity")
        close(sensitivity["max_tap_change_v"], np.max(np.abs(longer["taps_v"] - taps)), "40 ns tap sensitivity")
        snapshot = snapshot_reference(folder, case["snapshot"], bits, response, p, dt,
                                      phase, taps, decisions, scored, truth, payload_decisions)
        require((folder / case["eye_file"]).stat().st_size > 1000, "Missing required eye figure")
        require(case["verification"]["histogram_samples"] == count * p, "Eye must include the complete payload UI count")
        report["cases"].append({"tap_count": ntaps, "decisions_verified": count,
                                "post_channel_samples_max_abs_v": raw_error,
                                "post_dfe_samples_max_abs_v": summer_error,
                                "errors_verified": metrics["errors"], "raw_errors_verified": raw_metrics["errors"],
                                "phase_candidates_verified": p, "snapshot": snapshot})
        references[ntaps] = (raw, summer, decisions, locations, scored, expected)
    raw, summer, decisions, locations, scored, design = references[3]
    taps, phase = design["taps_v"], design["phase_samples"]
    stress = summary["stress"]
    threshold = design["threshold_v"] + 0.6 * design["main_cursor_v"]
    close(stress["threshold_v"], threshold, "Stress threshold")
    stress_summer, stress_decisions = history_feedback(raw, taps, threshold)
    stress_metrics = metric_reference(stress_summer[locations], truth, threshold, stress_decisions[locations])
    verify_metrics(stress["metrics"], stress_metrics, "Stress metrics")
    error_index = stress_metrics["first_error_payload_bit"]
    example = stress["first_error"]
    if error_index is None:
        require(example is None, "Stress first-error example must be absent when no errors occur")
    else:
        require(example["payload_bit"] == error_index and example["true_bit"] == int(truth[error_index]) and
                example["decision"] == int(stress_decisions[locations[error_index]]), "Incorrect stress first-error bits")
        close(example["post_channel_v"], raw[locations[error_index]], "Stress first-error raw voltage")
        close(example["summer_v"], stress_summer[locations[error_index]], "Stress first-error summer voltage")
        close(example["threshold_v"], threshold, "Stress first-error threshold")
        start = max(0, error_index - 5)
        require(example["source_context_start"] == start, "Incorrect first-error source context offset")
        exact(example["source_context"], truth[start:error_index + 6], "Stress source context")
    report["stress"] = {"errors_verified": stress_metrics["errors"],
                         "snapshot": snapshot_reference(folder, stress["snapshot"], bits, response, p, dt,
                                                        phase, taps, stress_decisions, scored, truth,
                                                        stress_decisions[locations])}
    injection = summary["injected_error"]
    payload_index = int(injection["payload_bit"])
    force_at = int(locations[payload_index])
    changed, changed_decisions = history_feedback(raw, taps, design["threshold_v"], force_at)
    require(injection["forced_from"] == int(decisions[force_at]) and
            injection["forced_to"] == int(changed_decisions[force_at]), "Incorrect injected flip description")
    rows = csv_values(folder, injection["csv"],
                      ["payload_bit", "truth", "baseline_decision", "injected_run_decision",
                       "baseline_summer_v", "injected_summer_v", "voltage_difference_v"])
    row_indices = np.arange(payload_index, payload_index + 12)
    row_clocks = locations[row_indices]
    expected_rows = np.column_stack([row_indices, truth[row_indices], decisions[row_clocks],
                                    changed_decisions[row_clocks], summer[row_clocks], changed[row_clocks],
                                    changed[row_clocks] - summer[row_clocks]])
    injection_error = close(rows, expected_rows, "Every persisted injected-error diagnostic row")
    extra_changes = int(np.count_nonzero(changed_decisions[locations[payload_index + 1:]] != decisions[locations[payload_index + 1:]]))
    require(extra_changes == injection["extra_decision_changes_after_injection"], "Incorrect subsequent decision-change count")
    close(injection["next_three_voltage_changes_v"],
          changed[locations[payload_index + 1:payload_index + 4]] - summer[locations[payload_index + 1:payload_index + 4]],
          "Injected error's next three voltage effects")
    report["injected_error"] = {"rows_verified": len(rows), "maximum_abs_v": injection_error,
                                 "subsequent_decision_changes_verified": extra_changes}
    report["passed"] = True
    report["limitations"] = [
        "Finite noiseless model verification; observed error fractions are not physical BER estimates.",
        "Worst-case residual bounds assume correct previous decisions.",
        "Eye image existence and complete sample count are checked; raster pixels and every histogram bin are not independently regenerated.",
    ]
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--write", type=Path, help="JSON evidence path; defaults to results/verification.json")
    args = parser.parse_args(argv)
    report = verify(args.results_dir)
    destination = args.write if args.write else args.results_dir / "verification.json"
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
