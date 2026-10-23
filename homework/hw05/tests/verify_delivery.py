"""Independently verify the default HW5 numerical outputs and anonymous PDF.

Run from any working directory. Requires NumPy and pypdf. This verifier uses
integer sample-grid coordinates, per-cycle local elapsed time, and histogram2d;
the solution uses global crossing times, searchsorted, and flattened bincount.
The reference is intentionally limited to the supplied 1 ps data and default
25 ps UI, 1 ps time bins, -0.1..1.1 V range, and 0.01 V voltage bins.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

import numpy as np
from pypdf import PdfReader

EXPECTED_INPUT_SHA256 = "5e72b7e14bab20ee171442ae58f01c3594c656d71baf48926378c6d4c297b615"


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_source_coverage(base):
    manifest = json.loads((base / "source_manifest.json").read_text(encoding="utf-8"))
    expected = {"instructions/HW5_ADS_waveforms.csv",
                "instructions/HW5_stat_eye_construction.pdf",
                "instructions/HW5_stateye_construction_sample.py"}
    entries = manifest["instruction_files"]
    require({entry["path"] for entry in entries} == expected and len(entries) == 3,
            "Manifest must account for all three supplied instruction files")
    for entry in entries:
        path = base / entry["path"]
        require(path.is_file() and path.stat().st_size == entry["bytes"]
                and sha256(path) == entry["sha256"],
                "Supplied instruction bytes do not match the source manifest")
    rewrites = manifest["submission_rewrites"]
    require(len(rewrites) == 3, "Manifest must account for all three supplied submission files")
    for entry in rewrites:
        require(all((base / path).is_file() for path in entry["replacements"]),
                "Missing anonymous submission replacement")
    return {"instruction_files_verified": len(entries), "anonymous_source_replacements_verified": len(rewrites)}


def independent_reference(data):
    """Return default-grid counts, crossings, and exclusion statistics."""
    require(data.shape == (1_639_001, 4), "Unexpected assignment data shape")
    voltage = data[:, 1]
    grid = np.arange(len(data), dtype=np.int64) + 1000
    require(np.allclose(data[:, 0] * 1e12, grid, atol=3e-10, rtol=0),
            "Reference requires the supplied integer 1 ps sample grid")
    require(np.all(np.isfinite(data)), "Nonfinite input data")
    require(voltage.min() >= -0.1 and voltage.max() <= 1.1, "Voltage outside reference grid")
    voltage_edges = np.arange(-10, 111, dtype=float) / 100
    time_edges = np.arange(51, dtype=float)
    fixed = np.histogram2d(voltage, grid % 50, bins=(voltage_edges, time_edges))[0].astype(np.int64)

    clock = data[:, 2] - data[:, 3]
    right = np.where((clock[:-1] < 0) & (clock[1:] >= 0))[0] + 1
    fraction = -clock[right - 1] / (clock[right] - clock[right - 1])
    crossings = grid[right - 1] + fraction
    # Fill each disjoint cycle from the first sample after its crossing to the
    # sample immediately before the next. Local subtraction avoids large times.
    age = np.full(len(data), np.nan)
    for j, point in enumerate(right):
        stop = int(right[j + 1]) if j + 1 < len(right) else len(data)
        age[point:stop] = np.arange(stop - point, dtype=float) + (1 - fraction[j])
    valid = np.isfinite(age) & (age >= 0) & (age < 50)
    triggered = np.histogram2d(voltage[valid], age[valid], bins=(voltage_edges, time_edges))[0].astype(np.int64)
    statistics = {
        "input_samples": len(data),
        "positive_crossings": len(right),
        "fixed_total": int(fixed.sum()),
        "triggered_total": int(triggered.sum()),
        "samples_before_first_trigger": int(right[0]),
        "samples_excluded_after_2UI": int(np.sum(np.isfinite(age) & (age >= 50))),
        "fixed_column_count_range": [int(fixed.sum(axis=0).min()), int(fixed.sum(axis=0).max())],
        "triggered_column_count_range": [int(triggered.sum(axis=0).min()), int(triggered.sum(axis=0).max())],
    }
    return fixed, triggered, crossings, statistics


def verify_synthetic(source):
    """Compare the solution with explicit, hand-counted small examples."""
    spec = importlib.util.spec_from_file_location("hw5_independent_check_target", source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    config = module.EyeConfig(ui_ps=2, bins_per_ui=2, vmin=-0.1, vmax=0.1, dv=0.1)
    time = np.arange(13, dtype=float) * 1e-12
    voltage = np.array([-0.1, -0.05, 0, 0.05, 0.1, 0.05, 0, -0.05, -0.1, 0, 0.1, 0.05, 0])
    clock = np.array([-1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1], dtype=float)
    arrays, summary = module.analyze(np.column_stack((time, voltage, clock, np.zeros(13))), config)
    fixed_expected = np.array([[2, 1, 0, 1], [2, 2, 3, 2]])
    triggered_expected = np.array([[1, 1, 0, 0], [1, 1, 2, 2]])
    require(np.array_equal(arrays["fixed_counts"], fixed_expected), "Hand-counted fixed histogram mismatch")
    require(np.array_equal(arrays["triggered_counts"], triggered_expected), "Hand-counted triggered histogram mismatch")
    require(np.allclose(arrays["trigger_times_ps"], [1.5, 6.5], atol=1e-13, rtol=0), "Synthetic crossings mismatch")
    require(summary["trigger"]["samples_before_first_trigger"] == 2, "Synthetic pre-trigger exclusion mismatch")
    require(summary["trigger"]["samples_excluded_after_2UI"] == 3, "Synthetic late exclusion mismatch")
    test_times = np.array([0, 1, 1.5, 2, 5.5, 6.5, 7, 10.5])
    _, valid = module.triggered_time_bin_indices(test_times, np.array([1.5, 6.5]), config)
    require(np.array_equal(valid, [False, False, True, True, False, True, True, False]),
            "Exact 2UI endpoints must be excluded, never wrapped")
    probability, totals = module.vertical_normalize(np.array([[0, 1], [0, 3]]))
    require(np.array_equal(totals, [0, 4]) and np.array_equal(probability, [[0, 0.25], [0, 0.75]]),
            "Empty columns or hand-computed probabilities mismatch")
    return {"hand_counted_histograms": True, "exact_2UI_exclusion": True, "empty_column_normalization": True}


def verify_report(path):
    reader = PdfReader(path)
    require(len(reader.pages) > 0, "Report is empty")
    metadata = reader.metadata or {}
    author = str(metadata.get("/Author", "")).strip().lower()
    require(author in {"", "anonymous", "anonymous student"}, "Report author metadata is not anonymous")
    xmp = reader.xmp_metadata
    if xmp is not None:
        for creator in xmp.dc_creator or []:
            require(str(creator).strip().lower() in {"", "anonymous", "anonymous student"},
                    "Report XMP creator is not anonymous")
    text = "\n".join(page.extract_text() or "" for page in reader.pages)
    require(not re.search(r"(?<!\d)\d{9}(?!\d)", text), "Report contains a possible nine-digit student identifier")
    require(not re.search(r"[A-Za-z]:[\\/]", text), "Report contains an absolute Windows path")
    root = reader.trailer["/Root"]
    names = root.get("/Names")
    if names is not None:
        require("/EmbeddedFiles" not in names.get_object(), "Report contains embedded files")
    for page in reader.pages:
        for annotation in page.get("/Annots", []):
            require("/T" not in annotation.get_object(), "Report contains annotation author metadata")
    return {"pages": len(reader.pages), "neutral_or_blank_author": True, "no_embedded_files": True,
            "no_annotation_authors": True, "generic_identity_pattern_checks": True,
            "scope": "Metadata and generic patterns; private source-token scan is separate."}


def main(argv=None):
    base = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=base / "instructions" / "HW5_ADS_waveforms.csv")
    parser.add_argument("--results", type=Path, default=base / "results")
    parser.add_argument("--source", type=Path, default=base / "src" / "hw5.py")
    parser.add_argument("--report", type=Path, default=base / "report" / "hw05_report.pdf")
    args = parser.parse_args(argv)
    require(sha256(args.data) == EXPECTED_INPUT_SHA256, "Assignment CSV checksum mismatch")
    data = np.loadtxt(args.data, delimiter=",", skiprows=1)
    fixed, triggered, crossings, statistics = independent_reference(data)
    npz_path = args.results / "histograms.npz"
    summary_path = args.results / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    configuration = summary["configuration"]
    for key, expected in {"ui_ps": 25, "bins_per_ui": 25, "vmin": -0.1, "vmax": 1.1,
                          "dv": 0.01, "origin_ps": 0, "dt_ps": 1, "window_ps": 50,
                          "time_bins": 50, "voltage_bins": 120}.items():
        require(configuration[key] == expected, f"Reference requires default configuration: {key}")
    require(summary["input"]["sha256"] == EXPECTED_INPUT_SHA256, "Summary input checksum mismatch")
    matrices = {}
    with np.load(npz_path, allow_pickle=False) as arrays:
        require(np.allclose(arrays["time_edges_ps"], np.arange(51), atol=1e-12, rtol=0), "Time edges mismatch")
        require(np.allclose(arrays["voltage_edges_V"], np.arange(-10, 111) / 100, atol=1e-14, rtol=0), "Voltage edges mismatch")
        require(arrays["trigger_times_ps"].shape == crossings.shape, "Crossing count mismatch")
        crossing_error = float(np.max(np.abs(arrays["trigger_times_ps"] - crossings)))
        require(crossing_error < 3e-10, "Interpolated crossing times mismatch")
        for prefix, reference in [("fixed", fixed), ("triggered", triggered)]:
            actual = arrays[prefix + "_counts"]
            probability = arrays[prefix + "_probability"]
            totals = reference.sum(axis=0)
            require(actual.dtype.kind in "iu", f"{prefix} counts must be integers")
            require(np.array_equal(actual, reference), f"{prefix} full histogram mismatch")
            require(np.array_equal(arrays[prefix + "_column_counts"], totals), f"{prefix} column count mismatch")
            expected_probability = np.stack([reference[:, j] / totals[j] if totals[j] else np.zeros(120)
                                             for j in range(50)], axis=1)
            require(probability.shape == (120, 50) and np.all(np.isfinite(probability)), f"{prefix} invalid probabilities")
            require(np.allclose(probability, expected_probability, atol=1e-15, rtol=0), f"{prefix} normalization mismatch")
            error = float(np.max(np.abs(probability.sum(axis=0) - 1)))
            require(error < 1e-14, f"{prefix} occupied columns must sum to one")
            summary_key = "fixed_rate" if prefix == "fixed" else "triggered"
            require(summary[summary_key]["total_samples"] == int(reference.sum()), f"{prefix} summary total mismatch")
            require(summary[summary_key]["minimum_column_count"] == int(totals.min()), f"{prefix} summary minimum count mismatch")
            require(summary[summary_key]["maximum_column_count"] == int(totals.max()), f"{prefix} summary maximum count mismatch")
            matrices[prefix] = {"bins_compared": int(reference.size), "differing_bins": 0,
                                "total_samples": int(reference.sum()), "maximum_probability_error": float(np.max(np.abs(probability - expected_probability))),
                                "maximum_column_sum_error": error,
                                "count_array_sha256": hashlib.sha256(reference.astype("<i8").tobytes()).hexdigest()}
    require(summary["trigger"]["crossings"] == statistics["positive_crossings"], "Summary crossings mismatch")
    for key in ["samples_before_first_trigger", "samples_excluded_after_2UI"]:
        require(summary["trigger"][key] == statistics[key], f"Summary {key} mismatch")
    intervals = np.diff(crossings)
    phase = crossings - (crossings[0] + np.arange(len(crossings)) * 50)
    timing_metrics = {
        "first_crossing_ps": float(crossings[0]), "last_crossing_ps": float(crossings[-1]),
        "minimum_interval_ps": float(intervals.min()), "mean_interval_ps": float(intervals.mean()),
        "maximum_interval_ps": float(intervals.max()), "minimum_phase_error_ps": float(phase.min()),
        "maximum_phase_error_ps": float(phase.max()), "peak_to_peak_phase_error_ps": float(np.ptp(phase)),
    }
    for key, expected in timing_metrics.items():
        require(abs(summary["trigger"][key] - expected) < 1e-8, f"Summary timing metric mismatch: {key}")
    result = {
        "status": "passed",
        "source_coverage": verify_source_coverage(base),
        "independent_method": "Integer 1 ps grid; per-cycle fractional ages; numpy.histogram2d; hand-counted synthetic examples",
        "statistics": statistics,
        "independent_received_clock_timing_ps": timing_metrics,
        "matrices": matrices,
        "maximum_trigger_time_error_ps": crossing_error,
        "synthetic_checks": verify_synthetic(args.source),
        "report_checks": verify_report(args.report),
        "sha256": {"input_csv": sha256(args.data), "solution_source": sha256(args.source),
                   "histograms_npz": sha256(npz_path), "summary_json": sha256(summary_path),
                   "report_pdf": sha256(args.report)},
    }
    output = args.results / "verification.json"
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "passed", "full_histogram_bins_compared": 12_000,
                      "maximum_trigger_time_error_ps": crossing_error,
                      "verification": output.name}, indent=2))


if __name__ == "__main__":
    main()
