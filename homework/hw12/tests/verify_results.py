#!/usr/bin/env python3
"""Independently verify saved HW12 numerical artifacts without rerunning them.

Default: require the assignment's >=2**20 evaluated bits. Use --allow-short
only for a smoke run, which is explicitly marked insufficient for submission.
This verifier does not import the solver or trust its self-reported checks.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CASE_KEYS = ("case1_no_tx_jitter", "case2_sj", "case3_rj", "case4_sj_rj")
MINIMUM_BITS = 2**20


def require(condition, description):
    if not bool(condition):
        raise AssertionError(description)


def close(actual, expected, description, rtol=2e-11, atol=0.0):
    require(np.allclose(actual, expected, rtol=rtol, atol=atol), description)


def read_distributions(path):
    groups = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["case", "side", "time_seconds", "pmf", "pdf_per_second"],
                f"Unexpected distribution CSV columns: {path.name}")
        for row in reader:
            key = (row["case"], row["side"])
            groups.setdefault(key, []).append([float(row[name]) for name in reader.fieldnames[2:]])
    require(set(groups) == {(case, side) for case in CASE_KEYS for side in ("tx", "rx")},
            f"Missing, duplicate-name, or unexpected cases: {path.name}")
    return {key: np.asarray(rows) for key, rows in groups.items()}


def verify_distribution(table, raw, stats, points, dt, description):
    require(table.shape == (points, 3), f"{description}: must contain exactly one UI")
    require(np.all(np.isfinite(table)), f"{description}: nonfinite CSV values")
    require(np.all(table[:, 1:] >= 0), f"{description}: negative probability")
    time_axis = (np.arange(points) - points // 2) * dt
    close(table[:, 0], time_axis, f"{description}: incorrect or endpoint-inclusive time axis", atol=1e-26)
    close(np.sum(table[:, 1]), 1, f"{description}: PMF not normalized", atol=1e-13)
    close(table[:, 2], table[:, 1] / dt, f"{description}: PDF != PMF / dt", atol=1e-15)
    close(np.sum(table[:, 2]) * dt, 1, f"{description}: PDF integral not one", atol=1e-13)
    raw = np.asarray(raw)
    require(raw.shape == (points,), f"{description}: raw histogram shape")
    require(np.all(raw >= 0) and np.all(np.isfinite(raw)), f"{description}: invalid raw counts")
    total = float(np.sum(raw))
    require(total > 0, f"{description}: no observations")
    reference = stats["reference_index"]
    require(isinstance(reference, int) and 0 <= reference < points, f"{description}: invalid reference")
    expected = np.roll(raw.astype(float), points // 2 - reference) / total
    close(table[:, 1], expected, f"{description}: CSV differs from raw counts or wrap reference", atol=1e-15)
    require(np.count_nonzero(table[:, 1]) == np.count_nonzero(raw), f"{description}: zero cells changed")
    mean = float(np.dot(time_axis, expected))
    variance = float(np.dot((time_axis - mean)**2, expected))
    close(stats["observation_count"], total, f"{description}: observation count")
    close(stats["mean_seconds"], mean, f"{description}: incorrect mean", atol=1e-26)
    close(stats["variance_seconds2"], variance, f"{description}: incorrect centered variance", atol=1e-38)
    close(stats["rms_seconds"], math.sqrt(variance), f"{description}: RMS != standard deviation", atol=1e-26)


def verify(results=ROOT / "results", allow_short=False):
    results = Path(results)
    summary = json.loads((results / "summary.json").read_text(encoding="utf-8"))
    params, cases, channel = summary["parameters"], summary["cases"], summary["channel"]
    bits, points = params["bits"], params["points_per_ui"]
    require(isinstance(bits, int) and bits > 0, "Evaluated bit count must be a positive integer")
    require(isinstance(points, int) and points > 0, "Invalid samples/UI")
    sufficient = bits >= MINIMUM_BITS
    require(sufficient or allow_short, "Assignment requires at least 2**20 evaluated bits; this is a smoke run")
    require(set(cases) == set(CASE_KEYS), "Four required cases must be present")
    for key, expected in {"ui_seconds": 25e-12, "rise_seconds": 5e-12,
                          "vpp_volts": 1.0, "sj_peak_seconds": 3e-12,
                          "sj_frequency_hz": 1e6, "rj_sigma_seconds": 1e-12,
                          "rs_ohms": 20.0, "rl_ohms": 85.0,
                          "cs_farads": 0.3e-12, "cl_farads": 0.3e-12,
                          "dv_volts": 1e-3, "kaiser_a": 0.5,
                          "truncation_seconds": 40e-9}.items():
        close(params[key], expected, f"Assignment parameter mismatch: {key}")
    dt = params["dt_seconds"]
    close(dt * points, params["ui_seconds"], "dt * samples/UI != one UI")
    require(params["warmup_bits"] * points >= channel["fir_taps"] - 1,
            "Warmup shorter than the retained FIR history")
    require(params["block_bits"] > 0, "Invalid block length")
    require(params["requested_bits"] == bits and params["analyzed_bits"] == bits,
            "Requested and analyzed bit counts differ")
    require(params["guard_bits_each_side"] >= 2, "Missing waveform guard bits")
    close(params["evaluated_sj_cycles"], bits * params["ui_seconds"] * params["sj_frequency_hz"],
          "Incorrect evaluated SJ cycle count")
    close(params["threshold_tx_volts"], params["vpp_volts"] / 2, "Tx threshold")
    close(params["threshold_rx_volts"], params["vpp_volts"] * channel["dc_gain"] / 2, "Rx threshold")
    require(abs(channel["fir_dc_relative_error"]) < 0.01, "FIR truncation loses more than 1% DC gain")
    eye_csv = read_distributions(results / "jitter_distributions.csv")
    crossing_csv = read_distributions(results / "crossing_distributions.csv")
    references = {side: cases[CASE_KEYS[0]][side]["reference_index"] for side in ("tx", "rx")}
    crossing_references = {side: cases[CASE_KEYS[0]]["crossing_" + side]["reference_index"] for side in ("tx", "rx")}
    # A full, non-streamed regeneration is an independent oracle for section
    # boundaries, warmup and RNG consumption; it never generates waveforms.
    total = params["warmup_bits"] + bits
    bit_rng = np.random.default_rng(np.random.SeedSequence([params["seed"], 11]))
    normal_rng = np.random.default_rng(np.random.SeedSequence([params["seed"], 29]))
    pattern = bit_rng.integers(0, 2, total + 4, dtype=np.int64)
    normal = normal_rng.normal(size=total + 4)
    indices = np.arange(-2, total + 2)
    changes = np.r_[False, pattern[1:] != pattern[:-1]]
    selected = changes & (indices >= params["warmup_bits"]) & (indices < total)
    edge_indices, edge_normals = indices[selected], normal[selected]
    per_case = {}
    for case in CASE_KEYS:
        with np.load(results / f"{case}_eyes.npz", allow_pickle=False) as saved:
            require(int(saved["evaluated_bits"]) == bits, f"{case}: saved bit count")
            require(int(saved["points_per_ui"]) == points, f"{case}: saved samples/UI")
            close(float(saved["dt"]), dt, f"{case}: saved dt")
            close(float(saved["dv"]), params["dv_volts"], f"{case}: saved dv")
            for side in ("tx", "rx"):
                eye = saved[f"{side}_counts"]
                require(eye.ndim == 2 and eye.shape[1] == points, f"{case}/{side}: eye shape")
                require(np.issubdtype(eye.dtype, np.integer), f"{case}/{side}: eye contains pseudocounts")
                require(np.all(eye >= 0), f"{case}/{side}: negative counts")
                require(np.all(eye.sum(axis=0) == bits), f"{case}/{side}: exact bit count / phase accounting")
                stats = cases[case][side]
                vmin = float(saved[f"{side}_vmin"])
                row = math.floor((params[f"threshold_{side}_volts"] - vmin) / params["dv_volts"])
                require(0 <= row < eye.shape[0], f"{case}/{side}: threshold outside eye")
                require(stats["voltage_bin_index"] == row, f"{case}/{side}: incorrect threshold row")
                close(stats["voltage_bin_edges"], [vmin + row * params["dv_volts"], vmin + (row + 1) * params["dv_volts"]],
                      f"{case}/{side}: threshold bin edges")
                require(stats["reference_index"] == references[side], f"{case}/{side}: reference changed across cases")
                verify_distribution(eye_csv[(case, side)], eye[row], stats, points, dt, f"{case}/{side}/eye")
                crossing = cases[case]["crossing_" + side]
                require(crossing["reference_index"] == crossing_references[side], f"{case}/{side}: crossing reference changed")
                verify_distribution(crossing_csv[(case, side)], saved[f"{side}_crossing_counts"], crossing,
                                    points, dt, f"{case}/{side}/crossings")
        per_case[case] = {f"{side}_sigma_ps": cases[case][side]["rms_seconds"] * 1e12 for side in ("tx", "rx")}
        sj = params["sj_peak_seconds"] if case in ("case2_sj", "case4_sj_rj") else 0
        rj = params["rj_sigma_seconds"] if case in ("case3_rj", "case4_sj_rj") else 0
        requested = sj * np.sin(2 * np.pi * params["sj_frequency_hz"] * edge_indices * params["ui_seconds"]) + rj * edge_normals
        for label, samples in (("requested", requested), ("applied", np.rint(requested / dt) * dt)):
            stats = cases[case]["input_jitter"][label]
            require(stats["observation_count"] == len(samples), f"{case}/{label}: RNG, warmup or true-edge count")
            close(stats["mean_seconds"], float(np.mean(samples)), f"{case}/{label}: regenerated input mean", atol=1e-26)
            close(stats["variance_seconds2"], float(np.var(samples)), f"{case}/{label}: regenerated input variance", atol=1e-38)
            close(stats["rms_seconds"], float(np.std(samples)), f"{case}/{label}: regenerated input RMS", atol=1e-26)
    require(cases[CASE_KEYS[0]]["tx"]["variance_seconds2"] == 0, "Ideal Tx must have exactly zero variance")
    theory_checks = {}
    if sufficient:
        for case, expected_ps2 in zip(CASE_KEYS[1:], [4.5, 1.0, 5.5]):
            measured = cases[case]["tx"]["variance_seconds2"] * 1e24
            # Finite SJ coverage, time quantization and finite counts justify
            # a tolerance; Rx is deliberately not forced to additive theory.
            relative_error = (measured - expected_ps2) / expected_ps2
            require(abs(relative_error) < 0.05, f"{case}: Tx variance differs from theory by >=5%")
            theory_checks[case] = {"measured_ps2": measured, "theoretical_ps2": expected_ps2,
                                   "relative_error": relative_error}
        tx2, tx3, tx4 = [cases[case]["tx"]["variance_seconds2"] for case in CASE_KEYS[1:]]
        independence_error = (tx4 - tx2 - tx3) / (tx2 + tx3)
        require(abs(independence_error) < 0.05, "Independent Tx SJ/RJ variance addition differs by >=5%")
        theory_checks["sj_rj_variance_addition_relative_error"] = independence_error
    return {"successful": True, "assignment_sufficient": sufficient, "evaluated_bits": bits,
            "points_per_ui": points, "distributions_checked": 16,
            "per_case": per_case, "tx_theoretical_checks": theory_checks,
            "scope": "Independent saved-artifact consistency, normalization, bit accounting, and Tx physical checks"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--allow-short", action="store_true")
    parser.add_argument("--output", type=Path, help="Optional JSON verification record")
    args = parser.parse_args(argv)
    try:
        record = verify(args.results, args.allow_short)
    except (AssertionError, KeyError, OSError, ValueError) as error:
        record = {"successful": False, "error": str(error)}
    text = json.dumps(record, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if record["successful"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
