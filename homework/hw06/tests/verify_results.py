"""Independently verify saved HW6 results without importing the implementation.

This uses an explicit Hermitian spectrum, polyphase cursor sums, and
bit-rate convolution, different algorithms from src/hw6.py. Run from any
directory: ``python homework/hw06/tests/verify_results.py``. NumPy only.
"""

import hashlib
import json
import math
from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parents[1]


def independent_boundaries(pulse, n, tolerance):
    peak = int(pulse.argmax())
    start = peak - n
    main = pulse[start:start + 2 * n]
    phase = (start + np.arange(2 * n)) % n
    positive = np.array([np.maximum(pulse[r::n], 0).sum() for r in range(n)])[phase]
    negative = np.array([np.minimum(pulse[r::n], 0).sum() for r in range(n)])[phase]
    kept_positive = np.array([np.where(pulse[r::n] >= tolerance, pulse[r::n], 0).sum()
                              for r in range(n)])[phase]
    kept_negative = np.array([np.where(pulse[r::n] <= -tolerance, pulse[r::n], 0).sum()
                              for r in range(n)])[phase]
    # Sum each phase including the main cursor, then remove its contribution.
    zero_exact = positive - np.maximum(main, 0)
    one_exact = main + negative - np.minimum(main, 0)
    zero_filtered = kept_positive - np.where(main >= tolerance, main, 0)
    one_filtered = main + kept_negative - np.where(main <= -tolerance, main, 0)
    return {"main": main, "wc0_exact": zero_exact, "wc1_exact": one_exact,
            "wc0": zero_filtered, "wc1": one_filtered,
            "omitted_positive": zero_exact - zero_filtered,
            "omitted_negative": one_filtered - one_exact,
            "outer_min": negative, "outer_max": positive}, start


def verify_saved_results():
    summary = json.loads((BASE / "results" / "summary.json").read_text(encoding="utf-8"))
    manifest = json.loads((BASE / "source_manifest.json").read_text(encoding="utf-8"))
    for item in manifest["provided_instruction_files"]:
        data = (BASE / "instructions" / item["file"]).read_bytes()
        assert len(data) == item["bytes"]
        assert hashlib.sha256(data).hexdigest() == item["sha256"], item["file"]
    arrays = np.load(BASE / "results" / "arrays.npz", allow_pickle=False)
    cfg = summary["config"]
    n, dt = cfg["samples_per_ui"], float(arrays["dt_s"])
    table = np.loadtxt(BASE / "instructions" / "HW6_channel.s2p", comments=("!", "#"))
    f = table[:, 0]
    s21 = table[:, 3] * np.exp(1j * np.deg2rad(table[:, 4]))
    np.testing.assert_array_equal(arrays["frequency_hz"], f)
    # libm complex-exponential results can differ by a few ulps across systems.
    # The original channel bytes are checked exactly by SHA-256 above.
    np.testing.assert_allclose(arrays["s21"], s21, rtol=0,
                               atol=8 * np.finfo(np.float64).eps)
    df, m = f[1] - f[0], len(arrays["full_sbr"])
    tb, tr = cfg["bit_period_ps"] * 1e-12, cfg["rise_time_ps"] * 1e-12
    window = np.i0(np.pi * cfg["kaiser_alpha"] * np.sqrt(1 - (f / f[-1]) ** 2))
    window /= np.i0(np.pi * cfg["kaiser_alpha"])
    transform = (cfg["vpp_v"] * tb * np.sinc(f * tb) * np.sinc(f * tr)
                 * np.exp(-1j * np.pi * f * (tb + tr)) * s21 * window / 2)
    spectrum = np.zeros(m, dtype=complex)
    spectrum[:len(f)] = transform
    spectrum[-len(f) + 1:] = transform[:0:-1].conj()
    pulse = np.fft.ifft(spectrum).real * (m * df)
    np.testing.assert_allclose(arrays["full_sbr"], pulse, rtol=2e-12, atol=2e-15)
    sample_indices = np.array([0, 7990, 8011, 8020, 10000, 30000])
    direct = np.array([df * (transform[0].real + 2 * np.real(np.sum(
        transform[1:] * np.exp(2j * np.pi * f[1:] * k * dt)))) for k in sample_indices])
    np.testing.assert_allclose(pulse[sample_indices], direct, rtol=2e-11, atol=2e-13)
    pulse = pulse[:cfg["sbr_samples"]]
    np.testing.assert_allclose(arrays["sbr"], pulse, atol=2e-15)
    expected, start = independent_boundaries(pulse, n, cfg["tolerance_v"])
    for key, value in expected.items():
        np.testing.assert_allclose(arrays[f"pda_{key}"], value, atol=3e-15, err_msg=key)
    assert start == summary["sbr"]["pda_start_index"]
    np.testing.assert_array_equal(arrays["pda_time_ui"], np.arange(2 * n) / n)
    assert int(pulse.argmax()) == summary["sbr"]["peak_index"]
    # Check every saved truncation study independently, including raw bounds.
    convergence = []
    for entry in summary["truncation_study"]:
        cut = round(entry["truncation_ns"] * 1e-9 / dt) + 1
        b, _ = independent_boundaries(arrays["full_sbr"][:cut], n, cfg["tolerance_v"])
        exact_gap = float(np.max(b["wc1_exact"] - b["wc0_exact"]))
        filtered_gap = float(np.max(b["wc1"] - b["wc0"]))
        np.testing.assert_allclose(exact_gap, entry["exact"]["maximum_vertical_gap_v"], atol=3e-15)
        np.testing.assert_allclose(filtered_gap, entry["filtered"]["maximum_vertical_gap_v"], atol=3e-15)
        convergence.append({"truncation_ns": entry["truncation_ns"],
                            "exact_gap_v": exact_gap, "filtered_gap_v": filtered_gap})
    binary = "".join((BASE / "instructions" / "HW6_wcbitpattern.txt").read_text().split())
    np.testing.assert_array_equal(arrays["pattern_bits"], np.array(list(binary), dtype=np.uint8))
    np.testing.assert_array_equal(arrays["random_bits"], np.random.default_rng(cfg["random_seed"]).integers(
        0, 2, len(binary), dtype=np.uint8))
    report = {"status": "passed", "instruction_hashes_verified": len(manifest["provided_instruction_files"]),
              "full_sbr_hermitian_ifft_verified": True, "direct_frequency_sum_samples_verified": len(sample_indices),
              "polyphase_pda_phases_verified": 2 * n, "truncation_study": convergence, "cases": {}}
    for case in ("random", "pattern"):
        print(f"Independently checking every {case} eye sample...", flush=True)
        bits = arrays[f"{case}_bits"]
        guard = math.ceil(len(pulse) / n) + 2
        first, stop = guard, len(bits) - guard
        labels = bits[first:stop].astype(bool)
        size = 1 << (len(bits) + math.ceil(len(pulse) / n) - 2).bit_length()
        bits_fft = np.fft.rfft(bits.astype(float), n=size)
        ones, zeros = np.empty(2 * n), np.empty(2 * n)
        histogram = np.empty_like(arrays[f"{case}_histogram"])
        for q in range(n):
            phase, offset = (start + q) % n, (start + q) // n
            taps = pulse[phase::n]
            # Each waveform phase is an ordinary symbol-rate FIR convolution.
            output_phase = np.fft.irfft(bits_fft * np.fft.rfft(taps, n=size), n=size)
            a = output_phase[first + offset:stop + offset]
            b = output_phase[first + offset + 1:stop + offset + 1]
            ones[q], ones[q + n] = a[labels].min(), b[labels].min()
            zeros[q], zeros[q + n] = a[~labels].max(), b[~labels].max()
            histogram[:, q] = np.histogram(a, bins=arrays["voltage_edges_v"])[0]
        np.testing.assert_allclose(ones, arrays[f"{case}_observed_min_one"], atol=3e-14)
        np.testing.assert_allclose(zeros, arrays[f"{case}_observed_max_zero"], atol=3e-14)
        np.testing.assert_array_equal(histogram, arrays[f"{case}_histogram"])
        np.testing.assert_array_equal(histogram.sum(axis=0), stop - first)
        probability = histogram / (stop - first)
        np.testing.assert_allclose(probability, arrays[f"{case}_probability"], atol=1e-15)
        np.testing.assert_allclose(probability.sum(axis=0), 1, atol=3e-14)
        assert np.min(ones - expected["wc1_exact"]) >= -3e-14
        assert np.min(expected["wc0_exact"] - zeros) >= -3e-14
        # A filtered bound permits the sum of every omitted tap, rather than
        # only one tolerance value, as the error margin at each phase.
        assert np.min(ones - (expected["wc1"] - expected["omitted_negative"])) >= -3e-14
        assert np.min(expected["wc0"] + expected["omitted_positive"] - zeros) >= -3e-14
        report["cases"][case] = {"input_bits_verified": len(bits), "stable_windows_verified": stop - first,
                                 "waveform_samples_verified": int((stop - first) * 2 * n),
                                 "histogram_count_conserved": True, "raw_and_tolerance_bounds_verified": True,
                                 "maximum_extrema_difference_v": float(max(
                                     np.max(abs(ones - arrays[f"{case}_observed_min_one"])),
                                     np.max(abs(zeros - arrays[f"{case}_observed_max_zero"]))))}
    (BASE / "results" / "verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    result = verify_saved_results()
    print(json.dumps(result, indent=2))
