"""Independently reproduce all six saved HW4 spectra and numerical checks.

Run from any directory:
    python homework/hw04/tests/verify_reference.py
    python homework/hw04/tests/verify_reference.py --results-dir PATH

The full-record reference is a causal boxcar convolution, followed by SciPy's
periodogram implementation. It never calls the rewritten waveform or spectral
estimators. The ordinary unit suite is also run once. Only after every check
passes is verification.json written to the selected results directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import scipy
from scipy.signal import periodogram
from scipy.signal.windows import kaiser

HOME = Path(__file__).resolve().parents[1]
EXPECTED_BIT_COUNT = 200_000
SPECTRAL_RELATIVE_TOLERANCE = 1e-12
POWER_RELATIVE_TOLERANCE = 1e-12


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_default_configuration(config: dict) -> None:
    """Reject partial or differently configured results rather than guessing."""
    expected = {
        "bit_duration": 100e-12,
        "rise_time": 20e-12,
        "dt": 2e-12,
        "vpp": 1.0,
        "z0": 50.0,
        "block_bits": 2_000,
        "nfft": 2**19,
        "kaiser_beta": float(np.pi),
        "averages": [10, 100],
    }
    for name, value in expected.items():
        require(config.get(name) == value,
                f"This independent audit requires the default {name}={value!r}.")


def load_reference_bits(path: Path) -> tuple[bytes, np.ndarray]:
    """Validate the unchanged source format before signed numerical conversion."""
    raw = path.read_bytes()
    require(len(raw) == EXPECTED_BIT_COUNT,
            f"{path.name} must contain exactly {EXPECTED_BIT_COUNT:,} binary bytes.")
    values = np.frombuffer(raw, dtype=np.uint8)
    require(bool(np.all((values == ord("0")) | (values == ord("1")))),
            f"{path.name} contains a nonbinary byte.")
    return raw, values.astype(np.int16) - ord("0")


def run_unit_tests() -> unittest.TestResult:
    suite = unittest.defaultTestLoader.discover(
        start_dir=str(HOME / "tests"), pattern="test_hw4.py"
    )
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    require(result.wasSuccessful(),
            "Unit tests must pass before accepting the generated spectra.")
    require(result.testsRun > 0, "The numerical unit suite was not discovered.")
    return result


def verify(results_dir: Path) -> dict:
    result = run_unit_tests()
    summary = json.loads((results_dir / "summary.json").read_text(encoding="utf-8"))
    config = summary["configuration"]
    validate_default_configuration(config)

    dt, z0, vpp = config["dt"], config["z0"], config["vpp"]
    rise_samples = round(config["rise_time"] / dt)
    samples_per_bit = round(config["bit_duration"] / dt)
    block_samples = samples_per_bit * config["block_bits"]
    window = kaiser(block_samples, config["kaiser_beta"], sym=True)
    window_energy = np.sum(window**2)

    audit = {
        "title": "Independent numerical verification of HW4 results",
        "status": "passed",
        "unit_test_count": result.testsRun,
        "unit_test_failures": len(result.failures),
        "unit_test_errors": len(result.errors),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "independent_method": (
            "Construct complete bipolar NRZ stream by causal boxcar convolution "
            "with preceding zero; split its samples into blocks; use "
            "scipy.signal.periodogram (detrend=False, density, symmetric Kaiser) "
            "and average power spectra. No rewritten estimator is called."
        ),
        "relative_error_policy": (
            "Peak pointwise error is computed only at bins above 1e-14 times "
            "the maximum reference spectral density. Relative L2 and maximum "
            "absolute errors include all bins."
        ),
        "tolerances": {
            "spectral_relative_L2": SPECTRAL_RELATIVE_TOLERANCE,
            "full_band_power_relative": POWER_RELATIVE_TOLERANCE,
        },
        "spectral_comparisons": [],
    }

    with np.load(results_dir / "spectra.npz", allow_pickle=False) as spectra:
        for number in (1, 2, 3):
            path = HOME / "instructions" / f"HW4_bit_pattern_{number}.txt"
            raw, bits = load_reference_bits(path)
            ideal = np.repeat((bits - 0.5) * vpp, samples_per_bit)
            # The first sample at a transition retains the preceding level.
            # Using an exclusive causal average gives j/r at sample j of a ramp.
            prehistory = np.full(rise_samples, -vpp / 2)
            waveform = np.convolve(
                np.r_[prehistory, ideal],
                np.ones(rise_samples) / rise_samples,
                mode="valid",
            )[:-1]
            require(waveform.size == bits.size * samples_per_bit,
                    "Independent convolution produced an incorrect record length.")
            rows = waveform.reshape(-1, block_samples)
            require(len(rows) == 100, "The audit must cover all 100 full segments.")

            total = np.zeros(config["nfft"] // 2 + 1)
            weighted_powers, raw_powers = [], []
            for k, row in enumerate(rows, 1):
                frequency, density = periodogram(
                    row, fs=1 / dt, window=window, nfft=config["nfft"],
                    detrend=False, return_onesided=True, scaling="density",
                )
                total += density / z0
                weighted_powers.append(
                    float(np.sum((row * window)**2) / window_energy / z0)
                )
                raw_powers.append(float(np.mean(row**2) / z0))
                if k not in config["averages"]:
                    continue

                expected = total / k
                saved = spectra[f"pattern_{number}_K{k}_w_per_hz"]
                require(saved.shape == frequency.shape,
                        f"Pattern {number}, K={k}: spectral shape mismatch.")
                require(bool(np.all(np.isfinite(saved)) and np.all(saved >= 0)),
                        f"Pattern {number}, K={k}: invalid linear power density.")
                np.testing.assert_allclose(
                    spectra["frequency_hz"], frequency, rtol=0, atol=0
                )

                difference = saved - expected
                significant = expected > expected.max() * 1e-14
                relative_l2 = float(np.linalg.norm(difference) / np.linalg.norm(expected))
                relative_peak = float(np.max(
                    np.abs(difference[significant]) / expected[significant]
                ))
                maximum_absolute = float(np.abs(difference).max())
                df = frequency[1] - frequency[0]
                power = float(saved.sum() * df)
                expected_power = float(np.mean(weighted_powers))
                power_error = abs(power - expected_power) / expected_power
                require(relative_l2 < SPECTRAL_RELATIVE_TOLERANCE,
                        f"Pattern {number}, K={k}: independent spectral mismatch.")
                require(power_error < POWER_RELATIVE_TOLERANCE,
                        f"Pattern {number}, K={k}: full-band power mismatch.")
                audit["spectral_comparisons"].append({
                    "pattern": number,
                    "K": k,
                    "source_sha256": hashlib.sha256(raw).hexdigest(),
                    "fft_bins_compared": len(frequency),
                    "maximum_absolute_error_w_per_hz": maximum_absolute,
                    "relative_L2_error": relative_l2,
                    "maximum_relative_error_above_threshold": relative_peak,
                    "saved_full_band_power_w": power,
                    "independent_window_weighted_power_w": expected_power,
                    "independent_unwindowed_power_w": float(np.mean(raw_powers)),
                    "relative_full_band_power_error": power_error,
                })
                print(f"Pattern {number}, K={k}: relative L2={relative_l2:.3e}; "
                      f"relative full-band power error={power_error:.3e}")

        # Exact ensemble power for the discrete 50-sample hold and 10-sample ramp.
        frequency = spectra["frequency_hz"]
        df = frequency[1] - frequency[0]
        expected_power = vpp**2 / (4 * z0) * (
            1 - (rise_samples - 1 / rise_samples) / (3 * samples_per_bit)
        )
        sampled_reference = spectra["sampled_random_reference_w_per_hz"]
        require(sampled_reference.shape == frequency.shape,
                "The sampled random reference has an invalid shape.")
        require(bool(np.all(np.isfinite(sampled_reference))
                     and np.all(sampled_reference >= 0)),
                "The sampled random reference must be finite linear power density.")
        integrated = float(sampled_reference.sum() * df)
        reference_error = abs(integrated - expected_power) / expected_power
        audit["sampled_random_reference_check"] = {
            "exact_expected_power_w": expected_power,
            "integrated_reference_power_w": integrated,
            "relative_error": reference_error,
        }
        require(reference_error < POWER_RELATIVE_TOLERANCE,
                "The sampled random reference fails its exact total-power check.")

    require(len(audit["spectral_comparisons"]) == 6,
            "All six independently reconstructed spectra must be checked.")
    output = results_dir / "verification.json"
    output.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(f"Independent verification passed; wrote {output.name}.")
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir", type=Path, default=HOME / "results",
        help="Results containing summary.json and spectra.npz (default: package results/).",
    )
    args = parser.parse_args()
    verify(args.results_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
