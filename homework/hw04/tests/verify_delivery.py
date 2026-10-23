"""Check published source coverage, saved spectra, and anonymous PDF metadata."""
from pathlib import Path
import hashlib
import json

import numpy as np
from pypdf import PdfReader

HOME = Path(__file__).resolve().parents[1]


def verify():
    manifest = json.loads((HOME / "source_manifest.json").read_text(encoding="utf-8"))
    for item in manifest["instruction_files"]:
        data = (HOME / item["path"]).read_bytes()
        assert len(data) == item["bytes"], item["path"]
        assert hashlib.sha256(data).hexdigest() == item["sha256"], item["path"]
    assert len(manifest["instruction_files"]) == 5

    summary = json.loads((HOME / "results/summary.json").read_text(encoding="utf-8"))
    config, derived = summary["configuration"], summary["derived_parameters"]
    assert config["averages"] == [10, 100]
    assert derived["unwindowed_resolution_scale_hz"] < 10e6
    assert derived["fft_grid_step_hz"] <= 1e6
    assert derived["window_enbw_hz"] < 10e6
    assert len(summary["patterns"]) == 3
    assert {pattern["pattern"] for pattern in summary["patterns"]} == {1, 2, 3}

    with np.load(HOME / "results/spectra.npz") as spectra:
        frequency = spectra["frequency_hz"]
        df = frequency[1] - frequency[0]
        assert frequency[0] == 0
        assert frequency[-1] >= 50e9
        for pattern in summary["patterns"]:
            number = pattern["pattern"]
            for k in (10, 100):
                # Names are deliberately explicit to support external reuse.
                key = f"pattern_{number}_K{k}_w_per_hz"
                density = spectra[key]
                assert density.shape == frequency.shape
                assert np.all(np.isfinite(density)) and np.all(density >= 0)
                actual = float(density.sum() * df)
                expected = pattern["estimates"][str(k)]["integrated_full_band_psd_w"]
                np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=1e-15)
                assert (HOME / f"results/pattern_{number}_K{k}.png").is_file()

    report = PdfReader(HOME / "report/hw04_report.pdf")
    metadata = report.metadata or {}
    assert metadata.get("/Author", "") in ("", "Anonymous"), "PDF author must be anonymous."
    assert len(report.pages) >= 6
    assert all((page.extract_text() or "").strip() for page in report.pages)
    independent = json.loads((HOME / "results/verification.json").read_text(encoding="utf-8"))
    assert independent["status"] == "passed"
    assert independent["unit_test_count"] > 0
    assert independent["unit_test_failures"] == independent["unit_test_errors"] == 0
    assert len(independent["spectral_comparisons"]) == 6
    assert {(item["pattern"], item["K"]) for item in independent["spectral_comparisons"]} == {
        (pattern, k) for pattern in (1, 2, 3) for k in (10, 100)
    }
    for comparison in independent["spectral_comparisons"]:
        assert comparison["relative_L2_error"] < independent["tolerances"]["spectral_relative_L2"]
        assert comparison["relative_full_band_power_error"] < independent["tolerances"]["full_band_power_relative"]
    print(json.dumps({"instruction_files_verified": 5,
                      "saved_spectra_verified": 6,
                      "report_pages": len(report.pages),
                      "pdf_author_anonymous": True,
                      "independent_verification_present": bool(independent)}, indent=2))


if __name__ == "__main__":
    verify()
