# HW13 - Jitter fitting and BER extrapolation

This analysis continues the verified HW12 baseline: 40 Gbps NRZ, 25 ps/UI,
200 samples/UI, and 1,048,576 evaluated bits per case. It computes all eight
empirical Tx/Rx CDFs and bathtub curves, then fits two Gaussian components to
the **tails only** of Case 3 and Case 4 receive jitter.

The [anonymous report](report/report.pdf) and [editable report](report/report.md)
explain the method, corrections, results, and limitations.

## Main results

| Model prediction | Case 3 Rx: RJ + channel DDJ | Case 4 Rx: SJ + RJ + channel DDJ |
| --- | ---: | ---: |
| Minimum BER | 2.95e-13 | 2.30e-8 |
| Optimum sampling time | 24.0038 ps | 24.2310 ps |
| Eye width at BER 1e-6 | 8.2470 ps | 3.1901 ps |
| Eye width at BER 1e-9 | 4.1779 ps | No opening |
| Eye width at BER 1e-12 | 0.8300 ps (0.0332 UI) | No opening |

These are **extrapolated model predictions**, not measured low-BER results.
Case 3 is sensitive to the fitting window: the nine documented windows predict
minimum BERs from about 1.57e-14 to 2.20e-11, including windows with no 1e-12
opening. Case 4 has no 1e-12 opening in all nine examined windows. Window
sensitivity is not a confidence interval or a certification of link performance.

## Reproduce and verify

From the repository root, with Python 3.12 or newer:

```sh
python -m pip install -r requirements.txt
python homework/hw13/verify_inputs.py
python -m unittest discover -s homework/hw13/tests -v
python homework/hw13/src/hw13.py
python homework/hw13/tests/verify_results.py
```

The analysis also works from another current directory. Input paths are resolved
relative to the script. To keep the committed results intact:

```sh
python homework/hw13/src/hw13.py --output work/hw13-experiment
python homework/hw13/tests/verify_results.py --results work/hw13-experiment
```

The primary input is a frozen, checksummed copy of the verified HW12 export.
HW13 reuses that saved full-length simulation; it does not rerun the channel.
The supplied older CSV is separately corrected and analyzed on every run for
comparison. Both datasets' empirical distributions are exported, along with
legacy fit metrics and the predictions of both original parameter sets.

## Method

All analysis times and Gaussian parameters are in ps. PMF values sum to one;
the PDF is PMF / 0.125 ps and has units 1/ps. A histogram value represents mass
in the cell centered on its physical sample time. CDFs are evaluated at upper
cell edges, survival functions (SFs) at lower edges, and interpolated assuming
uniform density within a cell. This makes the half-bin convention explicit.

With nominal boundaries at 12.5 and 37.5 ps and transition density 0.5:

`BER(t) = 0.5 * SF(t - 12.5) + 0.5 * CDF(t - 37.5)`.

The model has five independent parameters: mu1, sigma1, mu2, sigma2, and tau1;
tau2 = 1 - tau1. Positive sigma and valid weights are enforced. Direct Gaussian
CDF/SF evaluations in log space avoid subtracting a nearly-one CDF from one.
Empirical SFs use reverse cumulative sums for the same reason.

The baseline fit uses single-sided cumulative tail probabilities from 3e-4 to
3e-2. This includes at least about 41 accumulated eye-slice observations in
each fitted tail. Six deterministic starting points minimize unweighted log10
left-CDF and right-SF residuals, with bounded means, widths, and mixture weight.
The central DJ region is excluded. The selected tails' PDF agreement is checked
separately. Cumulative samples are correlated; the residual RMSE is descriptive,
not an independent-sample statistical error bar.

The minimum is located by a grid search followed by local refinement. Low-BER
crossings are solved with bracketed roots, including a refined minimum so an
opening narrower than the scan spacing is not missed. A target below the model
minimum is reported explicitly as no opening with width zero.

## Corrections and provenance

- The two supplied legacy programs have different manual fit parameters. Both
  are preserved anonymously as reference text, with their model predictions
  recorded in `results/summary.json`; neither is treated as a reference answer.
- The older export inherited a known 1e-15 plotting offset from the supplied
  HW12 routine. The legacy import subtracts that offset and preserves zero mass.
  The primary verified HW12 export contains no artificial offset.
- The legacy 200-point endpoint-inclusive time labels have spacing 25/199 ps.
  The physical 200-bin grid is instead `(index - 100) * 0.125 ps`.
- The rewrite fixes tail cancellation, inconsistent time units, invalid log-axis
  limits, missing numerical eye-width extraction, and current-directory paths.
- No artificial probabilities are added for plots. Empty cells are masked on
  log axes, and zero empirical BER is not interpreted as zero physical BER.

Eye-slice occupancy can weight slower edges differently from counting each
threshold crossing once. The one-UI histogram also folds extreme jitter. The
assignment's common-distribution, two-boundary approximation and assumed 0.5
transition density are retained. No ADS or independent low-BER measurement is
claimed.

## Files

- `instructions/`: both supplied instruction files, byte-for-byte, and the input
  SHA-256 manifest. The starter Python file is intentionally incomplete and is
  not executable analysis code.
- `inputs/`: verified HW12 snapshot and provenance, plus the supplied legacy CSV.
- `reference/`: both supplied programs under neutral filenames as reference text.
- `src/hw13.py`: portable analysis, tail fitting, root finding, and figure export.
- `tests/`: numerical unit checks and an independent saved-results verifier.
- `results/`: all eight CDFs/PDFs/bathtubs, two fit comparisons, window sensitivity,
  CSV curves, JSON metrics, and verification evidence.
- `report/`: a newly authored anonymous PDF, editable source, and report builder.

The original personally identifying reports and source archive names are
replaced by the new anonymous report and neutral reference names. Source
materials retain their course provenance. Repository ownership remains visible
on GitHub; anonymity applies to the submission content and new commit metadata.

Report regeneration:

```sh
python -m pip install -r homework/hw13/report/requirements.txt
python homework/hw13/report/build_report.py --refresh
```
