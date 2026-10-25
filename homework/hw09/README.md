# HW9 - DFE simulation and design

Verified 40 Gbps NRZ simulation for the supplied 50 cm stripline. The matched
source/load transfer is **H = S21/2** with a 0-2 V source, 5 ps rise time and
25 ps UI. One-, two- and three-tap 0/1 decision-feedback equalizers are designed
automatically and evaluated on the same 65,536-bit record.

[Anonymous report](report/report.pdf) · [Numerical summary](results/summary.json) ·
[Independent verification](results/verification.json)

## Reproduce

From the repository root, using Python 3.12 or newer:

~~~sh
python -m pip install -r requirements.txt
python homework/hw09/verify_inputs.py
python -m unittest discover -s tests -p test_hw09.py -v
python homework/hw09/analyze.py
python homework/hw09/verify_results.py
~~~

The analysis takes `--bits`, `--seed` and `--output` options. Its default seed is
20261025; the input is seeded independent Bernoulli bits, not an LFSR PRBS.
The report uses the committed default results. To rebuild it:

~~~sh
python -m pip install -r homework/hw09/report/requirements.txt
python homework/hw09/report/build_report.py --refresh
~~~

## What is included

- `instructions/`: both files supplied in the HW9 Instruction directory, plus
  the original legacy reader and note as reference material. The template is
  intentionally incomplete; it is not imported or run by the verified solver.
- `data/`: the unchanged channel Touchstone file.
- `analyze.py`: corrected DFE simulation, automatic design, figures and diagnostics.
- `channel.py`: channel/pulse/overlap-add utilities adapted from this repository's
  verified HW8 implementation; no transmitter FFE is applied.
- `results/`: three eyes, three 1 ns waveform comparisons, every scored decision,
  all 50 phase candidates, a mistuned-threshold example and verification evidence.
- `report/`: anonymous PDF, readable report source and reproducible report builder.
- `../../tests/test_hw09.py`: independent mathematical regression tests.

All five supplied files are listed with their SHA256 digests in
`instructions/manifest.json`. Original personally identifying student code and
reports are replaced by the rewritten code and anonymous report. Generated
bytecode and local paths are excluded.

## Design and conventions

The discrete source uses 50 samples/UI and a 10-sample linear transition.
The measured response (DC through 80 GHz) is Kaiser-windowed with alpha 0.5,
zero above 80 GHz, transformed with Hermitian endpoint accounting, and truncated
to a 20 ns FIR. A single-bit response is computed from the exact same sampled
pulse and FIR used for random-data simulation.

For full sampling latency L, cursor a_j = g[L+jP], where P=50. Positive j is a
postcursor. Additive 0/1 feedback uses c_j=-a_j in volts. Coefficients are not
normalized. After cancelling the selected postcursors, retain every other
precursor and postcursor in the finite model. If r denotes those residuals:

~~~text
zero_upper = sum(max(r, 0))
one_lower  = a0 + sum(min(r, 0))
threshold  = (zero_upper + one_lower) / 2
eye_bound  = a0 - sum(abs(r))
~~~

Search exactly 50 consecutive L values near the SBR maximum and select the
largest bound; the earliest L wins ties. The midpoint is still a minimax
threshold when the bound is nonpositive, but then no all-pattern separation
is established. Bounds assume previous decisions are correct.

The slicer uses `z[k] >= threshold`. A decision's first feedback tap is held from
`k+P//2` through `k+P//2+P-1`, covering the next sampling instant. The first
postcursor is `taps[0]`. Source bit k is scored at `k*P+L`, preserving full
propagation latency, not just L modulo P. The eye is the summer output before
the slicer.

## Default result

| Taps | Phase / 50 | Full latency (ps) | Threshold (V) | Model eye bound (V) | Observed eye (V) | Errors / bits |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 2900.5 | 0.389388930 | -0.000824118 | 0.073723539 | 0 / 65,536 |
| 2 | 0 | 2900.0 | 0.343890061 | 0.086382991 | 0.160946396 | 0 / 65,536 |
| 3 | 49 | 2899.5 | 0.316287348 | 0.136673606 | 0.206683012 | 0 / 65,536 |

The one-tap result has a negative worst-case bound despite no errors in this
finite record. It therefore does **not** establish separation for every possible
bit history. No physical BER=0 claim is made. Noise, jitter and circuit
nonidealities are not modeled.

All waveform CSVs include both endpoints: 2,001 samples from 0 through 1,000 ps.
Each eye folds all 65,536 payload UIs. Sampling appears at 0.5 UI and 1.5 UI;
eye voltage ranges adapt to each case and share the same logarithmic color scale.

The separately labeled three-tap stress run raises the threshold to
0.547061883 V and produces 14,439 errors. Its first payload error is bit 1:
true bit 1 is detected as 0 at 0.508835022 V. A separate forced flip at payload
bit 512 changes the next three summer voltages by the corresponding tap values;
no further bit-decision changes occur in that particular diagnostic.

## Verification and limitations

Twenty independent unit tests cover scalar versus held-block feedback,
threshold equality, signed taps, delayed truth alignment, exhaustive mixed-sign
residual bounds, phase selection, channel scaling and direct convolution.

The solver also compares every decision against an independent scalar recurrence
and reconstructs every post-channel payload sample from source bits and channel
cursors. `verify_results.py` independently checks the saved evidence against the
input data and seed. Errors and floating-point tolerances are recorded in JSON.

The continuous-time analytic pulse spectrum uses
`exp(-j*pi*f*(Tb+Tr))`; it is only a cross-check, not mixed into the discrete tap
design. Increasing FIR length from 20 to 40 ns changes the eye bound by about
-0.50 mV at the selected latency. ADS ChannelSim is optional and was not run.

## Provenance and privacy

The assignment, starter template, channel and legacy helper are supplied course
materials. DFE conventions follow the course DFE handout, pages 4-10. Automatic
threshold/phase selection and verification are derived and implemented here.
The original local sources are preserved. No blanket license or ownership over
supplied course materials is asserted.

The rewritten source, report text, report metadata and this contribution's commit
metadata omit personal names, student IDs, personal email addresses and local
machine paths. GitHub repository ownership and account activity remain visible.

## Course attribution

Course materials, supplied simulation data, and the supplied `snpreader.py`
helper are credited to **Prof. Chiu-Chih Chou (周求致), Department of Electrical
Engineering, National Central University, Taiwan**. See the
[full attribution](../../ATTRIBUTION.md) and [input inventory](../../SIMULATION_INPUTS.md).
