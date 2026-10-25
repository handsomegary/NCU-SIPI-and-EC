# NCU SIPI and EC

Coursework for **National Central University, EE8091: Introduction to Signal
Integrity and Electromagnetic Compatibility**.

Verified coursework for the assignments listed below. Each completed homework
includes supplied reference material, portable rewritten code, numerical checks,
and an anonymous report. HW1-HW8 are consolidated here from both course
repositories, preserving their Git histories. The same repository is intended
to hold HW1-HW14 as the remaining assignments are completed.

| Homework | Topic | Contents |
| --- | --- | --- |
| [HW1](homework/hw01/) | Spectrum analysis with a Kaiser window | Assignment, data, rewritten analysis, verification, figures, anonymous report |
| [HW2](homework/hw02/) | Step response via IDFT with a Kaiser window | Complete supplied files, rewritten code, independent SPICE and circuit checks, figures, anonymous report |
| [HW3](homework/hw03/) | Sectional convolution with overlap-add | Assignment, full input/reference data, rewritten solver, numerical checks, four ADS comparisons, anonymous report |
| [HW4](homework/hw04/) | Spectrum of digital signals | All instruction files, verified analysis, six PSD plots, independent checks, anonymous report |
| [HW5](homework/hw05/) | Fixed-rate and triggered statistical eye diagrams | Complete supplied files, rewritten analysis, independent histogram verification, regenerated figures, anonymous report |
| [HW6](homework/hw06/) | Peak distortion analysis | Complete instructions, rewritten solver, random and supplied-pattern eyes, independent full-record verification, anonymous report |
| [HW7](homework/hw07/) | CTLE design and simulation | All supplied instruction files, portable rewritten solver, independently checked eye height, three waveform comparisons, anonymous report |
| [HW8](homework/hw08/) | Transmitter FFE design and simulation | Assignment, template, channel data, five ZF configurations, waveform snapshots, numerical checks, anonymous report |
| [HW9](homework/hw09/) | Decision-feedback equalization | Complete instructions, channel data, 1/2/3-tap automatic designs, numerical verification, anonymous report |
| [HW10](homework/hw10/) | Differential channel simulation and equalizer design | Complete instructions, rewritten CTLE/FFE/DFE solver, four eyes, independent full-record verification, anonymous report |
| [HW11](homework/hw11/) | Crosstalk across eight coupled channels | Complete supplied inputs, corrected multi-input simulation, all eight receive eyes, independent numerical checks, anonymous report |
| [HW12](homework/hw12/) | Jitter simulation and basic analysis | Complete instructions, corrected streaming simulation, four full-length cases, independent checks, anonymous report |
| [HW13](homework/hw13/) | Jitter fitting and BER extrapolation | Complete instructions, verified HW12 inputs, stable Gaussian tail fits, all eight bathtub curves, sensitivity study, anonymous report |
| [學期總回顧](homework/semester-review/) | SIPI semester review | 原始重點回顧問題、89 題完整繁體中文解答、公式與講義頁碼對照 |

## Reproduce

Use Python 3.12 or newer to reproduce all included assignments.

```sh
python -m pip install -r requirements.txt
python homework/hw01/analyze.py
python homework/hw02/src/hw2.py
python homework/hw03/analyze.py --verification-study
python homework/hw04/src/hw4.py
python homework/hw05/src/hw5.py
python homework/hw06/src/hw6.py
python homework/hw07/src/hw7.py
python homework/hw07/tests/verify_reference.py
python homework/hw08/analyze.py
python homework/hw09/analyze.py
python homework/hw09/verify_results.py
python homework/hw10/src/hw10.py --validation
python homework/hw10/tests/verify_reference.py
python homework/hw11/analyze.py
python homework/hw12/src/hw12.py
python homework/hw13/src/hw13.py
python homework/hw13/tests/verify_results.py
python -m unittest discover -s tests -v
```

Each analysis reads the supplied data relative to its script and writes figures
and numerical results to its own `results/` directory. A graphical desktop is not
required. See each script's `--help` and homework README for options. The HW3 ADS
reference is included as lossless gzip because its original plain text exceeds
GitHub's ordinary per-file size limit; the solver reads it directly.

## Run HW6

```sh
python -m pip install -r homework/hw06/requirements.txt
python homework/hw06/src/hw6.py
python -m unittest discover -s homework/hw06/tests -v
python homework/hw06/tests/verify_results.py
```

The [HW6 guide](homework/hw06/) includes the anonymous report and all four
supplied instruction files. Its corrected causal pulse spectrum and matched
S21/2 transfer give a 0.436735919 V finite-model worst-case opening. Both 64,802-bit
records are processed completely; independent reconstruction verifies every
statistical-eye bin and the conditional boundaries.

## Run HW2

Python 3.12 or newer is required for HW2.

```sh
python -m pip install -r homework/hw02/requirements.txt
python homework/hw02/src/hw2.py
python -m unittest discover -s homework/hw02/tests -v
```

The default run uses the saved independent circuit reference. Add
`--ngspice ngspice` to run a fresh circuit simulation. See the
[HW2 guide](homework/hw02/) for the full source inventory, verification method,
and report rebuilding instructions.

## Run HW4

```sh
python -m pip install -r homework/hw04/requirements.txt
python homework/hw04/src/hw4.py
python homework/hw04/tests/verify_reference.py
```

See the [HW4 guide](homework/hw04/) for all six PSD plots, the independent
verification, and the anonymous report.

## Run HW5

```sh
python -m pip install -r homework/hw05/requirements.txt
python homework/hw05/src/hw5.py
python -m unittest discover -s homework/hw05/tests -p test_hw5.py -v
python homework/hw05/tests/verify_delivery.py
```

The [HW5 guide](homework/hw05/) documents the fixed-rate and triggered methods,
complete input files, regenerated eye diagrams, and anonymous report.

## Run HW11

HW11 can also be verified separately:

```sh
python -m pip install -r homework/hw11/requirements.txt
python homework/hw11/verify_inputs.py
python -m unittest discover -s homework/hw11/tests -v
```

The [HW11 guide](homework/hw11/) documents the port mapping, external RC
terminations, independent input streams, numerical approximations, and report
rebuilding procedure. Its report contains the requested single-transmitter and
all-transmitter comparisons. Verification uses independent circuit equations
and direct convolution; no ADS execution or BER extrapolation is claimed.

## Files

HW10's [guide and anonymous report](homework/hw10/) include a reproducible
two-tap FFE and five-tap DFE design. Three separate 262,144-symbol runs pass
the required EH > 120 mV and EW > 15 ps; the main run gives 144.871 mV and
15.811 ps. A separate full-waveform implementation verifies all four stages,
and halving the time step confirms numerical convergence. These are finite,
noiseless simulation results; the optional ADS comparison was not run.

HW12 is documented in [its guide](homework/hw12/). Its full baseline evaluates
2^20 bits per case at 200 samples/UI using bounded-memory filtering. Independent
tests and saved-result checks can be run with:

```sh
python -m unittest discover -s homework/hw12/tests -v
python homework/hw12/tests/verify_results.py
```

- `homework/hw01/instructions/`: supplied assignment, waveform data, reference
  example, and checksums for the uploaded copies.
- `homework/hw01/analyze.py`: rewritten, portable analysis and parameter estimation.
- `homework/hw01/results/`: reproducible figures and numerical results.
- `homework/hw01/report/`: anonymous report and its source.
- `homework/hw03/instructions/`: supplied HW3 assignment, bitstream, channel,
  compressed ADS reference, and legacy reference helpers.
- `homework/hw03/analyze.py`: rewritten full-pattern overlap-add solver.
- `homework/hw03/results/`: verification metrics and comparison figures.
- `homework/hw03/report/`: anonymous report and its reproducible source.
- `homework/hw06/`: complete HW6 inputs, portable PDA and sectional-convolution
  solver, five regenerated figures, saved numerical arrays, independent checks,
  and a six-page anonymous report with editable source.
- `homework/hw07/`: supplied instructions, verified CTLE analysis, numerical checks,
  generated figures, and an anonymous report with editable source.
- `homework/hw08/instructions/`: supplied FFE assignment, template, channel data,
  and legacy Touchstone helper.
- `homework/hw08/analyze.py`: rewritten five-case transmitter FFE simulation.
- `homework/hw08/results/`: eye diagrams, waveform snapshots, tap values, and checks.
- `homework/hw08/report/`: anonymous report and its reproducible source.
- `homework/hw09/`: supplied DFE instructions and channel data, automatic tap/phase/threshold design,
  three verified summer-output eyes, 1 ns snapshots, independent checks and an anonymous report.
- `tests/`: automated numerical checks.
- `homework/hw02/`: complete HW2 material, its verified implementation and tests,
  numerical results, and anonymous report with a full source appendix.
- `homework/hw05/`: complete HW5 instructions and waveform data, rewritten
  fixed-rate and triggered analysis, tests, results, and anonymous report.

## Method and provenance

The assignments, data, and reference examples were supplied as course materials.
They are distinguished from the rewritten analyses. Reference templates and legacy
helpers are retained for comparison; they are not the verified solutions. Source
files in the original local folders are preserved. Screen recordings and the
original personally identifying reports and student source headers are excluded.

The HW1 analysis implements the assignment's Kaiser equation directly, compensates
the coherent window gain, and compares several alpha values. It then checks the
estimated sinusoidal model against the original samples. Zero padding refines
the frequency grid; observation time and the window determine spectral resolving
power. HW3 reconstructs the finite-rise-time NRZ input, derives a sampled channel
impulse response from S21, convolves each section with zero padding, and adds the
overlapping tails at their original sample positions. Its report documents
truncation and sampling approximations and measures agreement with ADS.

HW2 reconstructs the reflected step from S11 with a correctly scaled Hermitian
IFFT. It solves the Kaiser parameter for the requested 40 ps rise, distinguishes
the transmission line's one-way and round-trip delay, and validates against both
fresh SPICE output and independent circuit equations. The report records source
data limitations and the difference between Kaiser and piecewise-linear inputs.

HW5 folds the received voltage into fixed nominal two-UI intervals and compares
it with alignment to interpolated positive differential-clock zero crossings.
Independent full-data histograms verify every count, the unique sample assignment,
and the per-time-column probability normalization.

HW8 designs transmitter taps by a finite zero-forcing system, normalizes their
absolute sum to the available voltage swing, constructs the post-FFE bit levels
and finite-rise-time waveform explicitly, and applies channel overlap-add
convolution. Its tests separate exact finite-system cancellation from residual
ISI in the full channel response. No ADS validation or BER estimate is claimed
for HW8.

The report, rewritten source, generated results, and commit metadata omit personal
names, student IDs, personal email addresses, and local machine paths. The public
GitHub owner account remains visible; anonymizing report content does not conceal
repository ownership. Original course materials retain their existing provenance;
no ownership or blanket license is claimed over supplied reference materials.

Use **NCU-SIPI-and-EC** as the consolidated course repository. The combined
verification includes the coursework workflow for HW1, HW3, HW6, HW7, HW8, and HW9, plus
the dedicated HW2, HW4, and HW5 workflows.

HW9 evaluates one-, two- and three-tap 0/1 DFE with matched S21/2 scaling.
Every saved decision is checked against the source bits at full channel latency.
Its one-tap worst-case bound is nonpositive even though no errors occur in the
65,536-bit noiseless record; the report distinguishes this from the positive
two-/three-tap bounds and from deliberately mistuned error examples.

## 學期總回顧

[89 題完整解答與講義對照](homework/semester-review/)涵蓋 13 個主題，從傅立葉轉換、
眼圖與等化器，到差動訊號、串音、抖動、PAM4、Scrambler 與 SSC。
[原始重點回顧問題](homework/semester-review/instructions/重點回顧問題.txt)保留原文，
解答依原題順序編號，可在完成各次作業後作為學期總複習。
