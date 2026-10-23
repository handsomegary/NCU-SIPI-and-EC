# NCU SIPI and EC

Coursework for **National Central University, EE8091: Introduction to Signal
Integrity and Electromagnetic Compatibility**.

Verified coursework for the assignments listed below. Each completed homework
includes supplied reference material, portable rewritten code, numerical checks,
and an anonymous report. Further assignments can be added as they are completed.

| Homework | Topic | Contents |
| --- | --- | --- |
| [HW1](homework/hw01/) | Spectrum analysis with a Kaiser window | Assignment, data, rewritten analysis, verification, figures, anonymous report |
| [HW3](homework/hw03/) | Sectional convolution with overlap-add | Assignment, full input/reference data, rewritten solver, numerical checks, four ADS comparisons, anonymous report |

## Reproduce

Python 3.10 or newer is required.

```sh
python -m pip install -r requirements.txt
python homework/hw01/analyze.py
python homework/hw03/analyze.py --verification-study
python -m unittest discover -s tests -v
```

Each analysis reads the supplied data relative to its script and writes figures
and numerical results to its own `results/` directory. A graphical desktop is not
required. See each script's `--help` and homework README for options. The HW3 ADS
reference is included as lossless gzip because its original plain text exceeds
GitHub's ordinary per-file size limit; the solver reads it directly.

## Files

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
- `tests/`: automated numerical checks.

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

The report, rewritten source, generated results, and commit metadata omit personal
names, student IDs, personal email addresses, and local machine paths. The public
GitHub owner account remains visible; anonymizing report content does not conceal
repository ownership. Original course materials retain their existing provenance;
no ownership or blanket license is claimed over supplied reference materials.
