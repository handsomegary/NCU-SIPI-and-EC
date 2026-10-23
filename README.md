# SIPI and EC

Coursework for **National Central University, EE8091: Introduction to Signal
Integrity and Electromagnetic Compatibility**.

HW1 is included now. The same repository is intended to hold HW1-HW14 as the
remaining assignments are completed.

| Homework | Topic | Contents |
| --- | --- | --- |
| [HW1](homework/hw01/) | Spectrum analysis with a Kaiser window | Assignment, data, rewritten analysis, verification, figures, anonymous report |

## Run HW1

Python 3.10 or newer is required.

```sh
python -m pip install -r requirements.txt
python homework/hw01/analyze.py
python -m unittest discover -s tests -v
```

The analysis reads the supplied data using a path relative to the script and
writes figures and numerical results to `homework/hw01/results/`. A graphical
desktop is not required. See `python homework/hw01/analyze.py --help` for options.

## Files

- `homework/hw01/instructions/`: supplied assignment, waveform data, reference
  example, and checksums for the uploaded copies.
- `homework/hw01/analyze.py`: rewritten, portable analysis and parameter estimation.
- `homework/hw01/results/`: reproducible figures and numerical results.
- `homework/hw01/report/`: anonymous report and its source.
- `tests/`: automated numerical checks.

## Method and provenance

The assignment, waveform data, and reference example were supplied as course
materials. Their role is distinguished from the rewritten analysis. The reference
example is retained for comparison; the verified estimates come from the data and
the rewritten program. Source files in the original local folders are preserved.

The analysis implements the assignment's Kaiser equation directly, compensates
the coherent window gain, and compares several alpha values. It then checks the
estimated sinusoidal model against the original samples. Zero padding refines
the frequency grid; observation time and the window determine spectral resolving
power.

The report, rewritten source, generated results, and commit metadata omit personal
names, student IDs, email addresses, and local machine paths. Repository ownership
is determined by the GitHub account hosting it.
