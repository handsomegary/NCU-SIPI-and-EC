# SIPI and EC

Coursework for **National Central University, EE8091: Introduction to Signal
Integrity and Electromagnetic Compatibility**.

HW1 and HW2 are included now. The same repository is intended to hold HW1-HW14 as the
remaining assignments are completed.

| Homework | Topic | Contents |
| --- | --- | --- |
| [HW1](homework/hw01/) | Spectrum analysis with a Kaiser window | Assignment, data, rewritten analysis, verification, figures, anonymous report |
| [HW2](homework/hw02/) | Step response via IDFT with a Kaiser window | Complete supplied files, rewritten code, independent SPICE and circuit checks, figures, anonymous report |

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

## Files

- `homework/hw01/instructions/`: supplied assignment, waveform data, reference
  example, and checksums for the uploaded copies.
- `homework/hw01/analyze.py`: rewritten, portable analysis and parameter estimation.
- `homework/hw01/results/`: reproducible figures and numerical results.
- `homework/hw01/report/`: anonymous report and its source.
- `tests/`: automated numerical checks.
- `homework/hw02/`: complete HW2 material, its verified implementation and tests,
  numerical results, and anonymous report with a full source appendix.

## Method and provenance

The assignment, waveform data, and reference example were supplied as course
materials. Their role is distinguished from the rewritten analysis. The reference
example is retained for comparison; the verified estimates come from the data and
the rewritten program. Source files in the original local folders are preserved.

The HW1 analysis implements the assignment's Kaiser equation directly, compensates
the coherent window gain, and compares several alpha values. It then checks the
estimated sinusoidal model against the original samples. Zero padding refines
the frequency grid; observation time and the window determine spectral resolving
power.

HW2 reconstructs the reflected step from S11 with a correctly scaled Hermitian
IFFT. It solves the Kaiser parameter for the requested 40 ps rise, distinguishes
the transmission line's one-way and round-trip delay, and validates against both
fresh SPICE output and independent circuit equations. The report records source
data limitations and the difference between Kaiser and piecewise-linear inputs.

The report, rewritten source, generated results, and commit metadata omit personal
names, student IDs, email addresses, and local machine paths. Repository ownership
is determined by the GitHub account hosting it.
