# HW2 Step response using IDFT and a Kaiser window

National Central University, EE8091 Introduction to Signal Integrity and
Electromagnetic Compatibility.

This assignment reconstructs a reflected step response from measured frequency
samples. Part I relates the Kaiser parameter to the 10-90% input rise time. Part
II applies that window to the supplied S11 of a Chebyshev low-pass filter and
compares the reconstruction with an independently simulated circuit.

## Run the verified analysis

From the repository root, with Python 3.12 or newer:

```sh
python -m pip install -r homework/hw02/requirements.txt
python homework/hw02/src/hw2.py
python -m unittest discover -s homework/hw02/tests -v
```

The default run uses the saved, high-precision ngspice reference in `simulation/`.
To regenerate that reference from the standalone netlist, use an installed
ngspice executable:

```sh
python homework/hw02/src/hw2.py --ngspice ngspice
```

`--output PATH` sends all new results to another directory. No graphical desktop
is required. The supplied Windows ngspice distribution is also retained under
`instructions/ngspice-42_64/`; an absolute path to its console executable can be
passed to `--ngspice`.

## Contents

| Directory or file | Purpose |
| --- | --- |
| `instructions/` | Every file from the supplied Instruction folder, including example files, simulator archives, the expanded distribution, and its licenses |
| `reference/` | Anonymous copies of the historical source variants, data, and ZIP archives; retained as reference rather than the verified implementation |
| `data/chebyshev_lpf.s2p` | S-parameter records used by the rewrite, with a generic header |
| `src/hw2.py` | Portable analysis, parser, numerical methods, figures, and verification |
| `simulation/` | Standalone circuit and freshly computed high-precision reference |
| `results/` | Numerical CSV files, verification JSON, and figures |
| `tests/` | Numerical and input-validation tests |
| `report/` | Anonymous report, Markdown source, and PDF builder |
| `source_manifest.json` | Source coverage, anonymous destination names, sizes, and SHA-256 checksums |

The original identified report is replaced by the rewritten anonymous report.
Regenerable Python bytecode caches are omitted. Archived source member names and
source text are anonymized as needed. The source folders remain unchanged.
Third-party simulator credits and licenses remain intact.

## Main corrections and verification

- The Kaiser parameter is solved numerically for a 40 ps 10-90% rise time,
  rather than chosen from a rough plot. The solution is approximately 1.26911
  for the supplied 20 GHz one-sided bandwidth, with F = 40 GHz.
- The continuous inverse transform uses the correct normalized sinc argument.
  The report explains the inconsistent printed expression in the assignment.
- The IFFT uses explicit conjugate symmetry, continuous-transform scaling,
  endpoint weights, and a centered time axis that includes negative time.
- The SPICE reference is rerun with high-precision output. Incident voltage is
  removed using V(input) minus V(source)/2.
- A 250 ps one-way transmission line gives a 500 ps reflection round trip.
- Kaiser and piecewise-linear inputs have different shapes even when their
  10-90% rise times match. The report separates this comparison from independent
  checks using matching inputs.
- The provided S11 matches the stated circuit. Isolated S21 phase anomalies are
  reported and preserved; S21 is not used to reconstruct the reflection.

Numerical tolerances and the actual errors are recorded in
`results/verification.json`. The program exits unsuccessfully if its checks fail.

## Rebuild the report

```sh
python -m pip install -r homework/hw02/report/requirements.txt
python homework/hw02/report/build_report.py
```

The report includes the complete rewritten analysis source and circuit netlist
as an appendix. Names, student IDs, and personal author metadata are omitted.

## Provenance

The assignment, example files, simulator package, and supplied S-parameter data
are course reference material. The verified implementation and report are
rewritten from those inputs. Historical reference scripts are not used as
numerical verification or represented as new code. Course and third-party
material retain their original ownership and license terms.

## Course attribution

Course materials, supplied simulation data, and the supplied `snpreader.py`
helper are credited to **Prof. Chiu-Chih Chou (周求致), Department of Electrical
Engineering, National Central University, Taiwan**. See the
[full attribution](../../ATTRIBUTION.md) and [input inventory](../../SIMULATION_INPUTS.md).
