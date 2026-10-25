# HW8 - Transmitter FFE design and simulation

NCU EE8091: design zero-forcing (ZF) transmitter feed-forward equalizers for the
supplied 50 cm stripline channel and compare their received eyes at 40 Gbps NRZ.
Construct the post-FFE waveform in the time domain, including finite rise/fall
time, then convolve it with the channel response.

| Case | Precursor taps | Main taps | Postcursor taps |
| --- | ---: | ---: | ---: |
| 1 - Unequalized | 0 | 1 | 0 |
| 2 - Precursor only | 1 | 1 | 0 |
| 3 - Postcursor only | 0 | 1 | 1 |
| 4 - Precursor and postcursor | 1 | 1 | 1 |
| 5 - Precursor and two postcursors | 1 | 1 | 2 |

## Reproduce

From the repository root, with Python 3.10 or newer:

```sh
python -m pip install -r requirements.txt
python homework/hw08/verify_inputs.py
python homework/hw08/analyze.py --bits 65536 --seed 0
python -m unittest discover -s tests -v
```

The default calculation uses a 25 ps UI, 5 ps rise/fall time, 0.5 ps samples
(50 samples/UI), and a 0 to 2 V input before equalization. The channel transfer
function is `H = S21 / 2`; the matched-source factor is applied once. Each FFE tap
vector is normalized so that the sum of its absolute values is one. All cases
use the same seeded 65,536-bit input, so their eyes can be compared directly.

The solver reads the bundled channel by default. To use an explicit input path
and a separate result directory:

```sh
python homework/hw08/analyze.py --help
python homework/hw08/analyze.py --data homework/hw08/instructions/HW7_stripline_50cm.s2p --output work/hw08-results --bits 65536 --seed 0
```

`--output-dir` is accepted as an alias for `--output`. Change the output directory
when experimenting to preserve the committed baseline results.

## Results and interpretation

The `results/` directory contains five individual received-eye figures and an
overview, channel and single-bit responses, tap values, equalized single-bit
responses, `metrics.csv`, and a machine-readable `summary.json`. All eye plots
use the same voltage limits, histogram bin edges, and color scale. Case 4 also
contains two Tx snapshots and two Rx snapshots, as requested by the assignment,
with a CSV of the plotted samples for each. Snapshot axes show time relative to
each window's start. The summary records their absolute simulation start times
and source payload positions; the Rx windows include the modeled channel delay.

ZF solves a finite cursor system. Its selected matrix residual measures agreement
with that system; it does not imply that every cursor of the physical channel has
vanished. The summary reports residual physical ISI separately. The eye metric
is the empirical vertical gap for the simulated bit sequence at the chosen
sampling phase, rather than a bit-error-rate estimate. Finite channel bandwidth,
windowing, impulse truncation, and the finite input sequence limit the model.
The report explains the numerical checks and resulting comparisons. The optional
ADS ChannelSim exercise is outside this implementation.

## Anonymous report

See [the PDF report](report/report.pdf). The rewritten source, report, generated
figures, and result tables contain no student name or student ID. The original
personally identifying solution and report are excluded, along with recordings
and Python cache files.

To rebuild the report after recomputing the results:

```sh
python -m pip install -r homework/hw08/report/requirements.txt
python homework/hw08/report/build_report.py --refresh
```

The assignment template and legacy helper in `instructions/` are supplied
references. Use `analyze.py` for the rewritten executable solution. The reference
manifest and `verify_inputs.py` record and check the original input bytes.

## Course attribution

Course materials, supplied simulation data, and the supplied `snpreader.py`
helper are credited to **Prof. Chiu-Chih Chou (周求致), Department of Electrical
Engineering, National Central University, Taiwan**. See the
[full attribution](../../ATTRIBUTION.md) and [input inventory](../../SIMULATION_INPUTS.md).
