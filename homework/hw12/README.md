# HW12 - Jitter simulation and basic analysis

Four 40 Gbps NRZ simulations compare no transmitter jitter, sinusoidal jitter
(SJ), Gaussian random jitter (RJ), and their combination through the supplied
loaded two-port channel. The full baseline evaluates **1,048,576 steady-state
bits per case**, at 200 samples/UI and a 1 mV eye-bin spacing.

## Files

- `instructions/`: all three supplied assignment files, preserved byte for byte,
  plus their SHA-256 manifest. The instructor Python file is an incomplete
  teaching template, not the executable solution.
- `src/hw12.py`: corrected, portable implementation with a command-line interface.
- `reference/legacy_reference.py.txt`: anonymized supplied implementation for
  comparison. It retains known errors and is not the verified solution.
- `tests/`: independent circuit, waveform, convolution, statistics, and saved
  result checks.
- `results/`: regenerated eyes, PMF/PDF plots, numerical CSVs, raw eye counts,
  full-run parameters, and verification evidence.
- `report/`: newly written anonymous PDF report, editable Markdown, and builder.

The original reports containing personal identifiers are replaced by the new
report. No student name, student number, or personal filesystem path is included
in the rewritten solution or its report.

## Reproduce

Run from the repository root with Python 3.12 or newer:

```sh
python -m pip install -r requirements.txt
python homework/hw12/verify_inputs.py
python -m unittest discover -s homework/hw12/tests -v
python homework/hw12/src/hw12.py
python homework/hw12/tests/verify_results.py
```

To refresh the figures from the saved counts and CSVs without rerunning the
simulation, use `python homework/hw12/src/hw12.py --render-only`.

The default writes `homework/hw12/results/`. To preserve the committed baseline
while experimenting, provide a separate output directory:

```sh
python homework/hw12/src/hw12.py --bits 1048576 --points-per-ui 200 --block-bits 4096 --seed 0 --output work/hw12-experiment
```

A short run such as `--bits 4096` is useful for a smoke check but does not meet
the assignment's sample-count requirement or adequately cover the 1 MHz SJ
cycle. A full baseline covers approximately 26.2 SJ cycles. No ADS validation
or low-BER extrapolation is claimed.

## Physical model

| Quantity | Value |
| --- | ---: |
| Bit rate / unit interval | 40 Gbps / 25 ps |
| Source swing / linear ramp duration | 1 V / 5 ps |
| Source resistor / shunt capacitor | 20 ohm / 0.3 pF |
| Load resistor / shunt capacitor | 85 ohm / 0.3 pF |
| SJ peak amplitude / frequency | 3 ps / 1 MHz |
| RJ mean / standard deviation | 0 ps / 1 ps |
| Time step / voltage bin | 0.125 ps / 1 mV |
| Kaiser parameter / retained impulse | 0.5 / 40 ns |

The Touchstone reader preserves the standard S11,S21,S12,S22 order, including
nonreciprocal networks. A modified nodal solve includes both terminations and
capacitors; the loaded transfer function is not replaced by S21. The
DC-correct inverse transform produces a discrete impulse kernel, so an extra
time-step factor must not be applied during discrete convolution.

The source waveform uses linear ramps shifted by the sum of SJ and RJ, rounded
to the nearest sample with `np.rint`. The maximum rounding error is half a time
step. Requested and applied jitter moments are recorded separately. Transitions
that would physically overlap are rejected rather than silently overwritten.
Independent seeded random streams provide the same data and RJ realization
across the four cases and across processing block sizes.

The simulator uses bounded-memory overlap-save filtering and incremental eye
histograms. It discards a full FIR-length warmup and evaluates exactly the
requested number of bits, without padding a zero tail into the statistics.

## Jitter distributions and variance

For each case, the Tx threshold is 0.5 V and the Rx threshold is half the loaded
DC output swing. A horizontal row of the statistical eye is extracted for one
UI, circularly aligned to the baseline case on the same side, and normalized to
a PMF. Zero cells remain zero. The PDF equals PMF divided by the physical time
step, and its integral is one. Plotting masks zero cells only for logarithmic
display.

The CSVs use seconds and density per second; figure axes convert these to ps and
density per ps. Variance is computed after subtracting the mean. Reported RMS
jitter is the centered RMS, equal to the standard deviation.

An independently interpolated threshold-crossing histogram is also saved.
It counts crossing events equally, while an eye slice measures voltage-bin
occupancy and can weight slower edges more strongly. These estimators are
reported separately rather than treated as identical.

The continuum Tx standard-deviation benchmarks are 0, 2.1213, 1, and 2.3452 ps.
Finite phase coverage, finite random samples, and nearest-grid quantization
explain small differences from those benchmarks. Rx variance addition is not
an unconditional identity: the channel can change how jitter affects crossings.
The report answers all four questions and explicitly discusses the apparent
Tx/Rx typo in the fourth question.

## Anonymous report

Read [the report](report/report.pdf) or its [editable source](report/report.md).
After a verified full simulation, rebuild it with:

```sh
python -m pip install -r homework/hw12/report/requirements.txt
python homework/hw12/report/build_report.py
```

The report builder requires the full baseline and successful independent
verification evidence. See `results/test_results.json` and
`results/verification.json` for the recorded checks. The additional
`results/block_size_check.json` records identical eye and crossing counts for
two different processing block sizes in an independent reproducibility run.
