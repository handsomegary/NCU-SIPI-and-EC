# HW10 - Differential channel simulation and equalizer design

A portable, independently verified solution for the supplied 40 Gb/s four-port
channel. The final CTLE/FFE/DFE design meets the assignment's **EH > 120 mV** and
**EW > 15 ps** requirements on three separate full-length validation records.

Read the [anonymous eight-page report](report/HW10_report.pdf). Its figures and
tables are generated from the [saved numerical results](results/metrics.json).
No student name, student number, personal source filename, or local account path
is included in the report or rewritten source. The PDF Author field is empty.

## Results

Each row uses 262,144 retained symbols, preceded and followed by 1,024 guard
symbols. The time step is 0.5 ps unless indicated otherwise. Every retained
symbol contributes to the extrema; there is no percentile clipping.

| Validation record | Post-DFE EH (mV) | Post-DFE EW (ps) | Decision errors |
| --- | ---: | ---: | ---: |
| Seed 0 | 144.871238 | 15.810803 | 0 |
| Seed 1 | 141.157745 | 15.966811 | 0 |
| Seed 2 | 143.729051 | 15.881126 | 0 |
| Seed 0, 0.25 ps time step | 144.851627 | 15.811108 | 0 |

The independent implementation reconstructs the transmitter waveform, applies
overlap-add convolution, and runs a serial full-waveform bipolar DFE. It agrees
with the production implementation's four stage measurements, waveform excerpts,
and every retained DFE decision. The [independent audit](results/independent_checks.json)
records 12 passing checks.

## Reproduce

From the repository root, with Python 3.12 or newer:

```sh
python -m pip install -r homework/hw10/requirements.txt
python homework/hw10/verify_inputs.py
python -m unittest discover -s homework/hw10/tests -v
python homework/hw10/src/hw10.py --validation
python homework/hw10/tests/verify_reference.py
```

The default `--validation` run regenerates all figures and numerical artifacts,
then checks two additional seeds and doubled time resolution. It exits with a
nonzero status if any reported compliance check fails. The independent audit
also exits nonzero on a failed check. Paths are resolved relative to the scripts;
no notebook, graphical desktop, or machine-specific path is needed.

For a quick smoke run that preserves the published results:

```sh
python homework/hw10/src/hw10.py --bits 4096 --validation --output /tmp/hw10-smoke
```

Choose any writable output directory on Windows. To rebuild the PDF from the
full saved results:

```sh
python -m pip install -r homework/hw10/report/requirements.txt
python homework/hw10/report/build_report.py
```

The optional design search can also be repeated:

```sh
python homework/hw10/design/optimize_ffe.py
```

It uses seed 70, with 16,384 generated and 14,336 retained training symbols. A
linear program searches five FFE coefficients and clock choices, assuming
correct feedback during design. The constraints require positive signed margins
at selected phase offsets. The selected solution has only two nonzero taps;
the released CTLE adds 0.01 dB headroom below the sampled 6 dB search boundary.
The actual solver uses its own recursive decisions and is checked on held-out
seeds 0, 1, and 2. See [design_search.json](results/design_search.json).

## Signal and measurement conventions

- Input pair: ports (1, 3); output pair: ports (2, 4). The mixed-mode order is
  `[c1, d1, c2, d2]` and
  `Sdd21 = (S21 - S23 - S41 + S43) / 2`.
- Each single-ended source has 2 V peak-to-peak swing. The source differential
  levels are -2 V and +2 V. Matched 50-ohm source loading is applied exactly once
  through `Hchannel = Sdd21 / 2`.
- The UI is 25 ps and the source rise time is 5 ps. The causal trapezoidal
  single-bit pulse, SBR, transmitted waveform, and DFE coefficients share the
  same amplitude, delay, and Kaiser-window convention.
- The supplied channel has 8,001 samples from DC to 80 GHz at 10 MHz spacing.
  The simulation uses Kaiser alpha 1 and a 20 ns causal impulse truncation.
- The DFE subtracts the listed postcursor coefficients multiplied by its past
  decisions. Its threshold is 0 V and its slicer levels are -1 and +1. Feedback
  changes halfway between decision instants, after the required previous
  decision is available.
- EH is the minimum logic-one voltage minus the maximum logic-zero voltage at
  the selected clock. Negative pre-equalization values represent overlapping
  conditional levels, i.e. a closed eye.
- EW is the contiguous interval around the clock for which all retained logic
  ones remain above zero and all retained logic zeros remain below zero.
  Crossings are linearly interpolated between samples. Histogram bins are used
  for visualization, not for calculating EH or EW.

## Equalizer parameters

| Component | Parameter | Selected value |
| --- | --- | --- |
| CTLE | Poles | 40, 40, 40 GHz |
| CTLE | Zero | 2.523829378 GHz |
| CTLE | DC gain | -9.751610265 dB |
| CTLE | Numerator gain parameter | 14.248389735 dB |
| CTLE | True peak gain | 5.990000170 dB at 28.114861 GHz |
| FFE | Causal tap order, delays 0 and 1 UI | `[-0.222313178938, 0.777686821062]` |
| FFE | Main tap / sum of absolute values | index 1 / 1.0 |
| DFE | Five subtraction coefficients (V) | `[0.0399107427444, -0.0110696383484, -0.00378981247593, 0.00230686699334, 0.00229698400827]` |
| DFE | First symbol's absolute clock time | 4069 ps |
| DFE | Clock phase within UI | 19 ps |

The CTLE numerator parameter is not the CTLE's maximum gain. The response is
`Aac*wp2*wp3*(s + (Adc/Aac)*wp1) / ((s+wp1)*(s+wp2)*(s+wp3))`, where the gain
parameters in this expression are linear amplitudes. Its continuous-frequency
maximum is independently checked against the 6 dB limit.

## Included files

- `instructions/`: all three original instruction files, the supplementary
  legacy Touchstone reader, and exact byte-count/SHA-256 provenance.
- `src/hw10.py`: rewritten solver with explicit matched loading, causal pulse
  timing, portable paths, deterministic bitstreams, and aligned DFE timing.
- `tests/test_hw10.py`: nine behavioral tests for mode conversion, source
  loading, pulse shape, linear convolution, FFE, DFE, file ordering, and eyes.
- `tests/verify_reference.py`: independent full-record numerical verification.
- `design/optimize_ffe.py`: reproducible finite-pattern design search.
- `results/`: S21/Sdd21 comparison, four eyes, four-output waveform comparison,
  CTLE/SBR figures, all coefficients and measurements, and compact numerical data.
- `report/`: newly generated anonymous PDF and its editable build source.

The original personal report is replaced by the new report. The supplied
template and legacy reader are reference material rather than the verified
implementation. The rewrite corrects the reference solution's dimensionally
inconsistent SBR phase, mismatched SBR/waveform window settings, dispersed
voltage scaling, and fixed DFE sample phase.

The verification establishes finite, noiseless simulation results under the
stated bandwidth and impulse truncation. It does not establish every possible
bit pattern, a physical BER, jitter tolerance, or nonlinear receiver behavior.
The assignment's optional ADS ChannelSim comparison was not run.
