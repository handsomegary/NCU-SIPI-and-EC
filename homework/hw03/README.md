# HW3 - Sectional convolution

NCU EE8091: simulate the received waveform for the complete 1,048,576-bit
10 Gbps input through a 20 cm channel, using FFT overlap-add and the supplied S21.
Compare the result with the supplied ADS reference at its actual timestamps.

## Reproduce

From the repository root, with Python 3.10 or newer:

```sh
python -m pip install -r requirements.txt
python homework/hw03/verify_inputs.py
python homework/hw03/analyze.py --verification-study
python -m unittest discover -s tests -v
```

The default calculation uses 50 samples/UI (2 ps), 512 bits per section,
20 ps rise/fall time, a 40 ns retained causal impulse response, and Kaiser
alpha = 0.5 (beta = pi * alpha). All 2^20 input bits are processed. Each block
is zero-padded to at least input-block length + FIR length - 1. Overlapping
output tails are added at the correct original sample positions.

The solver accepts plain or gzipped ADS input. Its output is disk-backed during
the computation; the large temporary waveform is removed on successful completion
unless `--keep-waveform` is selected. The committed results contain compact
figures, four 10 ns numerical comparison windows, and `summary.json`.

```sh
python homework/hw03/analyze.py --help
python homework/hw03/analyze.py --samples-per-ui 100 --output-dir /tmp/hw03-1ps
```

The second example performs a finer 1 ps run without overwriting the committed
2 ps results. Runtime and errors depend on the selected parameters and machine.

## Verification and interpretation

The tests compare overlap-add with independent direct convolution, check a final
partial block, preserve the preceding bit at block boundaries, check rise/fall
ramps and impulse scaling, and validate interpolation at nonuniform ADS times.
The input checksum tool also checks that gzip restores the byte-identical original
470,400,039-byte ADS table.

The supplied bit pattern ends at 104.8576 us. ADS continues to 105 us, beyond the
specified input. Only the interval with known input is used for ADS error metrics;
the report records how many later reference samples are excluded. No time shift,
gain, or voltage offset is fitted to make the comparison agree.

The channel data stop at 40 GHz. Higher frequencies are zero-filled, a frequency
window is applied, and the sampled impulse is truncated to a finite causal FIR.
These approximations introduce residual error. The report lists the checked tail
interval, its peak-relative level, retained DC gain, actual waveform errors, and
measured runtime. It does not claim an exact or identically zero-error ADS match.

The verified full-pattern baseline gives **22.27 mV RMS error** and **75.54 mV
maximum absolute error**, with no fitted time shift, gain, or offset. Convolution
took 4.27 s on the recorded machine. Over the four comparison windows, refining to
1 ps changes the computed waveform by at most 0.290 mV; extending the retained
tail to 60 ns changes it by at most 0.095 mV. These checks support numerical
convergence of the implemented model while leaving a material ADS discrepancy.
The optional study is recorded in `results/summary.json` and reproduced by the
command above; it does not establish the cause of the remaining disagreement.

## Anonymous report

See [the PDF report](report/report.pdf) or [the readable report](report/report.md).
The report and generated figures contain no student name or student ID. The
original personally identifying report is not included. Screen recordings are
excluded.

To rebuild the report after recomputing the results:

```sh
python -m pip install -r homework/hw03/report/requirements.txt
python homework/hw03/report/build_report.py --refresh
```

The assignment template and legacy helper scripts in `instructions/` are supplied
references. Use `analyze.py` for the rewritten executable solution.

## Course attribution

Course materials, supplied simulation data, and the supplied `snpreader.py`
helper are credited to **Prof. Chiu-Chih Chou (周求致), Department of Electrical
Engineering, National Central University, Taiwan**. See the
[full attribution](../../ATTRIBUTION.md) and [input inventory](../../SIMULATION_INPUTS.md).
