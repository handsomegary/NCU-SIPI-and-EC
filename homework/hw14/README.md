# HW14 - NRZ versus PAM4

Reproducible comparison of NRZ and PAM4 over the supplied 50 cm differential
channel, with transmitter FFE, receiver CTLE, decision-feedback equalization,
and the assignment's absolute sinusoidal and random timing jitter.

The [anonymous report](report/report.pdf) contains the four final eye diagrams,
equalizer settings, measured openings, verification results and limitations.
Its [Markdown source](report/report.md) and [report builder](report/build_report.py)
are included. No student name, identifier, email or machine-specific path is
included in the report, generated figures or report metadata.

## Recorded comparison

The default experiment uses 65,536 common information bits per case, seed 2025,
64 samples per symbol UI and a 40 ns retained channel kernel.

| Bit rate | NRZ worst EH | PAM4 worst EH | NRZ worst EW | PAM4 worst EW |
| --- | --- | --- | --- | --- |
| 20 Gbps | 840.467 mV | 444.136 mV | 35.401 ps | 35.001 ps |
| 50 Gbps | 16.532 mV | 67.334 mV | 6.883 ps | 7.490 ps |

All four baseline eyes are open, with no bit errors observed in their retained
records. The eye-height ordering persists under an independent random seed,
doubled time resolution and an 80 ns kernel: NRZ has the larger height at
20 Gbps, and PAM4 at 50 Gbps.

The 20 Gbps eye widths are nearly tied: their ordering changes with the second
seed. The small baseline NRZ width advantage is therefore not robust. At
50 Gbps PAM4 has the wider eye in **ps**, while NRZ has the wider eye in
**symbol UI** because PAM4's UI is twice as long. The report and JSON retain both
units, all three PAM4 openings, and the complete sensitivity results.

## Reproduce

Run from the repository root with Python 3.12 or newer:

```sh
python -m pip install -r homework/hw14/requirements.txt
python homework/hw14/verify_inputs.py
python homework/hw14/src/hw14.py --validation
python -m unittest discover -s homework/hw14/tests -v
python homework/hw14/tests/verify_results.py
python -m pip install -r homework/hw14/report/requirements.txt
python homework/hw14/report/build_report.py
```

The solver reads files relative to its own directory and uses a noninteractive
plotting backend. Use `python homework/hw14/src/hw14.py --help` for bit count,
bit rates, time resolution, impulse duration, seed and output-directory options.
The report builder requires successful saved-result verification before it
generates the PDF.

## Model and fair comparison

| Quantity | Shared convention |
| --- | --- |
| Source | 1 Vpp open-circuit differential voltage |
| Port terminations | 50 ohm at every single-ended source and load port |
| Transfer to loaded receiver | `H = Sdd21 / 2` |
| Transition | 10 ps linear 0%-100% ramp for every change in symbol level |
| Sinusoidal jitter | 1.5 ps peak at 10 MHz |
| Random jitter | Gaussian, 1 ps standard deviation |
| Transmitter FFE | Three symbol-spaced taps; sum of absolute taps = 1 |
| CTLE | One zero and two poles; absolute magnitude bounded by 18 dB |
| DFE | Five feedback taps; actual previous receiver decisions |

The same information bit rate is used for both formats. NRZ carries one bit per
symbol; PAM4 carries two. Therefore PAM4 has half the symbol rate and twice the
symbol UI. A common fraction of UI is **not** used for timing jitter: both formats
receive the same absolute jitter specification in ps.

The S4P ports are `(1, 3)` at the input and `(2, 4)` at the output. With
`Sij = outgoing wave i / incident wave j`,

```text
Sdd21 = (S21 - S23 - S41 + S43) / 2
H_loaded = Sdd21 / 2
```

The second factor of two belongs to the open-circuit Thevenin source convention,
not to the mixed-mode conversion. A matched ideal through channel must produce
half the open-circuit source voltage at the receiver.

PAM4 uses Gray labels `00, 01, 11, 10` for the ascending voltage levels
`-Vpp/2, -Vpp/6, +Vpp/6, +Vpp/2`. FFE filtering is applied to the symbol levels
before the resulting transitions are placed on the jittered transmitter clock.
The equalizer design uses the noiseless channel response before a separate
payload is evaluated. The DFE uses one threshold for NRZ and three for PAM4.
The same tap counts and CTLE resource limits apply to all four cases; coefficients
may differ, as permitted by the assignment.

The equalizer design tests 20 CTLE candidates and solves a worst-pattern linear
program at the nominal sampling phase and at offsets of +/-3 ps, with the
central DFE taps and decision thresholds frozen. This finite phase-stress
objective improves robustness; it does not bound unbounded Gaussian jitter.

## Measurements and interpretation

The four final eyes retain physical receiver-voltage units. No format-specific
rescaling makes one eye appear artificially larger. PAM4's three eyes are all
measured, and the smallest opening is used in the NRZ/PAM4 comparison.

- Eye height is the finite-record gap between adjacent conditional voltage
  envelopes at the design sampling phase. Overlap can produce a negative value.
- Eye width is the contiguous time interval around that phase for which the
  design threshold lies between the corresponding envelopes. Widths are recorded
  in both ps and **symbol** UI.
- The recorded symbol and bit error counts use actual recovered decisions,
  including Gray decoding for PAM4.
- Finite-record extrema are not 3-sigma metrics or BER contours. No observed
  errors does not demonstrate zero BER or a target such as `1e-12`.

The model includes the supplied linear differential channel and timing jitter.
Receiver voltage noise, clock recovery, device nonlinearity and FEC are outside
its scope. The lecture's approximate 9.54 dB loss-difference crossover is a useful
screening rule; the final eyes determine the conclusion for this experiment.

## Files

| Path | Purpose |
| --- | --- |
| `instructions/HW14_NRZ_vs_PAM4.pdf` | Supplied assignment |
| `instructions/HW14_diff_channel_50cm.s4p` | Channel data with identifying export comments removed |
| `instructions/input_manifest.json` | Input provenance and checksums |
| `verify_inputs.py` | Verify packaged inputs and numeric preservation |
| `src/hw14.py` | Portable simulation and equalizer design |
| `tests/` | Independent physical and numerical checks |
| `results/metrics.json` | All simulation settings, eye metrics and sensitivity evidence |
| `results/verification.json` | Independent saved-result verification |
| `results/eye_*.png` | Four final NRZ/PAM4 eye diagrams |
| `report/report.pdf` | Anonymous report |
| `report/build_report.py` | Rebuild report from recorded results |

The supplied assignment is preserved. Only comments carrying a local path and
project identity were removed from the exported channel copy; all numerical
tokens are unchanged. The manifest records that transformation explicitly.

## Course attribution

Course materials and supplied simulation data are credited to **Prof. Chiu-Chih
Chou (周求致), Department of Electrical Engineering, National Central University,
Taiwan**. See the [full attribution](../../ATTRIBUTION.md) and
[input inventory](../../SIMULATION_INPUTS.md).
