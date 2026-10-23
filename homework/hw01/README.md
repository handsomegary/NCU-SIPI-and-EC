# HW1 - Spectrum analysis with a Kaiser window

The task is to infer the number, amplitudes, and frequencies of sinusoidal
components from the supplied finite time record. The analysis compares multiple
Kaiser windows, fits the original samples, and verifies the reconstructed waveform.

## Reproduce

From the repository root:

```sh
python -m pip install -r requirements.txt
python homework/hw01/analyze.py
python -m unittest discover -s tests -v
```

The default input is `instructions/HW1_Kaiser.txt` relative to `analyze.py`.
Results are written to `results/` in the same homework directory. Custom input
and output paths can be selected using the options shown by `--help`.

## Corrected result

The supplied record contains three sinusoidal components:

| Component | Frequency (MHz) | Amplitude (V) | Amplitude (dB re 1 V) |
| --- | ---: | ---: | ---: |
| 1 | 525 | 1 | 0 |
| 2 | 602 | 1 | 0 |
| 3 | 888.5 | 0.003162277660 | -50 |

Small alpha values help distinguish the two nearby strong components. Larger
alpha values suppress their sidelobes enough to reveal the weak component, while
widening the strong components' main lobes. A single window setting is therefore
insufficient for interpreting every feature.

The program obtains its estimates from the data rather than using this table as
input. The numerical output includes the reconstruction residual and fit-order
comparison. Gain compensation calibrates isolated sinusoidal peaks; interference
between overlapping main lobes can still shift measured peak positions and heights.

## Read the work

- [Assignment and input data](instructions/)
- [Rewritten analysis](analyze.py)
- [Numerical output and figures](results/)
- [Anonymous PDF report](report/report.pdf)
- [Report source](report/)

The previous two-tone reconstruction is incomplete. This report and analysis
correct it by detecting and checking the weak third component against the original
samples. Supplied course references and the rewritten solution are kept in
separate directories for clear attribution.
