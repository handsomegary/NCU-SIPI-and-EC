# HW11 - Eight-channel crosstalk simulation

This package solves the supplied 16-port coupled-channel problem with the specified RC terminations, independent 10 Gbit/s NRZ transmitters, finite 20 ps transitions, and alternating 0 / 50 ps source delays. It includes the complete supplied instruction material, a completed solver, regenerated results, independent checks, and an anonymous English report.

## Files

| Path | Purpose |
| --- | --- |
| `instructions/hw11_assignment.pdf` | Original problem specification and instructor example figures. |
| `instructions/HW11_crosstalk_8lines.s16p` | Original 16-port network data, retained byte-for-byte. |
| `instructions/hw11_clarifications.txt` | Supplied clarification notes. |
| `instructions/hw11_incomplete_template.py` | Supplied incomplete template, preserved as a reference. |
| `instructions/legacy_snpreader.py` | Legacy reader reference; not used by the verified rewrite. |
| `instructions/legacy_run_notes.txt` | Legacy manual-switch notes; not current execution instructions. |
| `source_manifest.json` | Original/staged hashes, neutral source mapping, and source PDF metadata review. |
| `analyze.py` | Completed, reproducible solver and result generator. |
| `tests/` | Independent circuit, data-reader, convolution, source, and eye checks. |
| `results/metrics.json` | Configuration, results, and numerical diagnostics for the saved run. |
| `results/figures/` | Regenerated assignment plots and all-eight-eye overview. |
| `report/HW11_Crosstalk_Anonymous_Report.pdf` | Anonymous technical report containing every required output. |
| `report/build_report.py` | Editable report source; consumes the saved numerical results and figures. |

The four files from the original `Instruction` directory are all included under neutral names. The template is intentionally incomplete and contains placeholders and notebook magic; it is not a runnable solver. The original personal student report is excluded. The report in this package is newly written and contains only regenerated figures.

## Reproduce the results

Use Python 3.12 or newer. From this directory:

```text
python -m pip install -r requirements.txt
python verify_inputs.py
python analyze.py
python tests/verify_results.py
python report/build_report.py
```

The default analysis reproduces the saved full record with 262,144 bits per transmitter. Commands use paths relative to this package and do not require the original local download folder. The saved-artifact verifier also runs the unit tests and writes `tests/verification.json`. The report builder requires passing verification that matches the current results and fails if a required result or figure is absent, instead of substituting the instructor examples.

From the repository root, the complete audit can also be run as:

```text
python homework/hw11/tests/verify_results.py
```

To run only the component tests without the saved-artifact audit, use `python -m unittest discover -s homework/hw11/tests -v` from the repository root.

Refer to `python analyze.py --help` for the available analysis options. A short record can exercise the workflow, but its empirical voltage-cluster gaps will not be identical to the full saved run.

## Circuit and indexing

Channel labels in code and figures are zero-based. Physical Touchstone ports remain one-based:

| Channel | Tx physical port | Rx physical port | Source delay |
| --- | --- | --- | --- |
| 0 | 1 | 2 | 0 ps |
| 1 | 3 | 4 | 50 ps |
| 2 | 5 | 6 | 0 ps |
| 3 | 7 | 8 | 50 ps |
| 4 | 9 | 10 | 0 ps |
| 5 | 11 | 12 | 50 ps |
| 6 | 13 | 14 | 0 ps |
| 7 | 15 | 16 | 50 ps |

Every Tx has a 30 ohm source resistor and 1.5 pF shunt capacitance. Every Rx has a 65 ohm load resistor in parallel with 1 pF. The source levels are 0 and 1 V. Turning a Tx off sets its source voltage to zero while retaining its termination.

At each frequency, the solver combines KCL with the multiport relation

```text
(I - S) v - Z0 (I + S) i = 0
```

and solves all eight source right-hand sides together. This includes external loading, reflections, and crosstalk. Using `S/2` alone would describe matched resistive terminals, not the assigned RC circuit. Each loaded impulse response is convolved with its transmitter waveform; superposition produces the all-on outputs.

Separate `SeedSequence` child streams generate the transmitter data. Reusing one seed by resetting the generator for every Tx would produce identical patterns and violate the all-on experiment. The original legacy note suggesting a common seed is retained only for provenance and must not be followed for this test.

## Required outputs

| Assignment requirement | Regenerated figure |
| --- | --- |
| All sixteen `S1j` curves | `results/figures/s_parameters.png` |
| Rx0 with only Tx0 on | `results/figures/eye_rx0_only_tx0.png` |
| Rx0 with all Tx on | `results/figures/eye_rx0_all_on.png` |
| Rx4 with only Tx4 on | `results/figures/eye_rx4_only_tx4.png` |
| Rx4 with all Tx on | `results/figures/eye_rx4_all_on.png` |
| Four-nanosecond Rx1-Rx4 snapshot with only Tx4 on | `results/figures/tx4_only_snapshot.png` |
| Compute every Rx eye | `results/figures/eyes_all_on.png` |

The waveform plot labels any 5x display scaling of the coupled responses. That scaling does not affect the stored waveforms, eye histograms, or metrics.

## Interpretation and verification

The plots are finite-record eye histograms. Their color scale is probability mass per voltage bin and phase, not a probability density, BER, or bathtub extrapolation. An unvisited bin is not proof of zero error probability. The saved `opening_v` diagnostic is an empirical voltage-cluster quantile gap: samples are split at half the main-path DC gain, and the upper cluster's 0.1% quantile minus the lower cluster's 99.9% quantile is maximized over 50 phases, requiring at least 20% of samples in each cluster. It is not bit-conditioned or a worst-case eye opening. The histogram voltage grid is 0.5 mV. Eye phases remain absolute time modulo the UI, without per-receiver recentering.

The finite supplied bandwidth creates impulse pre-ringing and long tails. The rewrite keeps explicit negative and positive impulse support, performs linear FFT overlap-add with guards, and records support-truncation and larger-support comparison diagnostics. Independent tests check nonreciprocal Touchstone ordering, loaded scattering equations, an analytic RC divider, matched-load normalization, exact linear convolution and superposition, waveform delays and finite edges, source independence, and per-phase eye-histogram normalization.

These checks establish numerical consistency for the stated finite-bandwidth model. No ADS execution, hardware measurement, receiver-noise model, random jitter, clock recovery, equalizer optimization, or BER claim is included. The assignment's ADS images are instructor references only.

## Provenance and anonymity

`source_manifest.json` records SHA-256 values for the original and staged reference files without personal source-directory names or absolute local paths. The original assignment PDF was inspected for visible identity and document metadata; it has no author identity and is retained unchanged. The regenerated report has deliberately neutral PDF metadata and no named author, student identifier, email address, or local user path.

## Course attribution

Course materials, supplied simulation data, and the supplied `snpreader.py`
helper are credited to **Prof. Chiu-Chih Chou (周求致), Department of Electrical
Engineering, National Central University, Taiwan**. See the
[full attribution](../../ATTRIBUTION.md) and [input inventory](../../SIMULATION_INPUTS.md).
