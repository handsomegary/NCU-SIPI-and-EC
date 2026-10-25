# Supplied HW11 reference files

The four files from the supplied `Instruction` directory are retained byte-for-byte under neutral names. See `../source_manifest.json` for hashes and mapping. The assignment PDF contains the problem statement and instructor-provided example figures; those examples are not results of this repository.

| File | Purpose |
| --- | --- |
| `hw11_assignment.pdf` | Authoritative assignment, port map, circuit values, required outputs. |
| `HW11_crosstalk_8lines.s16p` | Original 16-port dataset; odd Tx ports and even Rx ports. |
| `hw11_clarifications.txt` | Original clarification notes, including multi-right-hand-side MNA and array dimensions. |
| `hw11_incomplete_template.py` | Incomplete supplied reference template. Deliberately not executable; use `../analyze.py`. |
| `legacy_snpreader.py` | Legacy support reader retained for provenance; not imported by the verified solver. |
| `legacy_run_notes.txt` | Legacy manual switches, not current run instructions. Its same-seed suggestion must not be used for the all-transmitter experiment. |

The original incomplete template contains syntax placeholders and notebook magic. The verified rewrite provides the completed method rather than silently filling or altering the reference. The legacy support files are retained only to make the source package complete. Follow the HW11 README for the reproducible workflow.

The original personal report is intentionally excluded. All report figures in `../report/` are regenerated from the solver outputs. No ADS execution is claimed.
