# Supplied HW1 material

- `HW1_time_freq_transf_(1)_Kaiser.pdf`: original assignment content with anonymous
  document metadata.
- `HW1_Kaiser.txt`: unchanged two-column input data, with time in ns and signal in V.
- `HW1_Kaiser_sample.py`: unchanged course reference example. It contains a
  machine-specific input path and a manually chosen two-tone reconstruction.
  Those choices are not the verified solution.
- `manifest.json`: SHA-256 hashes of the three uploaded reference files.

Run `../analyze.py` for the portable rewritten analysis and the complete
data-derived estimates. The report documents how window selection exposes the
weak component and why fitting the original samples is needed for verification.
