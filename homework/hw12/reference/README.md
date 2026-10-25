# Supplied reference implementation

`legacy_reference.py.txt` retains the supplied completed script for comparison.
Its filename is anonymous, and its machine-specific channel-data path was
replaced with a generic filename. It is reference text, not an executable or
verified solution. Use `../src/hw12.py` for the corrected implementation.

The original reports containing personal identifiers are not redistributed.
The new anonymous report in `../report/` is based on the regenerated results.
All three assignment files are preserved byte for byte in `../instructions/`.
Their SHA-256 checksums are recorded in the input manifest.

The rewrite addresses the reference script's hard-coded path, missing external
reader, incorrect first-edge boundary handling, padded zero tail, probability
floor in analysis arrays, inconsistent time-axis spacing, fragile one-UI
slicing, and memory use from retaining multiple full-length waveforms.
