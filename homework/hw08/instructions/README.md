# Supplied HW8 material

These reference files support NCU EE8091 HW8: transmitter FFE design and simulation.

- `HW8_FFE.pdf`: supplied assignment, including the five required equalizer cases
  and the optional ADS exercise.
- `HW8_FFE_sample.py`: incomplete course template, retained as reference material.
  Its blanks mean it is not an executable solution.
- `HW7_stripline_50cm.s2p`: unchanged 50 cm two-port stripline channel data reused
  from HW7. The option line specifies Hz, magnitude/angle, and 50 ohms.
- `snpreader.py`: supplied legacy Touchstone reader, retained unchanged for
  reference. The rewritten solver does not require this helper.
- `cautions.txt`: original usage note about changing precursor/postcursor counts
  in the legacy solution. The rewritten solver generates all five cases in one run.
- `manifest.json`: original byte sizes and SHA-256 checksums for all five files,
  with logical source labels that do not expose personal directory names.

The assignment and template were supplied together. The channel, helper, and
usage note were provided with the original solution. All five references are
preserved byte for byte; no student name or student ID was found in their text or assignment PDF
metadata, so no sanitization was needed. Existing attribution and any terms of
the supplied materials remain with their respective owners.

The original personally identifying solution and report are replaced by the
rewritten solver and anonymous report in the parent directory. Screen recordings,
Python caches, and unrelated local files are excluded. To check these references,
run `python homework/hw08/verify_inputs.py` from the repository root.
