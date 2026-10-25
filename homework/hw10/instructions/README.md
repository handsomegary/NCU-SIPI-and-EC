# Supplied HW10 material

The following three files are complete, byte-for-byte copies of the supplied
assignment directory:

- `HW10_diff_signaling.pdf`: assignment and optional ADS setup.
- `HW10_complicated_diff_channel.s4p`: four-port channel measurements, including DC.
- `HW10_differential_signaling_sample.py`: original notebook-oriented template.

`snpreader.py` is the supplementary legacy reader supplied with the reference
solution. The rewritten solver has its own validated reader and does not depend
on this helper or execute the incomplete template.

`manifest.json` records the byte count and SHA-256 digest of every supplied file.
From the repository root, run `python homework/hw10/verify_inputs.py` to check them.

The channel file uses Hz, magnitude/angle values, 50-ohm single-ended reference
impedances, and the row-major port ordering stated in its header. Its differential
pairs are ports (1, 3) at the input and (2, 4) at the output.

Personally identifying source filenames and the original student report are
replaced by the rewritten code and newly generated anonymous report in this
homework directory. Course materials remain distinct from the verified solution.
