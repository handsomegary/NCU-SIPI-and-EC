# Legacy analysis references

These two text files preserve both supplied HW13 Python variants under neutral
filenames. They are reference material, not the executable analysis for this
delivery. The `.py.txt` suffix deliberately keeps them out of Python test and
compilation discovery.

- `legacy_archived_analysis.py.txt`: the Python payload from the provided archive.
- `legacy_adjacent_analysis.py.txt`: the separate Python file supplied beside it.

The source bytes are preserved; original identifying archive and file names are
not retained. Original student PDF/DOCX reports and the original archive are not
included. See `../instructions/manifest.json` for SHA-256 hashes.

## Known issues in the legacy variants

- Both assume the current working directory contains the CSV.
- Both retain an artificial `1e-15` floor inherited from the former statistical
  eye construction. That floor is not observed tail probability.
- Both use 200 endpoint-inclusive time coordinates spanning 25 ps while using
  `25/200` ps as the sample interval, creating a coordinate mismatch.
- Both compute extreme Gaussian tails through subtraction from one, which can
  cancel small probabilities in floating-point arithmetic.
- Both mutate the time-step variable when changing from seconds to picoseconds.
- Both contain a zero upper limit for a logarithmic plot axis.
- Their manually selected Gaussian parameters differ; neither performs an
  objective tail-fit quality check or computes numerical eye-width/minimum-BER
  metrics.

The current analysis uses the copied, verified HW12 distribution as its primary
input. The originally supplied CSV remains available as a legacy comparison.
No claim is made that the legacy CSV was regenerated from either archived code
variant.
