# Supplied HW3 material

These reference files support NCU EE8091 HW3: sectional convolution.

- `HW3_time_freq_transf_(3)_sectional_convolution.pdf`: supplied assignment.
- `HW3_bitpattern_1048576.txt`: unchanged complete binary pattern.
- `HW3_mline_channel.s2p`: unchanged two-port channel data, 0 to 40 GHz.
- `HW3_Rx_waveform_ADS.txt.gz`: lossless compressed original ADS reference. The
  solver reads gzip directly. To restore the original text explicitly, run
  `python -m gzip -d HW3_Rx_waveform_ADS.txt.gz` on a copy.
- `HW3_sectional_convolution_sample.py`: incomplete course template, retained as
  reference material. Its blanks mean it is not an executable solution.
- `snpreader.py` and `raw.py`: supplied legacy helper scripts, retained unchanged
  for reference. The verified rewrite does not require these helpers.
- `manifest.json`: file sizes and SHA-256 checksums; the ADS entry contains both
  the stored gzip checksum and the decompressed original checksum.

The original personally identifying report and source header are replaced by the
anonymous report and rewritten solver in the parent directory. The two recordings
are excluded. Numerical reference data are unchanged. Attribution and any existing
terms of supplied materials remain with their respective owners.
