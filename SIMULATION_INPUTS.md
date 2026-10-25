# Homework simulation input inventory

The numerical inputs below are included in the repository. They were checked
against the supplied homework folders; the accompanying
[SHA-256 inventory](simulation_inputs_manifest.json) records the stored files
and source comparison. This inventory covers homework inputs and reference
helpers, rather than bundled third-party simulator examples.

Course materials and the supplied Touchstone reader are credited to
**Prof. Chiu-Chih Chou (周求致), Department of Electrical Engineering, National
Central University, Taiwan**. See [attribution and original source links](ATTRIBUTION.md).

## Inputs by homework

| Homework | Supplied numerical or simulator inputs |
| --- | --- |
| [HW1](homework/hw01/) | [HW1_Kaiser.txt](homework/hw01/instructions/HW1_Kaiser.txt) |
| [HW2](homework/hw02/) | [SPICE deck](homework/hw02/instructions/HW2_Chebyshev_LPF_S11_sample.sp), [second supplied deck](homework/hw02/instructions/HW2_Chebyshev_LPF_S11_sample_1.sp), and [original reference S2P](homework/hw02/reference/part2_reference/HW2_Chebyshev_LPF.s2p); the rewritten solver also has its [generated S2P](homework/hw02/data/chebyshev_lpf.s2p) |
| [HW3](homework/hw03/) | [channel S2P](homework/hw03/instructions/HW3_mline_channel.s2p), [bit pattern](homework/hw03/instructions/HW3_bitpattern_1048576.txt), and [ADS reference waveform, lossless gzip](homework/hw03/instructions/HW3_Rx_waveform_ADS.txt.gz) |
| [HW4](homework/hw04/) | [bit pattern 1](homework/hw04/instructions/HW4_bit_pattern_1.txt), [pattern 2](homework/hw04/instructions/HW4_bit_pattern_2.txt), and [pattern 3](homework/hw04/instructions/HW4_bit_pattern_3.txt) |
| [HW5](homework/hw05/) | [HW5_ADS_waveforms.csv](homework/hw05/instructions/HW5_ADS_waveforms.csv) |
| [HW6](homework/hw06/) | [channel S2P](homework/hw06/instructions/HW6_channel.s2p) and [worst-case bit pattern](homework/hw06/instructions/HW6_wcbitpattern.txt) |
| [HW7](homework/hw07/) | [HW7_stripline_50cm.s2p](homework/hw07/instructions/HW7_stripline_50cm.s2p) |
| [HW8](homework/hw08/) | [HW7_stripline_50cm.s2p](homework/hw08/instructions/HW7_stripline_50cm.s2p), reused by the assignment |
| [HW9](homework/hw09/) | [HW7_stripline_50cm.s2p](homework/hw09/data/HW7_stripline_50cm.s2p), reused by the assignment |
| [HW10](homework/hw10/) | [HW10_complicated_diff_channel.s4p](homework/hw10/instructions/HW10_complicated_diff_channel.s4p) |
| [HW11](homework/hw11/) | [HW11_crosstalk_8lines.s16p](homework/hw11/instructions/HW11_crosstalk_8lines.s16p) |
| [HW12](homework/hw12/) | [HW12_channel.s2p](homework/hw12/instructions/HW12_channel.s2p) |
| [HW13](homework/hw13/) | [provided reference CSV](homework/hw13/inputs/jitter_variables_for_hw13.csv); the rewritten analysis uses the [verified HW12 output CSV](homework/hw13/inputs/hw12_verified_jitter.csv), with its provenance documented in the HW13 guide |
| [HW14](homework/hw14/) | [HW14_diff_channel_50cm.s4p](homework/hw14/instructions/HW14_diff_channel_50cm.s4p) and [assignment PDF](homework/hw14/instructions/HW14_NRZ_vs_PAM4.pdf) |

## Reader and file formats

The supplied helper's filename is **`snpreader.py`**, not `snpreader.sp`.
Touchstone `.s2p`, `.s4p`, and `.s16p` files hold 2-, 4-, and 16-port S parameters;
`.sp` files are SPICE circuit decks. These are different types of simulation assets.

Reader copies are preserved with the relevant reference packages, including
HW2, HW3, HW6, HW7, HW8, HW9, HW10, HW11 (`legacy_snpreader.py`), and HW12. The rewritten
solvers may use their own checked readers; retaining the original helper does
not mean an incomplete teaching template is the executable solution.

The HW3 ADS reference is stored as `.txt.gz` because the original text exceeds
GitHub's ordinary file-size limit. Decompression preserves the complete original
waveform; it is not a shortened sample.

HW2's [additional supplied references](homework/hw02/reference/auxiliary/)
also retain `raw.py`, `last.txt`, and the separate `prbs.txt` waveform. These are
kept for completeness; the verified HW2 solver uses the inputs documented in
its guide. Its cleaned/generated S2P preserves the source numerical records,
while the untouched source S2P is also retained under `reference/`.

## HW14 channel preservation

The supplied HW14 assignment PDF is unchanged. The channel copy removes only
export comments containing a local path and project identifiers; its option line
and all **264,033 numerical tokens** match the original. The
[HW14 input manifest](homework/hw14/instructions/input_manifest.json) records
both the original and published hashes and the numerical-token digest.

The assignment states that no sample code is supplied. The repository now also
includes the [rewritten HW14 solver, comparisons, checks, and anonymous report](homework/hw14/).
