"""Build the anonymous HW3 report from the checked-in Markdown source.

Run from any directory: python homework/hw03/report/build_report.py
The report contains no student name, number, source date or machine path.
"""

from __future__ import annotations

import argparse
import html
import json
import math
from pathlib import Path
import re

from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.lib.pagesizes import A4
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, KeepTogether, PageBreak, PageTemplate,
    Paragraph, Preformatted, Spacer, Table, TableStyle,
)

HERE = Path(__file__).resolve().parent
INK = colors.HexColor("#12263A")
ACCENT = colors.HexColor("#006F7C")
MUTED = colors.HexColor("#516173")
PALE = colors.HexColor("#ECF4F5")
LINE = colors.HexColor("#D4DFE6")
PAGE_WIDTH, PAGE_HEIGHT = A4
MARGIN = 45
WIDTH = PAGE_WIDTH - 2 * MARGIN


def styles():
    specs = {
        "Body": dict(fontName="Helvetica", fontSize=10.0, leading=14.0,
                     textColor=INK, spaceAfter=7),
        "Title": dict(fontName="Helvetica-Bold", fontSize=27, leading=31,
                      textColor=INK, spaceBefore=3, spaceAfter=13),
        "H1": dict(fontName="Helvetica-Bold", fontSize=16.5, leading=21,
                   textColor=INK, spaceBefore=10, spaceAfter=9, keepWithNext=True),
        "H2": dict(fontName="Helvetica-Bold", fontSize=11.5, leading=16,
                   textColor=ACCENT, spaceBefore=9, spaceAfter=6, keepWithNext=True),
        "Caption": dict(fontName="Helvetica", fontSize=8.5, leading=11.5,
                        textColor=MUTED, spaceBefore=4, spaceAfter=9),
        "TableBody": dict(fontName="Helvetica", fontSize=8.6, leading=11.3,
                          textColor=INK),
        "TableHead": dict(fontName="Helvetica-Bold", fontSize=8.6, leading=11.3,
                          textColor=colors.white),
        "Kicker": dict(fontName="Helvetica-Bold", fontSize=9, leading=13,
                       textColor=ACCENT, spaceAfter=6),
        "Bullet": dict(fontName="Helvetica", fontSize=9.5, leading=13.2,
                       textColor=INK, leftIndent=10, firstLineIndent=-8,
                       spaceAfter=5),
        "Code": dict(fontName="Courier", fontSize=7.5, leading=10.3,
                     backColor=PALE, borderPadding=9, spaceBefore=4, spaceAfter=10),
    }
    return {key: ParagraphStyle(key, alignment=TA_LEFT, **values)
            for key, values in specs.items()}


STYLE = styles()


def rich(text):
    """Convert the deliberately small Markdown subset to ReportLab markup."""
    text = html.escape(text, quote=False)
    text = re.sub(r"`([^`]+)`", r'<font name="Courier" size="8.7">\1</font>', text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'\1', text)
    return text


def p(text, style="Body"):
    return Paragraph(rich(text), STYLE[style])


def table(rows):
    ncol = len(rows[0])
    # The parameter tables deliberately leave more room for their descriptions.
    widths = ([WIDTH * 0.31, WIDTH * 0.69] if ncol == 2 else
              [WIDTH * 0.28] + [WIDTH * 0.72 / (ncol - 1)] * (ncol - 1))
    cells = [[p(cell, "TableHead" if i == 0 else "TableBody")
              for cell in row] for i, row in enumerate(rows)]
    tbl = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LINEBELOW", (0, -1), (-1, -1), 0.5, LINE),
    ]))
    return tbl


def figure(path, caption):
    path = (HERE / path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    w, h = ImageReader(str(path)).getSize()
    image = Image(str(path), width=WIDTH, height=WIDTH * h / w)
    image.hAlign = "LEFT"
    return KeepTogether([image, p(caption, "Caption")])


def parse_markdown(source):
    lines = source.read_text(encoding="utf-8").splitlines()
    story = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line == "<!-- pagebreak -->":
            story.append(PageBreak())
            i += 1
            continue
        if line.startswith("<!--"):
            i += 1
            continue
        if line.startswith("```"):
            code = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code.append(lines[i])
                i += 1
            story.append(Preformatted("\n".join(code), STYLE["Code"]))
            i += 1
            continue
        match = re.match(r"!\[(.*)\]\(([^)]+)\)$", line)
        if match:
            story.append(figure(match.group(2), match.group(1)))
            i += 1
            continue
        if line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = [s.strip() for s in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", cell) for cell in row):
                    rows.append(row)
                i += 1
            story.extend([table(rows), Spacer(1, 9)])
            continue
        if line.startswith("# "):
            story.append(p(line[2:], "Title"))
        elif line.startswith("## "):
            story.append(p(line[3:], "H1"))
        elif line.startswith("### "):
            story.append(p(line[4:], "H2"))
        elif line.startswith("> "):
            story.append(p(line[2:], "Kicker"))
        elif line.startswith("- "):
            story.append(p("- " + line[2:], "Bullet"))
        else:
            # Wrapped Markdown paragraphs become one PDF paragraph.
            chunk = [line]
            while (i + 1 < len(lines) and lines[i + 1].strip()
                   and not re.match(r"^(#|\||!\[|```|<!--|- |> )", lines[i + 1].strip())):
                i += 1
                chunk.append(lines[i].strip())
            story.append(p(" ".join(chunk)))
        i += 1
    return story


def page(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(LINE)
    canvas.setLineWidth(0.6)
    canvas.line(MARGIN, PAGE_HEIGHT - 36, PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 36)
    canvas.setFont("Helvetica-Bold", 8)
    canvas.setFillColor(ACCENT)
    canvas.drawString(MARGIN, PAGE_HEIGHT - 25, "NATIONAL CENTRAL UNIVERSITY / EE8091")
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawRightString(PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 25, "HW03 / SECTIONAL CONVOLUTION")
    canvas.line(MARGIN, 33, PAGE_WIDTH - MARGIN, 33)
    canvas.drawString(MARGIN, 21, "Anonymous numerical report")
    canvas.drawRightString(PAGE_WIDTH - MARGIN, 21, f"{doc.page}")
    canvas.restoreState()


def build(source, output):
    rl_config.invariant = True
    document = BaseDocTemplate(
        str(output), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=51, bottomMargin=45, title="HW3 - Sectional Convolution",
        author="Anonymous", subject="NCU EE8091 numerical report",
        creator="Anonymous report build", pageCompression=1,
    )
    frame = Frame(MARGIN, 45, WIDTH, PAGE_HEIGHT - 96,
                  leftPadding=0, rightPadding=0, bottomPadding=0, topPadding=0)
    document.addPageTemplates(PageTemplate(id="Report", frames=frame, onPage=page))
    document.build(parse_markdown(source))
    print(f"Built {output.name}")


def report_markdown(summary):
    """Refresh the narrative and recorded values from the numerical run."""
    s = summary
    q, h, sim = s["parameters"], s["impulse"], s["simulation"]
    data, ads = s["input"], s["ads_comparison"]
    check, clock, env = s["direct_convolution_check"], s["timing"], s["environment"]
    windows = s["comparison_windows"]
    # Figure captions and the interpretation describe this documented baseline.
    # Refuse a misleading automatic refresh after changing analysis defaults.
    if not (data["bits"] == 2**20 and sim["block_bits"] == 512
            and math.isclose(q["dt_s"], 2e-12, rel_tol=1e-12, abs_tol=0)
            and math.isclose(h["tail_duration_s"], 40e-9, rel_tol=1e-12, abs_tol=0)
            and math.isclose(h["kaiser_alpha"], 0.5, rel_tol=1e-12, abs_tol=0)):
        raise ValueError("This report narrative describes the documented default configuration; update the narrative before refreshing a different configuration")
    us = lambda value: f"{value * 1e6:.6f}"
    mv = lambda value: f"{value * 1e3:.4f}"
    number = lambda value: f"{value:,}"
    rows = []
    for i, window in enumerate(windows, 1):
        rows.append(
            f"| {i}: {us(window['start_s'])}-{us(window['end_s'])} | "
            f"{number(window['sample_count'])} | {mv(window['rmse_v'])} | "
            f"{mv(window['max_absolute_error_v'])} |"
        )
    table_a = "\n".join(rows[:2])
    table_b = "\n".join(rows[2:])
    omitted = ads["excluded_sample_count"]
    ref_end = ads["file_end_s"]
    window_note = ("No additional taper is applied because alpha is zero; the measured band is still truncated at 40 GHz."
                   if h["kaiser_alpha"] == 0 else
                   "The moderate taper follows the assignment's guidance that the finite input rise time already suppresses high-frequency content.")
    study = s.get("verification_study")
    if study is None:
        raise ValueError("Regenerate analysis with --verification-study before refreshing the report")
    variant_names = {
        "half_sample_interval": "1 ps sampling",
        "no_kaiser_taper": "alpha = 0",
        "weak_kaiser_taper": "alpha = 0.1",
        "shorter_tail": "T2 = 20 ns",
        "longer_tail": "T2 = 60 ns",
    }
    study_rows = "\n".join(
        f"| {variant_names[v['variant']]} | {mv(v['ads_rmse_v'])} | "
        f"{mv(v['ads_max_absolute_error_v'])} | {mv(v['max_change_from_baseline_v'])} |"
        for v in study["variants"]
    )
    return f"""> National Central University | EE8091 | Signal Integrity and Electromagnetic Compatibility

# HW3: Sectional convolution

**Anonymous numerical report**

This study computes the receiver waveform of a 20 cm microstrip channel using its supplied S21 data and FFT overlap-add. All {number(data['bits'])} bits are processed at {q['dt_s'] * 1e12:g} ps spacing. The known input duration is **{us(data['duration_s'])} us**, rather than a rounded 105 us.

On the known input interval, comparison at {number(ads['sample_count'])} ADS timestamps gives an RMSE of **{mv(ads['rmse_v'])} mV** and a maximum absolute difference of **{mv(ads['max_absolute_error_v'])} mV**. The measured overlap-add loop takes **{sim['convolution_seconds']:.3f} s** for the complete pattern. These are results from this run; no shift, gain or voltage-offset fit is used.

## 1. Problem and numerical settings

| Quantity | Value and interpretation |
| --- | --- |
| Course / task | NCU EE8091 / time-frequency transform (3): sectional convolution |
| Channel | {q['line_length_m'] * 100:g} cm microstrip; S21 used with matched source and load |
| Source / load | {q['source_ohm']:g} ohm / {q['load_ohm']:g} ohm; S-parameter reference is 50 ohm |
| Bit rate / UI | {q['bit_rate_bps'] / 1e9:g} Gbps / {q['ui_s'] * 1e12:g} ps |
| Input levels | {q['incident_low_v']:g} to {q['incident_high_v']:g} V incident wave; source high level is 2 V before the matched divider |
| Rise and fall | {q['rise_fall_s'] * 1e12:g} ps, 0-100% linear transition starting at each changed bit boundary |
| Input length | 2^20 = {number(data['bits'])} bits; {number(data['ones'])} ones and {number(data['zeros'])} zeros |
| Input duration | {us(data['duration_s'])} us |
| Sample interval / rate | {q['dt_s'] * 1e12:g} ps / {1 / q['dt_s'] / 1e9:g} Gsample/s |
| Points per UI / ramp | {sim['samples_per_ui']} / {sim['rise_points']}; new level is reached after {sim['rise_points']} sample intervals |
| Bits / sections | {sim['block_bits']} bits per full section; {number(sim['block_count'])} sections; final section contains {sim['final_block_bits']} bits |
| Output storage | float32 disk-backed array, {number(sim['output_points'])} samples ({sim['output_bytes'] / 1e6:.3f} MB decimal) |

All physical signal settings come from the assignment. Sampling, section size, window strength and impulse cutoff are explicit implementation choices. Sections need not divide the bit count: the loop also handles a shorter final section.

<!-- pagebreak -->

## 2. Channel impulse response and cutoff

Because both terminations match the 50 ohm S-parameter reference, the incident-to-received voltage transfer function is S21. The real time-domain response is obtained by a Hermitian inverse real FFT. The transform supplies conjugate negative frequencies and counts DC once, equivalent to the lecture's twice-real one-sided construction with half-weight DC.

| Transform parameter | Value |
| --- | --- |
| Measured band / spacing | DC to {h['frequency_max_hz'] / 1e9:g} GHz / {h['frequency_step_hz'] / 1e6:g} MHz; {number(h['frequency_samples'])} frequency samples |
| Period from spacing | 1/df = {h['impulse_period_s'] * 1e9:g} ns; {number(h['impulse_period_points'])} time samples |
| Frequency interpolation | None within the supplied uniform grid; unmeasured frequencies above 40 GHz are zero-filled |
| Kaiser window | alpha = {h['kaiser_alpha']:g}; beta = pi*alpha = {h['kaiser_beta']:.6f}; W(0) = 1 |
| Retained FIR | 0 to {h['tail_duration_s'] * 1e9:g} ns inclusive; {number(h['fir_points'])} coefficients |
| Impulse peak | Delay {h['peak_delay_s'] * 1e9:.3f} ns; continuous amplitude {h['peak_continuous_impulse_per_s']:.6e} 1/s |

The window is W(f) = I0[pi*alpha*sqrt(1-(f/fmax)^2)] / I0(pi*alpha). {window_note} The discrete coefficients are `h_d[n] = dt*h_cont[n]`; convolving voltage samples with `h_d` needs no additional `dt` multiplier.

![Figure 1. The normalized impulse over one 100 ns period. The late rise near the period end is the wrapped negative-time part of the band-limited transform.](../results/impulse_response.png)

The selected cutoff is T2 = {h['tail_duration_s'] * 1e9:g} ns. On the explicitly checked positive-time interval {h['positive_tail_check_start_s'] * 1e9:.3f}-{h['positive_tail_check_end_s'] * 1e9:g} ns, the largest discarded magnitude is **{h['positive_tail_peak_relative_db']:.2f} dB** relative to the impulse peak, below the -80 dB guideline. This interval check is not an infinite-time bound and does not cover the wrapped part near 100 ns.

The supplied DC gain is {h['channel_dc_gain']:.9f}; the retained FIR sums to {h['retained_discrete_dc_gain']:.9f}, a difference of {h['dc_gain_difference']:.9f}. No DC renormalization hides this effect. Finite bandwidth, windowing and causal truncation remain approximations; a small individual tail coefficient does not guarantee a 0.1 mV bound on the accumulated waveform error.

<!-- pagebreak -->

## 3. Input waveform and overlap-add

At a bit change, the old level is retained at the boundary sample and then moves linearly to the new level over 20 ps. Equal adjacent bits retain a constant level. The previous bit of each section is taken from the full pattern; only the first section assumes a preceding zero. This preserves a physical ramp across a section boundary.

![Figure 2. First 20 input bits after construction. The waveform uses the 0-1 V incident levels, not the 0-2 V source levels.](../results/input_waveform.png)

For an input section of L samples and an FIR of M samples, linear convolution has L+M-1 samples. The FFT is zero-padded to at least that length. With the default choices, L = {sim['block_bits'] * sim['samples_per_ui']:,}, M = {h['fir_points']:,}, L+M-1 = {sim['maximum_linear_output_points_per_block']:,}, and NFFT = {sim['fft_points']:,}. The FIR spectrum is computed once and reused.

```text
H = rfft(h_discrete, n=NFFT)
for start_bit in range(0, number_of_bits, block_bits):
    take all remaining bits up to block_bits
    construct x using the real preceding bit
    z = irfft(rfft(x, n=NFFT) * H, n=NFFT)
    keep len(x) + len(h_discrete) - 1 points
    add z at output offset start_bit * samples_per_ui
```

The section input intervals are adjacent. Their **output tails overlap and are added** at the original sample offsets. This addition is the step that makes sectional convolution equivalent to whole-waveform linear convolution for the chosen FIR. FFT arithmetic uses NumPy's float64 transform; the saved voltage array and each update are float32. No full input waveform or full time vector is allocated.

### Coverage and implementation corrections

- Every supplied bit is consumed. Floor-dividing the original 100-bit section size would process 1,048,500 bits and drop the last 76; the rewritten loop uses start indices through the entire pattern.
- FFT padding uses the input **sample** count plus the FIR count minus one, avoiding circular wraparound and confusion between bits and samples.
- A preceding bit is carried across sections, so the first transition in a section has the correct ramp.
- Time is derived from integer sample indices with float64 arithmetic when comparing short windows, avoiding a long float32 time vector.
- The Touchstone reader selects S21 in the standard two-port order and supports RI, MA and DB complex formats. File locations are relative to the program or provided explicitly.

<!-- pagebreak -->

## 4. ADS comparison: startup and early waveform

The ADS export has nonuniform timestamps, with steps ranging from {ads['min_ads_step_s'] * 1e12:.6g} to {ads['max_ads_step_s'] * 1e12:.6g} ps. The computed voltage is linearly interpolated at each original ADS timestamp. Error is `overlap_add - ADS`. No time shift, amplitude scale or constant offset is fitted. RMSE is sample-weighted on these timestamps, not a time-integrated error norm.

![Figure 3. Startup comparison, 0-10 ns. The lower axis shows the signed voltage difference.](../results/comparison_1.png)

![Figure 4. Early-data comparison, 8.000000-8.010000 us.](../results/comparison_2.png)

| Window (us) | ADS samples | RMSE (mV) | Max abs. (mV) |
| --- | --- | --- | --- |
{table_a}

<!-- pagebreak -->

## 5. ADS comparison: middle and final known bits

![Figure 5. Mid-pattern comparison, 50.000000-50.010000 us.](../results/comparison_3.png)

![Figure 6. Final 10 ns of the supplied pattern, 104.847600-104.857600 us. The right endpoint is excluded from the metrics.](../results/comparison_4.png)

| Window (us) | ADS samples | RMSE (mV) | Max abs. (mV) |
| --- | --- | --- | --- |
{table_b}

Only timestamps before {us(data['duration_s'])} us enter the metrics. The later ADS continuation is outside the specified input; its {number(omitted)} rows are excluded. Section 7 explains this boundary.

<!-- pagebreak -->

## 6. Verification, runtime and reproducibility

### Numerical checks and overall accuracy

| Check or metric | Recorded result |
| --- | --- |
| Input coverage | {number(sim['processed_bit_count'])} of {number(data['bits'])} bits processed; {number(sim['block_count'])} sections |
| FFT padding | {sim['fft_points']:,} >= {sim['maximum_linear_output_points_per_block']:,}; full linear output retained per section |
| Independent convolution check | {check['bit_count']} bits in {check['block_bits']}-bit sections; final partial section {check['final_block_bits']} bits; FFT overlap-add compared with `numpy.convolve` of the whole input |
| Direct-check maximum error | {check['max_absolute_error_v']:.9e} V; passes tolerance {check['tolerance_v']:.1e} V |
| ADS samples compared | {number(ads['sample_count'])}; timestamp range {us(ads['start_s'])}-{us(ads['end_s'])} us |
| Overall RMSE / mean abs. error | {mv(ads['rmse_v'])} / {mv(ads['mean_absolute_error_v'])} mV |
| Overall max. abs. / mean signed | {mv(ads['max_absolute_error_v'])} / {mv(ads['mean_signed_error_v'])} mV |

The independent test exercises overlap addition, state transfer and a final partial section. It checks the FFT algorithm against a separate direct convolution. **Algorithmic correctness and ADS agreement are separate findings:** the algorithm passes its numerical check while the ADS residual above remains unresolved.

### Measured computation time

| Item | Time or environment |
| --- | --- |
| Overlap-add loop | {sim['convolution_seconds']:.3f} s, including section construction and output flush |
| ADS file streaming / comparison | {clock['ads_streaming_seconds']:.3f} s |
| Complete analysis run | {clock['total_seconds']:.3f} s, including figures and output processing measured before summary serialization |
| Software / platform | Python {env['python']}; NumPy {env['numpy']}; {env['platform']} {env['processor']} |

The overlap-add timer starts after output allocation and the FIR FFT. These times describe one execution on the available machine and storage system. ADS was not rerun here, so no measured speedup over ADS is claimed. The assignment's example PC timings are not substituted for a local benchmark.

<!-- pagebreak -->

## 7. Sensitivity, scope and reproduction

### Sampling, window and tail sensitivity

This study recomputes only the four 10 ns windows ({number(study['reference_samples'])} ADS timestamps), preserving all preceding input history required by each finite FIR. It does not rerun the full pattern for each variant. The baseline on this subset is {mv(study['baseline_ads_rmse_v'])} mV RMSE and {mv(study['baseline_ads_max_absolute_error_v'])} mV maximum absolute error; subset values differ from the whole-interval metrics.

| Variant | ADS RMSE (mV) | ADS max. (mV) | Max change (mV) |
| --- | --- | --- | --- |
{study_rows}

"Max change" compares the variant's voltage with the baseline at the same ADS timestamps. Halving the sample interval changes the waveform by at most 0.291 mV; extending T2 to 60 ns changes it by at most 0.096 mV. Removing the taper reduces the subset RMSE, but leaves a residual above 20 mV. These checks do not identify a unique cause of the ADS mismatch. The reported default remains alpha = 0.5, without fitted alignment or voltage corrections.

### Known input and reference boundary

ADS continues to {us(ref_end)} us; the supplied pattern specifies only `0 <= t < {us(data['duration_s'])} us`. Exactly {number(omitted)} later reference rows are excluded. The finite convolution's output tail assumes zero input after the pattern ends. The continued bit sequence used by ADS is unspecified, so that tail is not validated against the later ADS continuation.

### Limits and rerunning the work

The 40 GHz frequency limit, Kaiser weighting, 2 ps waveform sampling, finite causal FIR and float32 storage all contribute approximation error. The 40 ns tail test checks one finite interval; it is not an infinite-support proof. Comparison against the provided ADS waveform is the evidence for the reported configuration, rather than a claim that alternative samplings or cutoffs produce the same error.

From the repository root, run the following after installing `requirements.txt`; PDF generation additionally uses `homework/hw03/report/requirements.txt`. The numerical input hash verifier and automated solver tests are supplied beside the analysis.

```text
python homework/hw03/verify_inputs.py
python homework/hw03/analyze.py --verification-study
python homework/hw03/report/build_report.py --refresh
```

`analyze.py` writes the JSON summary, four compact waveform CSVs and all figures. The large disk-backed waveform is deleted by default; `--keep-waveform` retains it. `--refresh` updates this Markdown narrative from the new summary before rebuilding the PDF. The complete executable source is included in `analyze.py`; the short algorithm above describes its central loop.

**References.** NCU EE8091 HW3 assignment, including Appendix 2 (input timing) and Appendix 3 (Kaiser window), in `instructions/`; supplied bit pattern, channel S2P and ADS waveform, identified by the input manifest. Numerical values in this report are recorded in `results/summary.json`.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE / "report.md")
    parser.add_argument("--output", type=Path, default=HERE / "report.pdf")
    parser.add_argument("--refresh", action="store_true", help="Refresh Markdown values from results/summary.json")
    parser.add_argument("--summary", type=Path, default=HERE.parent / "results" / "summary.json")
    args = parser.parse_args()
    if args.refresh:
        summary = json.loads(args.summary.read_text(encoding="utf-8"))
        args.source.write_text(report_markdown(summary), encoding="utf-8")
    build(args.source, args.output)


if __name__ == "__main__":
    main()
