"""Build the anonymous HW6 report from the checked-in Markdown source.

Run from any directory: python homework/hw06/report/build_report.py
The report contains no student name, number, personal date or machine path.
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
    canvas.drawRightString(PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 25, "HW06 / PEAK DISTORTION ANALYSIS")
    canvas.line(MARGIN, 33, PAGE_WIDTH - MARGIN, 33)
    canvas.drawString(MARGIN, 21, "Anonymous numerical report")
    canvas.drawRightString(PAGE_WIDTH - MARGIN, 21, f"{doc.page}")
    canvas.restoreState()


def build(source, output):
    rl_config.invariant = True
    document = BaseDocTemplate(
        str(output), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=51, bottomMargin=45, title="HW6 - Peak Distortion Analysis",
        author="", subject="NCU EE8091 numerical report",
        creator="Anonymous report build", pageCompression=1,
    )
    frame = Frame(MARGIN, 45, WIDTH, PAGE_HEIGHT - 96,
                  leftPadding=0, rightPadding=0, bottomPadding=0, topPadding=0)
    document.addPageTemplates(PageTemplate(id="Report", frames=frame, onPage=page))
    document.build(parse_markdown(source))
    print(f"Built {output.name}")



def report_markdown(summary, verification):
    """Refresh the six-page report using the verified default analysis."""
    s, v = summary, verification
    c, spectrum, pulse, pda = s['config'], s['spectrum'], s['sbr'], s['pda']
    exact, filtered = pda['exact'], pda['filtered']
    random, pattern = s['cases']['random'], s['cases']['pattern']
    if not (c['samples_per_ui'] == 250 and c['bit_period_ps'] == 25
            and c['rise_time_ps'] == 5 and c['vpp_v'] == 2
            and c['kaiser_alpha'] == 1 and c['truncation_ns'] == 10
            and c['tolerance_v'] == .0002 and c['random_seed'] == 0
            and v['status'] == 'passed'):
        raise ValueError('Refresh requires the documented verified default configuration')
    rows = '\n'.join(
        f"| {item['truncation_ns']:g} | {item['exact']['maximum_vertical_gap_v']:.9f} | "
        f"{item['filtered']['maximum_vertical_gap_v']:.9f} |"
        for item in s['truncation_study'])
    largest_error = max(case['maximum_extrema_difference_v'] for case in v['cases'].values())
    error_ceiling = math.ceil(largest_error / 1e-18) * 1e-18
    return f"""> NCU EE8091 / Signal Integrity and Electromagnetic Compatibility

# HW6: Peak distortion analysis

Anonymous technical report. Rewritten calculations, regenerated figures and independently checked numerical results.

## 1. Assignment and signal model

The supplied two-port network carries a 40 Gbit/s unipolar NRZ source. This work computes the worst-case inner eye boundaries, overlays them with random-data and supplied-pattern eyes, and explains the difference. All four supplied instruction files are retained unchanged in `instructions/`.

| Parameter | Value or convention |
| --- | --- |
| Input levels and UI | 0 to 2 V; Tb = 25 ps |
| Rise and fall time | Equal 5 ps linear ramps |
| Matched transfer | H(f) = S21(f) / 2; 50 ohm ports |
| Supplied S21 grid | {spectrum['frequency_points']:,} points; DC to 80 GHz; 10 MHz spacing |
| Time grid and IDFT | 0.1 ps; 250 samples/UI; {spectrum['idft_points']:,} points |
| Kernel and blocks | 10 ns inclusive SBR; {c['sbr_samples']:,} taps; 1,024 bits per section |
| PDA tolerance | 0.2 mV per interfering tap; unfiltered bounds also retained |

The single-bit input is a causal trapezoid: it rises on [0, Tr], stays high on [Tr, Tb], and falls on [Tb, Tb + Tr]. With sinc(x) = sin(pi*x)/(pi*x), its physical Fourier transform is:

```text
Vin(f) = Vpp * Tb * sinc(f*Tb) * sinc(f*Tr)
         * exp(-j*pi*f*(Tb + Tr))
P(f)   = Vin(f) * S21(f)/2 * W(f)
```

Here W(f) = I0(pi*alpha*sqrt(1-(f/fmax)^2))/I0(pi*alpha), with alpha = 1 and fmax = 80 GHz. The causal phase uses the sum Tb + Tr. Frequencies beyond the supplied band are set to zero. The taper and finite band are part of the numerical channel model, rather than measured extensions of it.

The real IDFT uses Hermitian symmetry and the scale M*df. DC is counted once and positive/negative frequency pairs are counted once each. No additional factor of two or matched-source factor is applied during the later convolution.

**Recorded result.** The finite-model worst-case vertical opening is **{exact['maximum_vertical_gap_v']:.9f} V**. The supplied pattern approaches it much more closely than the equally long random record. These are deterministic ISI results; no BER, jitter or random-noise performance is inferred.

<!-- pagebreak -->

## 2. SBR and worst-case boundaries

![Figure 1. Reconstructed SBR and logarithmic tail. The 10 ns truncation is shown separately from the sharp main response.](../results/sbr.png)

The SBR peak is {pulse['peak_voltage_v']:.9f} V at {pulse['peak_time_ps']:.1f} ps. Define kmain = {pulse['peak_index']} and kstart = kmain - 250 = {pulse['pda_start_index']}. PDA evaluates the 500 points beginning at kstart; the SBR peak therefore lies at exactly 1 UI on the plotted axis.

At each time t, write the output as the reference bit's response plus all other symbol-spaced SBR contributions. The smallest possible one includes every negative interfering contribution; the largest possible zero includes every positive one:

```text
WC1(t) = p(t) + sum of all negative nonmain cursor values
WC0(t) =        sum of all positive nonmain cursor values
gap(t) = WC1(t) - WC0(t)
```

![Figure 2. Solid curves use every retained SBR tap. Dashed curves omit interfering taps with absolute value below 0.2 mV.](../results/pda.png)

The exact maximum gap occurs at {exact['optimum_relative_time_ui']:.3f} UI, corresponding to SBR time {exact['optimum_sbr_time_ps']:.1f} ps and nominal global phase {exact['optimum_global_phase_ps']:.1f} ps. WC1 = {exact['minimum_one_v']:.9f} V and WC0 = {exact['maximum_zero_v']:.9f} V. "Exact" here refers to the finite, windowed, bandwidth-limited SBR used by the solver; it is not an infinite-support physical-channel guarantee.

<!-- pagebreak -->

## 3. Purely random NRZ

The random case uses NumPy's default_rng with seed 0 and the same {random['input_bits']:,}-bit length as the supplied file. The equal-ramp input can be written as a sum of shifted single-bit trapezoids. Linearity therefore permits sectional convolution of the binary impulse train with the already shaped channel SBR.

Every FFT section is padded to at least its input length plus the SBR length minus one. Entire output tails are added at their original sample offsets. The final 290-bit partial section is retained. This produces the same finite-model output as direct linear convolution, without applying the pulse shaping twice.

![Figure 3. Random NRZ: single-color trace view and full-record statistical eye, both overlaid with unfiltered PDA bounds.](../results/random_eye.png)

The waveform is folded into one UI using the sample origin kstart, then that one-UI plot is duplicated. No arbitrary time offset is fitted. For steady-state eyes, {random['guard_bits_per_end']} reference bits at each record end are excluded; {random['stable_bits']:,} stable windows remain. The solver still convolves all {random['input_bits']:,} input bits and their complete finite-model tail.

The single-color view shows {c['max_displayed_traces']} uniformly selected stable traces for readability. The statistical eye and every reported extremum use all stable windows. With {c['voltage_bins']} voltage bins from -0.1 to 1.1 V, each time column contains {random['histogram_samples_per_phase']:,} observations. Color is probability per voltage bin conditional on time, not a voltage density or BER. Zero-count bins are masked.

| Metric at the exact PDA optimum | Random record |
| --- | --- |
| Smallest observed one | {random['observed_min_one_at_exact_optimum_v']:.9f} V |
| Largest observed zero | {random['observed_max_zero_at_exact_optimum_v']:.9f} V |
| Observed vertical gap | {random['observed_gap_at_exact_optimum_v']:.9f} V |
| Distance above WC1 / below WC0 | {random['one_clearance_at_exact_optimum_v']*1e3:.6f} / {random['zero_clearance_at_exact_optimum_v']*1e3:.6f} mV |
| Maximum exact-bound violation | 0 V within numerical precision |

The random eye is wider than the worst-case prediction because this finite record need not contain the particular bit combinations that maximize every ISI contribution at a chosen phase.

<!-- pagebreak -->

## 4. Supplied pattern and comparison

The supplied bit file contains {pattern['input_bits']:,} bits, of which {pattern['input_ones']:,} are ones. It is processed with the identical SBR, block size, guard, time origin, voltage bins and plotting scales used for random NRZ.

![Figure 4. Supplied pattern: sampled single-color traces and the statistical eye using all stable windows. Solid overlays include every retained SBR tap.](../results/pattern_eye.png)

| At 0.976 UI | PDA bound | Supplied pattern |
| --- | --- | --- |
| Minimum one (V) | {exact['minimum_one_v']:.9f} | {pattern['observed_min_one_at_exact_optimum_v']:.9f} |
| Maximum zero (V) | {exact['maximum_zero_v']:.9f} | {pattern['observed_max_zero_at_exact_optimum_v']:.9f} |
| Vertical gap (V) | {exact['maximum_vertical_gap_v']:.9f} | {pattern['observed_gap_at_exact_optimum_v']:.9f} |

![Figure 5. Full-data statistical eyes on shared time, voltage and probability scales: random NRZ at left, supplied pattern at right.](../results/eye_comparison.png)

At this optimum, the supplied pattern is only {pattern['one_clearance_at_exact_optimum_v']*1e6:.3f} microvolts above WC1 and {pattern['zero_clearance_at_exact_optimum_v']*1e6:.3f} microvolts below WC0. Its observed gap is {pattern['observed_gap_at_exact_optimum_v']-exact['maximum_vertical_gap_v']:.9f} V larger than the bound, whereas the random gap exceeds it by {random['observed_gap_at_exact_optimum_v']-exact['maximum_vertical_gap_v']:.9f} V. This supports the assignment's intended distinction: a deliberately adverse pattern reveals rare ISI combinations more efficiently than an equally long random record. It does not establish that every boundary point is reached by one common bit pattern.

<!-- pagebreak -->

## 5. Tolerance, truncation and verification

A tolerance removes small taps; it does not make them random or physically absent. The filtered maximum gap is {filtered['maximum_vertical_gap_v']:.9f} V at {filtered['optimum_relative_time_ui']:.3f} UI. Ignoring many tiny contributions can cumulatively move an envelope by more than a single 0.2 mV threshold. The largest phasewise sum of omitted absolute contributions is {pda['maximum_sum_of_omitted_taps_v']*1e3:.6f} mV. Solid overlays therefore retain every tap, and tolerance checks use the accumulated omitted positive/negative sums.

| SBR truncation (ns) | Unfiltered maximum gap (V) | Filtered maximum gap (V) |
| --- | --- | --- |
{rows}

Extending the finite kernel from 10 to 20 ns changes the maximum gap by 0.218182 mV and either inner boundary by at most 0.218487 mV. The largest interior tail sample after 10 ns and before 90 ns is {pulse['maximum_interior_tail_after_truncation_v']*1e6:.6f} microvolts. The 100 ns periodic wraparound is excluded from this tail measurement. These sensitivity checks quantify finite-model dependence; they do not prove convergence to infinite support.

### Independent numerical checks

The saved-result verifier imports no functions from the analysis module. It rebuilds the full real SBR with a two-sided Hermitian IFFT and separately evaluates six time points by a direct frequency sum. It computes PDA using phase totals with the main contribution removed. To validate every eye sample without reusing overlap-add, it reconstructs each waveform phase as an independent symbol-rate FIR convolution.

| Check | Result |
| --- | --- |
| Synthetic unit tests | 9 passed; exhaustive bit combinations, tolerance signs, pulse gain and phase, direct convolution, temporal alignment, histogram and readers |
| Supplied-file integrity | All {v['instruction_hashes_verified']} instruction SHA-256 hashes match |
| PDA and truncation study | All {v['polyphase_pda_phases_verified']} conditional phases and all four truncation lengths agree |
| Full-record conditional values | {v['cases']['random']['waveform_samples_verified']:,} evaluations per case across two UI |
| Independent extrema difference | At most {error_ceiling:.3e} V |
| Histogram counts | Every bin agrees; {random['stable_bits']*c['samples_per_ui']:,} one-UI samples per case |
| Column normalization | Maximum error below 1.4e-15 |

The two-UI conditional evaluations overlap in time and retain the same reference-bit label; they are not distinct histogram observations. The displayed statistical eye is instead the duplicated one-UI histogram, as requested by the assignment.

<!-- pagebreak -->

## 6. Reproduce and interpret the results

Run these commands from the repository root with Python 3.10 or newer:

```text
python -m pip install -r homework/hw06/requirements.txt
python homework/hw06/src/hw6.py
python -m unittest discover -s homework/hw06/tests -v
python homework/hw06/tests/verify_results.py
```

To regenerate this report after the verified default run:

```text
python -m pip install -r homework/hw06/report/requirements.txt
python homework/hw06/report/build_report.py --refresh
```

The solver is portable: its default paths are resolved relative to the homework directory. It runs without a graphical desktop. `results/summary.json` records configuration and measured quantities; `results/arrays.npz` stores the SBR, boundaries, both input records, exact histogram counts, probabilities, trace subset and observed extrema. `results/verification.json` records the independent checks. `source_manifest.json` identifies all supplied instruction files by checksum.

### Interpretation and scope

- PDA directly predicts the deterministic finite-model inner eye limits from one SBR. A long random simulation is not needed to obtain those limits.
- The supplied pattern approaches the limits closely; the finite random record leaves a larger observed eye opening. Probabilities inside the eye still require a statistical calculation.
- Truncation controls the available ISI memory and computational cost. A per-tap tolerance must be reported separately from accumulated envelope relaxation.
- The one-UI folding origin follows kmain and kstart. Duplicating one UI samples the same time-phase distribution on both displayed UIs; directly grouping disjoint two-UI blocks can change finite-record pattern coverage.
- Linearity and equal rise/fall shaping are assumed. Random noise, jitter, nonlinear driver/receiver behavior and an explicit receiver threshold/error model are absent.
- The optional ADS task was not run. The validation here is independent numerical reconstruction, not comparison with a newly executed ADS simulation.

### Provenance and anonymity

The report and implementation are rewritten from the assignment and reviewed reference work, with new figures and verified calculations. The four original instruction files are preserved byte for byte. The channel and bit file mentioned separately are identical to the instruction copies, so they are stored once. The rewritten two-port reader is self-contained; the legacy `snpreader.py` is not required.

No original student report or personally identifying source file is included. This report has no personal author field, student number, email address, signature, or workstation path. Source files in the local reference folders are unchanged. The GitHub hosting account remains visible independently of the report's anonymous content.

**Reference.** NCU EE8091, HW6: Peak distortion analysis, including the circuit, required comparisons and implementation notes in `instructions/HW6_PDA.pdf`. The incomplete supplied starter remains in `instructions/HW6_PDA_sample.py`; the executable replacement is `src/hw6.py`.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true', help='Refresh report text from verified default results')
    args = parser.parse_args()
    source, output = HERE / 'report.md', HERE / 'hw06_report.pdf'
    if args.refresh or not source.exists():
        summary = json.loads((HERE.parent / 'results/summary.json').read_text(encoding='utf-8'))
        verification = json.loads((HERE.parent / 'results/verification.json').read_text(encoding='utf-8'))
        source.write_text(report_markdown(summary, verification), encoding='utf-8')
    build(source, output)


if __name__ == '__main__':
    main()
