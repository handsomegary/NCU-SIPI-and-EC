"""Build an anonymous HW8 report from the recorded numerical results.

The Markdown source, PDF metadata and generated figures contain no student
identity. Run with --refresh after regenerating the default simulation results.
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
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, KeepTogether, PageBreak, PageTemplate,
    Paragraph, Preformatted, Spacer, Table, TableStyle,
)

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
PAGE_WIDTH, PAGE_HEIGHT = A4
MARGIN = 44
WIDTH = PAGE_WIDTH - 2 * MARGIN
INK = colors.HexColor("#162D40")
ACCENT = colors.HexColor("#007A87")
MUTED = colors.HexColor("#526575")
LINE = colors.HexColor("#D9E3E8")
PALE = colors.HexColor("#EEF5F6")


def make_styles():
    definitions = {
        "Body": dict(fontName="Helvetica", fontSize=10, leading=14,
                     textColor=INK, spaceAfter=8),
        "Title": dict(fontName="Helvetica-Bold", fontSize=27, leading=32,
                      textColor=INK, spaceAfter=15),
        "Heading": dict(fontName="Helvetica-Bold", fontSize=17, leading=22,
                        textColor=INK, spaceAfter=12, keepWithNext=True),
        "Subheading": dict(fontName="Helvetica-Bold", fontSize=11,
                           leading=15, textColor=ACCENT, spaceBefore=7,
                           spaceAfter=7, keepWithNext=True),
        "Caption": dict(fontName="Helvetica", fontSize=8.5, leading=11.5,
                        textColor=MUTED, spaceBefore=3, spaceAfter=12),
        "TableBody": dict(fontName="Helvetica", fontSize=8.4, leading=11,
                          textColor=INK),
        "TableHead": dict(fontName="Helvetica-Bold", fontSize=8.4,
                          leading=11, textColor=colors.white),
        "Kicker": dict(fontName="Helvetica-Bold", fontSize=9, leading=13,
                       textColor=ACCENT, spaceAfter=7),
        "Code": dict(fontName="Courier", fontSize=8, leading=11,
                     textColor=INK, backColor=PALE, borderPadding=9,
                     spaceBefore=6, spaceAfter=13),
    }
    return {k: ParagraphStyle(k, alignment=TA_LEFT, **v)
            for k, v in definitions.items()}


STYLE = make_styles()


def rich(text):
    text = html.escape(str(text), quote=False)
    text = re.sub(r"`([^`]+)`", r'<font name="Courier" size="8">\1</font>', text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    return re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", text)


def p(text, style="Body"):
    return Paragraph(rich(text), STYLE[style])


def table(rows):
    n = len(rows[0])
    if n == 2 and rows[0][0].startswith("Numerical"):
        widths = [WIDTH * 0.72, WIDTH * 0.28]
    elif n == 2:
        widths = [WIDTH * 0.35, WIDTH * 0.65]
    elif n == 3:
        widths = [WIDTH * 0.21, WIDTH * 0.39, WIDTH * 0.40]
    else:
        widths = [WIDTH / n] * n
    cells = [[p(cell, "TableHead" if i == 0 else "TableBody")
              for cell in row] for i, row in enumerate(rows)]
    item = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
    item.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LINEBELOW", (0, -1), (-1, -1), .6, LINE),
    ]))
    return item


def figure(relative_path, caption):
    path = (HERE / relative_path).resolve()
    width, height = ImageReader(str(path)).getSize()
    scale = min(WIDTH / width, 259 / height)
    img = Image(str(path), width=width * scale, height=height * scale)
    img.hAlign = "CENTER"
    return KeepTogether([img, p(caption, "Caption")])


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
        elif line.startswith("<!--"):
            pass
        elif line.startswith("```"):
            block = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                block.append(lines[i])
                i += 1
            story.append(Preformatted("\n".join(block), STYLE["Code"]))
        elif (match := re.match(r"!\[(.*)\]\(([^)]+)\)$", line)):
            story.append(figure(match.group(2), match.group(1)))
        elif line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = [v.strip() for v in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", v) for v in row):
                    rows.append(row)
                i += 1
            story.extend([table(rows), Spacer(1, 8)])
            continue
        elif line.startswith("# "):
            story.append(p(line[2:], "Title"))
        elif line.startswith("## "):
            story.append(p(line[3:], "Heading"))
        elif line.startswith("### "):
            story.append(p(line[4:], "Subheading"))
        elif line.startswith("> "):
            story.append(p(line[2:], "Kicker"))
        else:
            paragraph = [line]
            while (i + 1 < len(lines) and lines[i + 1].strip()
                   and not re.match(r"^(#|\||!\[|```|<!--|> )", lines[i + 1].strip())):
                i += 1
                paragraph.append(lines[i].strip())
            story.append(p(" ".join(paragraph)))
        i += 1
    return story


def page(canvas, document):
    canvas.saveState()
    canvas.setStrokeColor(LINE)
    canvas.setLineWidth(.6)
    canvas.line(MARGIN, PAGE_HEIGHT - 36, PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 36)
    canvas.setFont("Helvetica-Bold", 8)
    canvas.setFillColor(ACCENT)
    canvas.drawString(MARGIN, PAGE_HEIGHT - 25, "NCU / EE8091")
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawRightString(PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 25,
                          "HW08 / TRANSMITTER FFE")
    canvas.line(MARGIN, 33, PAGE_WIDTH - MARGIN, 33)
    canvas.drawString(MARGIN, 21, "Anonymous numerical report")
    canvas.drawRightString(PAGE_WIDTH - MARGIN, 21, str(document.page))
    canvas.restoreState()


def build(source, output):
    rl_config.invariant = True
    doc = BaseDocTemplate(
        str(output), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=51, bottomMargin=45, title="HW8 - Transmitter FFE Design and Simulation",
        author="", subject="EE8091 anonymous numerical report",
        creator="Anonymous report builder", pageCompression=1,
    )
    frame = Frame(MARGIN, 45, WIDTH, PAGE_HEIGHT - 96,
                  leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates(PageTemplate(id="Report", frames=frame, onPage=page))
    doc.build(parse_markdown(source))
    print(f"Built {output.name}")


def markdown_table(rows):
    lines = ["| " + " | ".join(map(str, rows[0])) + " |",
             "| " + " | ".join("---" for _ in rows[0]) + " |"]
    lines += ["| " + " | ".join(map(str, row)) + " |" for row in rows[1:]]
    return "\n".join(lines)


def flatten_scalars(value, prefix=""):
    for key, item in value.items():
        label = f"{prefix} / {key}" if prefix else key
        if isinstance(item, dict):
            yield from flatten_scalars(item, label)
        elif isinstance(item, (str, int, float, bool)) or item is None:
            if isinstance(item, float):
                item = f"{item:.6g}"
            yield label.replace("_", " "), str(item)


def report_markdown(summary):
    cases = summary["cases"]
    if len(cases) != 5:
        raise ValueError("The report requires all five assignment configurations")
    params = summary["parameters"]
    for key, expected in {"tb_s": 25e-12, "tr_s": 5e-12, "dt_s": .5e-12,
                          "source_vpp_v": 2, "kaiser_alpha": .5,
                          "fir_truncation_s": 20e-9}.items():
        if not math.isclose(params[key], expected, rel_tol=1e-10, abs_tol=0):
            raise ValueError(f"Review the report narrative for changed {key}")
    # This report's narrative is for the assignment's physical parameters.
    # Recorded simulation size and measured numerical results are inserted below.
    tap_rows = [["Case (pre / post)", "Coefficients, left to right", "Design residual (V)"]]
    metric_rows = [["Case", "SBR main (mV)", "Empirical gap (mV)", "Tx swing (V)"]]
    for i, case in enumerate(cases, 1):
        tap_rows.append([
            f"{i} ({case['npre']} / {case['npost']})",
            ", ".join(f"{c:+.6f}" for c in case["taps"]),
            f"{case['zf_matrix_residual_v']:.3g}",
        ])
        metric_rows.append([
            str(i), f"{case['main_cursor_v'] * 1e3:.3f}",
            f"{case['empirical_sampled_vertical_gap_v'] * 1e3:.3f}",
            f"{case['tx_swing_v']:.6f}",
        ])
    validation = [["Numerical check / recorded quantity", "Value"]]
    check_labels = {
        "ola_vs_direct_max_abs_v": "Sectional vs direct convolution, max error (V)",
        "continuous_sbr_vs_sampled_pulse_max_abs_v": "Analytic vs sampled pulse, max error (V)",
        "continuous_sbr_vs_sampled_pulse_rms_v": "Analytic vs sampled pulse, RMS error (V)",
        "continuous_sbr_vs_sampled_pulse_relative_peak": "Analytic vs sampled pulse, error / peak",
        "sampled_pulse_fft_vs_folded_linear_filter_max_abs_v": "Sampled DFT vs folded linear convolution, max error (V)",
        "sampled_pulse_periodic_vs_linear_boundary_max_abs_v_before_t2": "Periodic vs unfurled pulse, max error before 20 ns (V)",
        "dt_convergence_main_window_coarse_max_abs_v": "Main-window sampled pulse error at 0.5 ps (V)",
        "dt_convergence_main_window_fine_max_abs_v": "Main-window sampled pulse error at 0.25 ps (V)",
        "dt_convergence_fine_dt_s": "Refined time spacing (s)",
    }
    for key, label in check_labels.items():
        if key in summary["verification"]:
            validation.append([label, f"{summary['verification'][key]:.6g}"])
    # Avoid a layout avalanche from accidentally dumping a new large schema.
    if len(validation) > 23:
        raise ValueError("Review new verification fields before rebuilding report")
    evidence_path = RESULTS / "test_results.json"
    if not evidence_path.exists():
        raise FileNotFoundError("Run and record the independent test suite before building the report")
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if not evidence.get("successful"):
        raise ValueError("Refusing a report that claims unsuccessful tests passed")
    bits = params["payload_bits"]
    seed = params.get("seed", 0)
    best = max(cases, key=lambda c: c["empirical_sampled_vertical_gap_v"])
    best_index = cases.index(best) + 1
    channel = summary["channel"]
    case4 = cases[3]
    snapshots = case4["snapshots"]
    first_tx = next(s for s in snapshots if s["signal"] == "tx" and s["snapshot"] == 1)
    first_rx = next(s for s in snapshots if s["signal"] == "rx" and s["snapshot"] == 1)
    snapshot_rx_offset_ns = (first_rx["start_time_s"] - first_tx["start_time_s"]) * 1e9
    environment = ", ".join(f"{key} {value}" for key, value in summary.get("environment", {}).items())
    return f"""# Transmitter FFE design and simulation

> HW8 - Introduction to Signal Integrity and Electromagnetic Compatibility

This report designs and evaluates a transmitter feedforward equalizer (FFE) for
the supplied 50 cm stripline channel. All five required configurations are
recomputed from the same random bit sequence. The post-FFE transmitter waveform
is constructed explicitly and propagated through the channel by sectional
convolution. The equivalent FFE frequency response is shown only for reference.

## Model and parameters

| Quantity | Setting |
| --- | --- |
| Channel | Supplied HW7_stripline_50cm.s2p, 50 ohm reference |
| Bit rate / modulation | 40 Gb/s, binary NRZ |
| Unit interval | 25 ps |
| Source bit levels | 0 V and 2 V before FFE |
| Source-to-load transfer | H(f) = S21(f) / 2, applied once |
| Rise and fall duration | 5 ps, 0 to 100%, for every transition size |
| Time sample spacing | 0.5 ps, 50 samples/UI |
| FFE tap spacing | 1 UI |
| Kaiser parameter | alpha = 0.5; beta = pi times alpha |
| Retained causal channel FIR | 20 ns |
| Random sequence | {bits:,} measured bits; seed {seed}; guard bits excluded |

### Calculation sequence

Read the Touchstone frequency grid and complex S21. Apply the matching factor and
Kaiser window, construct a Hermitian spectrum, and obtain a real sampled channel
impulse response. Compute the single-bit response (SBR), sample its cursors one
UI apart, solve the finite zero-forcing system, and normalize the taps.

Convolve the binary voltage levels with the taps. At each bit boundary, interpolate
from the previous level to the new level over exactly 5 ps. Convolve this actual
time waveform with the channel FIR in zero-padded sections, adding overlapping
tails. Reuse the same bit pattern for every configuration.

<!-- pagebreak -->

## Tap design and numerical results

Index channel cursors relative to the unequalized SBR maximum. For taps indexed
from -Npre through Npost, use the following finite design system:

```text
A[i,j] = a[r_i - n_j]         (unselected cursors set to zero)
A c_raw = e_main
c = c_raw / sum(abs(c_raw))
```

The main tap is at vector position Npre. A causal implementation stores the
leftmost tap first and adds Npre UI of latency. A precursor does not require
access to an unavailable future bit.

{markdown_table(tap_rows)}

The design residual measures A c minus the normalized target at the **retained
design cursors**. It is not the residual of the full physical channel. Omitted
channel cursors can contribute at these same times, and normalization reduces
the main-cursor voltage from the unnormalized unit target.

For example, case 4 retains a maximum selected off-main contribution of
{case4['physical_selected_residual_v'] * 1e3:.3f} mV at the original design phase.
The absolute sum of all other UI-spaced analytic cursors at its chosen sampling
phase is {case4['full_residual_isi_l1_v'] * 1e3:.3f} mV. These physical-response
quantities remain much larger than the finite matrix residual.

{markdown_table(metric_rows)}

The empirical gap is the minimum sampled high voltage minus the maximum sampled
low voltage in this finite, noiseless random run, at the phase selected by each
equalized SBR. A negative value means the two observed distributions overlap.
It is neither a BER estimate nor a 3-sigma eye-height measurement or a guarantee
for every possible bit sequence. No gain or vertical offset is fitted between
cases. Source swing is bounded by the absolute-sum tap constraint.

For this recorded run, case {best_index} has the largest empirical gap,
{best['empirical_sampled_vertical_gap_v'] * 1e3:.3f} mV. Cases 1 and 2 have
overlapping sampled distributions. The finite ZF system targets cancellation,
not direct maximization of this gap; increasing the tap count does not make the
gap improve monotonically.

<!-- pagebreak -->

## Channel response in the time domain

![Unequalized single-bit response and UI-spaced cursor samples. The pulse phase uses Tb + Tr, giving a causal trapezoid.](../results/single_bit_response.png)

![Full equalized SBRs. Finite tap designs reduce selected ISI, but residual precursor and postcursor contributions remain.](../results/equalized_sbr.png)

<!-- pagebreak -->

## Eye diagrams: cases 1 and 2

![Case 1: unequalized reference. All eye plots use voltage increasing upward and two repeated UI of the phase-folded probability histogram.](../results/eye_case1_no_ffe.png)

![Case 2: one precursor tap and no postcursor taps. The main cursor is aligned consistently using this case's SBR.](../results/eye_case2_pre1.png)

<!-- pagebreak -->

## Eye diagrams: cases 3 and 4

![Case 3: one postcursor tap and no precursor taps. Transmitter de-emphasis reduces low-frequency source voltage under a fixed swing budget.](../results/eye_case3_post1.png)

![Case 4: one precursor and one postcursor. Its two transmitter snapshots and two receiver snapshots follow.](../results/eye_case4_pre1_post1.png)

<!-- pagebreak -->

## Eye diagram: case 5 and frequency reference

![Case 5: one precursor and two postcursors. More taps address more selected cursors; they do not guarantee the largest sampled eye opening.](../results/eye_case5_pre1_post2.png)

![Normalized tap coefficients and equivalent FFE frequency responses. The waveform calculation uses explicit transmitter levels followed by channel convolution.](../results/ffe_taps.png)

<!-- pagebreak -->

## Case 4: transmitter waveforms

![Transmitter snapshot 1. All level transitions use the same 5 ps duration, including transitions between intermediate FFE levels.](../results/case4_tx_snapshot1.png)

![Transmitter snapshot 2, from a second part of the same deterministic random sequence. These two views satisfy the required Tx snapshots.](../results/case4_tx_snapshot2.png)

<!-- pagebreak -->

## Case 4: receiver waveforms

![Receiver snapshot 1. Time is relative to the recorded snapshot start, and all voltages retain the simulated source-to-load scale.](../results/case4_rx_snapshot1.png)

![Receiver snapshot 2, from a second interval. These two views satisfy the required Rx snapshots.](../results/case4_rx_snapshot2.png)

Tx and Rx axes each start at their own local zero. Each paired Rx slice starts
{snapshot_rx_offset_ns:.6f} ns after its Tx slice, using the SBR-derived sampling
latency minus half a UI. Their absolute start samples, start times and source
payload indices are recorded in summary.json; the four CSV files preserve the
plotted samples. The local zeros do not denote simultaneous physical times.

<!-- pagebreak -->

## Verification and interpretation

{markdown_table(validation)}

The independent regression run recorded **{evidence['tests_run']} passing tests**
({evidence['hw08_tests']} for HW8 and {evidence['tests_run'] - evidence['hw08_tests']} existing regression tests),
with {evidence.get('failures', 0)} failures and {evidence.get('errors', 0)} errors.
It covers transfer-function scaling, waveform construction, tap indexing,
section boundaries, eye normalization and existing repository analyses. The
recorded evidence is in results/test_results.json.

The analytic SBR and the sampled finite-ramp convolution use different
continuous-time versus discrete-time approximations. Their measured discrepancy
is retained in the verification results. Zero filling above the measured
bandwidth, Kaiser tapering, and finite causal impulse truncation also limit
accuracy. A histogram with no added noise or jitter cannot establish link BER.

The finite 20 ns FIR has DC gain {channel['finite_fir_dc_gain']:.9f}, compared
with the source-to-load DC value {channel['dc_gain']:.9f}; the relative DC deficit
is {channel['dc_gain_relative_truncation_error'] * 100:.4f}%. It is recorded rather
than removed by an amplitude renormalization. The refined pulse check reduces
the sampling discrepancy; it does not independently verify the random-data eye.

The optional ADS ChannelSim exercise is not run here. No numerical measured or ADS
waveform reference is supplied for this HW8 data set, so the checks establish numerical consistency
of the stated model, rather than an independent hardware or circuit-simulator
match. Tx FFE emphasizes transitions by sacrificing steady-state amplitude under
the swing constraint; an open eye need not have a larger absolute signal level.

### Reproduce

```sh
python -m pip install -r requirements.txt
python homework/hw08/verify_inputs.py
python homework/hw08/analyze.py
python -m unittest discover -s tests -v
python -m pip install -r homework/hw08/report/requirements.txt
python homework/hw08/report/build_report.py
```

The supplied assignment and incomplete template are retained under instructions/.
The channel and legacy Touchstone helper are supplied references; analyze.py is
the rewritten executable solution. The original student report and source are
not published. The report, PDF author field, figures and generated numerical data
contain no student name, student number, email address or local machine path.
Source-byte hashes are recorded in instructions/manifest.json.

Recorded runtime packages: {environment}.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="refresh report.md from results/summary.json")
    args = parser.parse_args()
    source = HERE / "report.md"
    if args.refresh or not source.exists():
        summary = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
        source.write_text(report_markdown(summary), encoding="utf-8")
    build(source, HERE / "report.pdf")


if __name__ == "__main__":
    main()
