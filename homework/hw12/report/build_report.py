"""Build the anonymous HW12 report from recorded, verified simulation outputs.

Run this after the full default simulation and independent verification suite.
No names, student numbers, machine paths, or user-specific metadata are added.
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
MARGIN = 43
WIDTH = PAGE_WIDTH - 2 * MARGIN
INK = colors.HexColor("#193044")
ACCENT = colors.HexColor("#087B8A")
MUTED = colors.HexColor("#506473")
LINE = colors.HexColor("#D4E1E6")
PALE = colors.HexColor("#F0F6F7")
CASE_KEYS = ("case1_no_tx_jitter", "case2_sj", "case3_rj", "case4_sj_rj")


def styles():
    definitions = {
        "Body": dict(fontName="Helvetica", fontSize=9.6, leading=13.2,
                     textColor=INK, spaceAfter=8),
        "Title": dict(fontName="Helvetica-Bold", fontSize=25, leading=29,
                      textColor=INK, spaceAfter=13),
        "Heading": dict(fontName="Helvetica-Bold", fontSize=16, leading=21,
                        textColor=INK, spaceAfter=11, keepWithNext=True),
        "Subheading": dict(fontName="Helvetica-Bold", fontSize=10.7, leading=14,
                           textColor=ACCENT, spaceBefore=6, spaceAfter=7,
                           keepWithNext=True),
        "Caption": dict(fontName="Helvetica", fontSize=8.4, leading=11.2,
                        textColor=MUTED, spaceBefore=4, spaceAfter=12),
        "TableBody": dict(fontName="Helvetica", fontSize=8.1, leading=10.7,
                          textColor=INK),
        "TableHead": dict(fontName="Helvetica-Bold", fontSize=8.1, leading=10.7,
                          textColor=colors.white),
        "Kicker": dict(fontName="Helvetica-Bold", fontSize=9, leading=12,
                       textColor=ACCENT, spaceAfter=10),
        "Code": dict(fontName="Courier", fontSize=7.8, leading=11,
                     textColor=INK, backColor=PALE, borderPadding=8,
                     spaceBefore=5, spaceAfter=12),
    }
    return {key: ParagraphStyle(key, **value) for key, value in definitions.items()}


STYLE = styles()


def rich(text):
    text = html.escape(str(text), quote=False)
    text = re.sub(r"`([^`]+)`", r'<font name="Courier" size="8.6">\1</font>', text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    return re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", text)


def paragraph(text, style="Body"):
    return Paragraph(rich(text), STYLE[style])


def table(rows):
    count = len(rows[0])
    if count == 2:
        widths = [WIDTH * .40, WIDTH * .60]
    elif count == 3:
        widths = [WIDTH * .42, WIDTH * .29, WIDTH * .29]
    elif count == 5:
        widths = [WIDTH * .19] + [WIDTH * .2025] * 4
    else:
        widths = [WIDTH / count] * count
    cells = [[paragraph(cell, "TableHead" if i == 0 else "TableBody")
              for cell in row] for i, row in enumerate(rows)]
    result = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LINEBELOW", (0, -1), (-1, -1), .6, LINE),
    ]))
    return result


def figure(relative_path, caption):
    path = (HERE / relative_path).resolve()
    width, height = ImageReader(str(path)).getSize()
    limit = 420 if path.name == "rx_eyes.png" else 245
    scale = min(WIDTH / width, limit / height)
    image = Image(str(path), width=width * scale, height=height * scale)
    image.hAlign = "CENTER"
    return KeepTogether([image, paragraph(caption, "Caption")])


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
                row = [value.strip() for value in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", value) for value in row):
                    rows.append(row)
                i += 1
            story.extend([table(rows), Spacer(1, 8)])
            continue
        elif line.startswith("# "):
            story.append(paragraph(line[2:], "Title"))
        elif line.startswith("## "):
            story.append(paragraph(line[3:], "Heading"))
        elif line.startswith("### "):
            story.append(paragraph(line[4:], "Subheading"))
        elif line.startswith("> "):
            story.append(paragraph(line[2:], "Kicker"))
        else:
            block = [line]
            while (i + 1 < len(lines) and lines[i + 1].strip()
                   and not re.match(r"^(#|\||!\[|```|<!--|> )", lines[i + 1].strip())):
                i += 1
                block.append(lines[i].strip())
            story.append(paragraph(" ".join(block)))
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
    canvas.drawRightString(PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 25, "HW12 / JITTER ANALYSIS")
    canvas.line(MARGIN, 33, PAGE_WIDTH - MARGIN, 33)
    canvas.drawString(MARGIN, 21, "Anonymous numerical report")
    canvas.drawRightString(PAGE_WIDTH - MARGIN, 21, str(document.page))
    canvas.restoreState()


def build(source, output):
    rl_config.invariant = True
    doc = BaseDocTemplate(
        str(output), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=51, bottomMargin=45,
        title="HW12 - Jitter Simulation and Basic Analysis", author="",
        subject="Anonymous signal-integrity numerical report",
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
    lines.extend("| " + " | ".join(map(str, row)) + " |" for row in rows[1:])
    return "\n".join(lines)


def fnum(value):
    return f"{value:.6g}"


def load_evidence():
    path = RESULTS / "test_results.json"
    if not path.exists():
        raise FileNotFoundError("Run the independent tests and record test_results.json first")
    evidence = json.loads(path.read_text(encoding="utf-8"))
    if evidence.get("successful") is not True:
        raise ValueError("Independent verification must pass before this report is built")
    return evidence


def report_markdown(summary):
    params = summary["parameters"]
    cases = [summary["cases"][key] for key in CASE_KEYS]
    evidence = load_evidence()
    verified = json.loads((RESULTS / "verification.json").read_text(encoding="utf-8"))
    if not verified.get("successful") or not verified.get("assignment_sufficient"):
        raise ValueError("Full saved-artifact verification must pass before publication")
    if verified.get("evaluated_bits") != params["bits"] or verified.get("points_per_ui") != params["points_per_ui"]:
        raise ValueError("Verification record does not match the simulation size or grid")
    for key, expected in {
        "ui_seconds": 25e-12, "rise_seconds": 5e-12, "sj_peak_seconds": 3e-12,
        "sj_frequency_hz": 1e6, "rj_sigma_seconds": 1e-12,
        "kaiser_a": .5, "truncation_seconds": 40e-9,
    }.items():
        if not math.isclose(params[key], expected, rel_tol=1e-10, abs_tol=0):
            raise ValueError(f"Review report narrative for changed {key}")
    if params["bits"] < 2**20:
        raise ValueError("The published report requires at least 2**20 evaluated bits")

    variances = [[case[side]["variance_seconds2"] * 1e24 for side in ("tx", "rx")]
                 for case in cases]
    rms = [[case[side]["rms_seconds"] * 1e12 for side in ("tx", "rx")]
           for case in cases]
    metric_rows = [["Case", "Tx variance (ps^2)", "Rx variance (ps^2)",
                    "Tx RMS (ps)", "Rx RMS (ps)"]]
    for i, (v, r) in enumerate(zip(variances, rms), 1):
        metric_rows.append([str(i), fnum(v[0]), fnum(v[1]), fnum(r[0]), fnum(r[1])])
    theory_rows = [["Tx case", "Theoretical variance (ps^2)", "Theoretical RMS (ps)"]]
    for label, variance in [("1: no jitter", 0), ("2: SJ", 4.5), ("3: RJ", 1), ("4: SJ + RJ", 5.5)]:
        theory_rows.append([label, fnum(variance), fnum(math.sqrt(variance))])

    q1_lhs, q1_rhs = variances[1][0] + variances[2][0], variances[3][0]
    q2_lhs, q2_rhs = variances[0][1] + variances[1][0], variances[1][1]
    q3_lhs, q3_rhs = variances[0][1] + variances[2][0], variances[2][1]
    q4_lhs, q4_rhs = variances[0][1] + variances[3][0], variances[3][1]
    equations = [["Comparison", "LHS (ps^2)", "RHS (ps^2)", "RHS - LHS (ps^2)", "Residual / LHS (%)"]]
    for label, lhs, rhs in [("Q1: Tx SJ + RJ", q1_lhs, q1_rhs),
                            ("Q2: DDJ + Tx SJ", q2_lhs, q2_rhs),
                            ("Q3: DDJ + Tx RJ", q3_lhs, q3_rhs),
                            ("Q4: intended Rx4", q4_lhs, q4_rhs)]:
        equations.append([label, fnum(lhs), fnum(rhs), f"{rhs-lhs:+.6g}", f"{100*(rhs-lhs)/lhs:+.5f}"])
    checks = summary.get("checks", {})
    check_rows = [["Recorded numerical check", "Value"]]
    for key, value in checks.items():
        if isinstance(value, (str, bool, int, float)):
            label = key.replace("_", " ")
            check_rows.append([label, fnum(value) if isinstance(value, float) else str(value)])
    if len(check_rows) > 13:
        check_rows = check_rows[:13]
    channel_rows = [["Channel diagnostic", "Recorded value"]]
    for key, value in summary.get("channel", {}).items():
        if isinstance(value, (str, bool, int, float)) and not any(x in key.lower() for x in ("path", "file")):
            channel_rows.append([key.replace("_", " "), fnum(value) if isinstance(value, float) else str(value)])
    if len(channel_rows) > 11:
        channel_rows = channel_rows[:11]
    bits = params["bits"]
    dt_ps = params["dt_seconds"] * 1e12
    warmup = params["warmup_bits"]
    tests_run = evidence.get("tests_run", "recorded")
    max_theory_percent = 100 * max(
        abs(value["relative_error"])
        for value in verified["tx_theoretical_checks"].values()
        if isinstance(value, dict) and "relative_error" in value
    )
    block_path = RESULTS / "block_size_check.json"
    block_section = ""
    if block_path.exists():
        block = json.loads(block_path.read_text(encoding="utf-8"))
        if not block.get("passed"):
            raise ValueError("Recorded block-size invariance check did not pass")
        block_section = (
            f"A separate end-to-end check evaluated {block['bits']:,} bits using block sizes "
            f"{block['block_sizes'][0]} and {block['block_sizes'][1]}. All four cases produced "
            "exactly equal integer Tx/Rx eye and crossing-count arrays; only floating-point "
            "moment sums use roundoff tolerances. This checks boundary handling independently "
            "of the full-run block choice. Evidence: `results/block_size_check.json`."
        )
    channel = summary["channel"]
    dc_error_percent = 100 * channel["fir_dc_relative_error"]
    rx_slice_counts = ", ".join(f"{int(case['rx']['observation_count']):,}" for case in cases)
    crossing_plot = RESULTS / "crossing_diagnostics.png"
    crossing_section = ""
    if crossing_plot.exists():
        crossing_section = """![Independent crossing-event diagnostic. Interpolated threshold crossings count edge events, while an eye slice counts time spent in a narrow voltage bin.](../results/crossing_diagnostics.png)

The threshold-slice method follows the assignment and lecture. A narrow voltage
bin weights a transition approximately in proportion to the inverse of its
local slope. Receiver edge slopes depend on the data pattern, so its normalized
eye slice need not exactly match the histogram of individual threshold-crossing
events. The additional event-based diagnostic exposes this estimator difference;
it is not substituted silently for the requested eye-slice result.
"""

    return f"""# Jitter simulation and basic analysis

> HW12 - Introduction to Signal Integrity and Electromagnetic Compatibility

This report recomputes all four required 40 Gb/s NRZ experiments using the supplied
two-port channel. The transmitter contains no jitter, sinusoidal jitter (SJ),
random jitter (RJ), or both. Receiver eye diagrams, one-UI transmitter and receiver
jitter distributions, and their variances are obtained from the same channel and
common random data. Every numerical value below comes from the saved run.

## Physical and numerical model

| Quantity | Setting |
| --- | --- |
| Bit rate / unit interval | 40 Gb/s NRZ / 25 ps |
| Source voltage levels | 0 and 1 V |
| Source resistance / capacitance | 20 ohm / 0.3 pF |
| Load resistance / capacitance | 85 ohm / 0.3 pF |
| Rise and fall duration | 5 ps |
| SJ | 3 ps peak amplitude, 1 MHz |
| RJ | Zero-mean Gaussian, 1 ps standard deviation |
| Evaluated data | {bits:,} bits per case; seed {params['seed']} |
| Excluded startup | {warmup:,} warmup bits per case |
| Time samples | {params['points_per_ui']} per UI; {dt_ps:g} ps spacing |
| Eye voltage bin | {params['dv_volts']*1000:g} mV |
| Kaiser window / retained impulse | alpha = 0.5 / 40 ns |

The SJ period is 1 microsecond, or 40,000 UI. The measured record spans
{bits / 40000:.4f} SJ cycles. The source resistance, load resistance and shunt
capacitances are included in the source-to-load transfer H(f); H is not simply
the channel's S21. There is no receiver additive noise or clock recovery model.

<!-- pagebreak -->

## Channel and waveform calculation

The Touchstone data use a 50 ohm reference impedance. At each frequency, modified
nodal analysis solves the external circuit together with the two-port relation
`(I-S)V = Z0(I+S)Iport`. The channel is oriented with S21 transferring excitation
from port 1 to port 2. The Kaiser-weighted positive-frequency data are completed
by conjugate symmetry to obtain a real impulse response on the simulation grid.
The causal 40 ns interval is retained as a finite impulse response (FIR).

![Loaded channel transfer and impulse response used in the four experiments. The original S21 is a reference; the source-to-load H includes both terminations and their shunt capacitances.](../results/channel_diagnostics.png)

{markdown_table(channel_rows)}

The retained FIR has DC gain {channel['fir_dc_gain']:.9f}, compared with
{channel['dc_gain']:.9f} for the frequency-domain circuit. The signed relative
difference is {dc_error_percent:.4f}%. The FIR is not rescaled to hide this
windowing/causal-truncation difference. The receiver threshold remains the
assignment's frequency-domain DC definition.

<!-- pagebreak -->

## Four receiver statistical eyes

![Receiver statistical eyes for cases 1-4. Every input includes the same finite rise and fall duration. Colors represent the per-time-column voltage probability; zero-probability bins are masked only for logarithmic display.](../results/rx_eyes.png)

Case 1 already has receiver timing spread because channel memory produces
data-dependent jitter (DDJ). The ideal transmitter does not remove channel ISI.
Case 2 adds the slowly varying SJ displacement, case 3 adds independent random
transition displacement, and case 4 combines both source impairments.

The simulation processes bounded blocks and carries the FIR state between them.
Guard samples are excluded before accumulating exactly the requested number of
payload bits. Every time column of each eye is normalized using its actual count.
No artificial probability is added to the analysis arrays; zero-count bins remain
zero. The two displayed UI repeat the same one-UI statistical eye.

Source edge displacements are rounded to the nearest time sample before waveform
construction. The maximum timing-rounding error is {dt_ps/2:g} ps. This matches
the assignment's sampled-edge model and avoids the original truncation-toward-zero
bias. It also prevents a narrow transmitter voltage slice from sparsely sampling
continuous sub-grid edge positions. The rounding is a numerical approximation,
not a physical bound on Gaussian RJ.

<!-- pagebreak -->

## One-UI jitter distributions

![Normalized transmitter and receiver jitter PMFs for all four cases. Each bin contains probability, and the probabilities sum to one.](../results/tx_rx_pmf.png)

![The same distributions expressed as PDFs. Each PMF bin is divided by the time-bin width; density times the bin width sums to one.](../results/tx_rx_pdf.png)

The transmitter threshold is 0.5 V. The receiver threshold is half the loaded
DC output, `real(H(0))/2`. A single one-UI interval is selected by circularly
centering the baseline transition region; the reference is reused across cases
at each observation side. This prevents losing bins when the transition straddles
the phase origin. A PDF value has inverse-time units; a PMF bin is dimensionless.

The four receiver slices contain {rx_slice_counts} occupied voltage-bin samples,
respectively. These counts, rather than the total transmitted bit count alone,
set the smallest nonzero increments in the normalized slice PMFs.

<!-- pagebreak -->

## Numerical statistics and transmitter theory

These are the assignment's **normalized eye-slice statistics**. The separate
crossing-event diagnostic uses a different estimator and is identified explicitly.

{markdown_table(metric_rows)}

For samples at bin times t_i and a normalized PMF p_i, the statistics are

```text
mean = sum(t_i * p_i)
variance = sum((t_i - mean)^2 * p_i)
centered RMS = sqrt(variance)
PDF_i = p_i / dt
```

The reported RMS is the standard deviation after removing the mean timing
offset. Raw RMS would also include the square of a fixed delay or phase offset.
Every time bin is spaced by the same dt used by the waveform simulator.

{markdown_table(theory_rows)}

For uniformly sampled sinusoidal phase, SJ has zero mean and variance A^2/2.
Its density is `1/(pi*sqrt(A^2-j^2))` for `abs(j)<A`, with more probability near
the turning points. Independent zero-mean Gaussian RJ has variance sigma_RJ^2.
Thus the combined transmitter variance is `3^2/2 + 1^2 = 5.5 ps^2`, and its
centered RMS is 2.345208 ps. Finite phase coverage, random sampling and time/voltage
binning account for small measured differences from these limiting values.

<!-- pagebreak -->

## Answers to the four assignment questions

{markdown_table(equations)}

### 1. Does Tx case 2 variance plus Tx case 3 variance equal Tx case 4 variance?

**Yes in the independent-source population model.** For any two random timing
errors, `Var(X+Y)=Var(X)+Var(Y)+2*Cov(X,Y)`. Here SJ phase and the Gaussian random
source are independent, so their covariance is zero in the ensemble. The recorded
comparison is {q1_lhs:.6g} versus {q1_rhs:.6g} ps^2, with a residual of
{q1_rhs-q1_lhs:+.6g} ps^2. Finite-record and histogram differences need not be zero.
Nearest-grid rounding also matters: rounding `SJ+RJ` does not always equal the
sum of separately rounded SJ and RJ. Thus the sampled four-case experiment is
not an exact algebraic sum even before finite-sample covariance is considered.

### 2. Does Rx case 1 variance plus Tx case 2 variance equal Rx case 2 variance?

**Approximately, when SJ acts as a common slow time shift through the channel.**
A constant displacement of a waveform passes through an LTI channel as the same
displacement. A 1 MHz SJ changes slowly relative to the 25 ps UI and the dominant
channel memory, so its effect is close to that common-delay model. If this shift
is independent of the data-dependent receiver timing, their variances add.
The recorded comparison is {q2_lhs:.6g} versus {q2_rhs:.6g} ps^2, residual
{q2_rhs-q2_lhs:+.6g} ps^2. The full channel plus threshold operation and histogram
estimator do not make exact equality a universal identity.

### 3. Does Rx case 1 variance plus Tx case 3 variance equal Rx case 3 variance?

**There is no general equality.** Independent edge-to-edge RJ changes the widths
and relative positions of neighboring pulses. Their overlapping channel responses
alter both the crossing voltage slope and the crossing time. Receiver timing
therefore cannot generally be represented as the baseline DDJ plus an unchanged
copy of the transmitter RJ. The recorded comparison is {q3_lhs:.6g} versus
{q3_rhs:.6g} ps^2, residual {q3_rhs-q3_lhs:+.6g} ps^2. This is not a failure of the
variance-of-a-sum theorem; the proposed decomposition is the issue.

<!-- pagebreak -->

## Combined jitter and measurement interpretation

### 4. Printed equation and the intended receiver comparison

The supplied question prints `Var(Rx1)+Var(Tx4)=Var(Tx4)`. Taken literally, the
right-hand side cancels the second term and the equation requires zero baseline
DDJ. The measured `Var(Rx1)={variances[0][1]:.6g} ps^2` is nonzero, so the literal
equation is false. This report preserves that reading explicitly.

The likely intended right-hand side is `Var(Rx4)`, continuing questions 2 and 3.
For that interpretation, the prediction is {q4_lhs:.6g} ps^2 and the measured
receiver case 4 value is {q4_rhs:.6g} ps^2, residual {q4_rhs-q4_lhs:+.6g} ps^2.
There is no general equality for the same reason as question 3: the fast RJ is
transformed together with the data pattern by the channel and the threshold
crossing operation. Approximate addition of slow SJ does not make the entire
combined transmitter distribution pass through unchanged.

{crossing_section}

Gaussian RJ is theoretically unbounded. The observed maximum and minimum in a
finite run are not intrinsic peak-to-peak RJ. With about one million bits, the
available low-probability evidence is around the finite-count scale; no BER
extrapolation to 1e-12 or smaller is claimed. The optional ADS ChannelSim exercise
has not been performed, and the assignment's example results are not presented
as fresh ADS verification of this run.

<!-- pagebreak -->

## Verification, corrections, and reproduction

The recorded independent suite completed **{tests_run} tests**, with
{evidence.get('failures',0)} failures and {evidence.get('errors',0)} errors.
It checks the circuit against an independent formulation, waveform boundary
handling, block convolution against a full calculation, histogram normalization,
circular phase selection and statistical calculations. Detailed results are kept
in `results/test_results.json` and the verification source under `tests/`.

{markdown_table(check_rows)}

The saved-artifact verifier independently checked {verified['distributions_checked']} distributions,
including the time grid, all PMF/PDF normalizations and moments, exact evaluated
bit counts and regenerated requested/applied jitter. The largest relative
transmitter-variance difference from ideal theory is {max_theory_percent:.4f}%
in this recorded run.

{block_section}

The rewritten solution removes hard-coded machine paths and an unavailable
legacy helper, avoids retaining all full-length waveforms, excludes startup and
padded samples, and uses one consistent time spacing. It preserves zero analysis
probabilities, computes real-valued thresholds, and treats PMF and PDF units
separately. Every eye uses its own defined voltage grid. The source, figures,
tables, and this report contain no student identity.

The calculation is a finite-bandwidth, windowed and truncated channel model.
The FIR and numerical bin choices are recorded so their effects can be assessed.
Independent numerical agreement validates the stated implementation; it is not
an independent hardware measurement or a new circuit-simulator result.

### Reproduce from the repository root

```sh
python -m pip install -r homework/hw12/requirements.txt
python homework/hw12/src/hw12.py
python -m unittest discover -s homework/hw12/tests -v
python homework/hw12/tests/verify_results.py --output homework/hw12/results/verification.json
python -m pip install -r homework/hw12/report/requirements.txt
python homework/hw12/report/build_report.py --refresh
```

### Sources and provenance

[1] Supplied HW12 assignment, `instructions/HW12_jitter_simulation.pdf`, including
the circuit, four experiments, questions and example answer. [2] Supplied
`instructions/HW12_channel.s2p`. [3] Course handout 12, Jitter Basics, pages 2-12,
for timing error, SJ/RJ distributions and independent-source variance addition.
[4] Course handout 13, Jitter Basics (2), page 2, for the statistical-eye slice
method. The instruction template is retained as a reference, not the executable
solution. The independently regenerated report replaces personally identifying
student reports; it does not claim that their original results were reverified.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="regenerate report.md from verified results")
    args = parser.parse_args()
    source = HERE / "report.md"
    if args.refresh or not source.exists():
        summary = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
        source.write_text(report_markdown(summary), encoding="utf-8")
    build(source, HERE / "report.pdf")


if __name__ == "__main__":
    main()
