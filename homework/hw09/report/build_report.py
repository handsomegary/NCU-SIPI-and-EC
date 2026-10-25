"""Build an anonymous HW9 report from reproducible numerical results.

Run the analysis and tests before this script. The report contains measured
values from results/summary.json; it does not execute or alter the solver.
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
    if n == 2:
        widths = [WIDTH * 0.42, WIDTH * 0.58]
    elif n == 3:
        widths = [WIDTH * 0.18, WIDTH * 0.55, WIDTH * 0.27]
    elif n == 5:
        widths = [WIDTH * v for v in (.10, .20, .20, .25, .25)]
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
    scale = min(WIDTH / width, 255 / height)
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
                          "HW09 / DECISION FEEDBACK EQUALIZER")
    canvas.line(MARGIN, 33, PAGE_WIDTH - MARGIN, 33)
    canvas.drawString(MARGIN, 21, "Anonymous numerical report")
    canvas.drawRightString(PAGE_WIDTH - MARGIN, 21, str(document.page))
    canvas.restoreState()


def build(source, output):
    rl_config.invariant = True
    doc = BaseDocTemplate(
        str(output), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=51, bottomMargin=45,
        title="HW9 - Decision Feedback Equalizer Simulation and Design",
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


def report_markdown(summary):
    """Render the recorded run; reject incompatible model settings."""
    model = summary["model"]
    if (model["ui_ps"], model["rise_time_ps"], model["sample_step_ps"],
            model["source_vpp_v"], model["kaiser_alpha"],
            model["impulse_truncation_ns"], model["channel_transfer"]) != (
            25, 5, .5, 2, .5, 20, "S21/2"):
        raise ValueError("Review narrative for changed physical model")
    params = {"dt_s": model["sample_step_ps"] * 1e-12,
              "payload_bits": model["payload_bits"], "seed": model["seed"],
              "guard_bits_each_end": model["guard_bits_each_side"]}
    cases = []
    for native in summary["cases"]:
        d, m, snap = native["design"], native["metrics"], native["snapshot"]
        if snap["samples"] != 2001 or not math.isclose(snap["duration_ns"], 1.0):
            raise ValueError("Each required snapshot must span exactly 1 ns")
        cases.append({"ntaps": d["tap_count"], "taps_v": d["taps_v"],
                      "threshold_v": d["threshold_v"],
                      "sample_phase_samples": d["phase_samples"],
                      "sample_index": d["latency_samples"],
                      "predicted_gap_v": d["worst_case_eye_v"],
                      "empirical_gap_v": m["observed_eye_v"],
                      "error_count": m["errors"],
                      "error_fraction": m["observed_error_fraction"],
                      "eye_file": native["eye_file"],
                      "waveform_file": snap["file"],
                      "snapshot_start_s": snap["start_time_s"]})
    if [c["ntaps"] for c in cases] != [1, 2, 3]:
        raise ValueError("Expected exactly the one-, two- and three-tap cases")
    design = [["Taps", "Signed feedback coefficients (V)", "Threshold (V)"]]
    timing = [["Taps", "Phase within UI (ps)", "Latency (ns)", "Bound gap (mV)",
               "Observed gap (mV)"]]
    decisions = [["Taps", "Payload errors / bits", "Observed error fraction"]]
    for c in cases:
        design.append([str(c["ntaps"]),
                       ", ".join(f"{v:+.6f}" for v in c["taps_v"]),
                       f"{c['threshold_v']:.6f}"])
        phase_ps = c["sample_phase_samples"] * params["dt_s"] * 1e12
        timing.append([str(c["ntaps"]), f"{phase_ps:.3f}",
                       f"{c['sample_index'] * params['dt_s'] * 1e9:.6f}",
                       f"{c['predicted_gap_v'] * 1e3:.3f}",
                       f"{c['empirical_gap_v'] * 1e3:.3f}"])
        decisions.append([str(c["ntaps"]),
                          f"{c['error_count']:,} / {params['payload_bits']:,}",
                          f"{c['error_fraction']:.6g}"])
    figure_pages = []
    for c in cases:
        n = c["ntaps"]
        figure_pages.append(f"""<!-- pagebreak -->

## {n}-tap DFE: eye and waveform

![Post-DFE summer-output eye for {n} feedback tap(s). The horizontal marker is the selected decision threshold; the sampling phase is centered in each UI. Each eye uses its own voltage range and 420 voltage bins; color limits are shared.](../results/{c['eye_file']})

![Exactly 1 ns (40 UI) of the same received waveform before and after the {n}-tap DFE. Decision samples and threshold retain their simulated voltage scale. Time is relative to this window's start.](../results/{c['waveform_file']})

The snapshot begins at simulation time {c['snapshot_start_s'] * 1e9:.6f} ns.
It has 2,001 waveform samples at 0.5 ps spacing, including both
endpoints from 0 to 1,000 ps. Its local time origin does not remove channel
latency. The full recorded run,
rather than this 40-bit view alone, supplies the decision-error count.
""")
    unavailable = [c for c in cases if c["predicted_gap_v"] <= 0]
    bound_note = ""
    if unavailable:
        labels = ", ".join(str(c["ntaps"]) for c in unavailable)
        bound_note = (f"For the {labels}-tap configuration, the worst-case bound "
                      "does not give a positive separating gap. Its selected "
                      "midpoint is therefore a heuristic operating threshold, "
                      "not a safe threshold for every binary pattern. Zero "
                      "observed errors in the finite random record do not "
                      "contradict the conservative bound.\n\n")
    evidence = verification_narrative(summary)
    decision_text = decision_narrative(summary)
    sensitivity = sensitivity_narrative(summary)
    return f"""# Decision feedback equalizer simulation and design

> HW9 - Introduction to Signal Integrity and Electromagnetic Compatibility

This report designs one-, two-, and three-tap receiver decision feedback
equalizers (DFEs) for the supplied 50 cm stripline channel. Automatic design
uses the channel's sampled single-bit response and residual ISI bounds.
A separate recursive time-domain simulation evaluates actual slicer decisions.
Each required configuration includes a post-DFE eye and a 1 ns comparison
against the same post-channel waveform.

## Model and parameters

| Quantity | Setting |
| --- | --- |
| Channel | Supplied HW7_stripline_50cm.s2p, 50 ohm reference |
| Bit rate / modulation | 40 Gb/s, binary NRZ; 25 ps/UI |
| Source voltage levels | 0 and 2 V |
| Matched source-to-load transfer | H(f) = S21(f) / 2, applied once |
| Transition duration | 5 ps, linear, 0 to 100 percent |
| Time spacing | 0.5 ps; 50 samples/UI |
| Kaiser window | alpha = 0.5; beta = pi times alpha |
| Retained causal channel FIR | 20 ns |
| Evaluated random payload | {params['payload_bits']:,} bits; seed {params['seed']} |
| Guard data | {params['guard_bits_each_end']:,} bits at each end, excluded from scoring |
| Equalization | Receiver DFE only; no Tx FFE or Rx CTLE |
| Noise / jitter | None added |

### Waveform construction and symbol association

Construct the actual NRZ voltage waveform using the finite linear ramps, then
apply the channel by zero-padded sectional convolution. The DFE receives this
waveform; it does not receive the transmitted bits as feedback. The source truth
is used only for evaluation. All three configurations reuse the same waveform.

A single source bit starts at sample kP, where P = 50. Let p be that bit's
sampled trapezoid and g = p * h its response through the same finite channel FIR.
For an absolute sampling latency L, define a[j] = g[L + jP]. The corresponding
source bit is b[k] at received sample kP + L. Positive j denotes a past bit
(postcursor), and negative j a future bit (precursor). L modulo P gives only
the clock phase; the complete L is retained for association with source truth.

<!-- pagebreak -->

## Automatic threshold and sampling-point design

### Question 1: select a good threshold without manual tuning

The single-ended slicer outputs 0 or 1. The summer adds signed feedback
coefficients in volts; a positive postcursor therefore requires a negative tap.
For N taps, cancel the first N postcursors at the candidate sampling phase:

```text
y[kP + L] = sum_j a[j] b[k-j]
c[j] = -a[j],                         j = 1,...,N
d[k] = y[kP + L] + sum_j c[j] bhat[k-j]
bhat[k] = 1 if d[k] >= Vt else 0
```

Assuming previous decisions are correct, combine each selected postcursor with
its feedback coefficient. Call the remaining off-main cursors r[j]. Include
all nonzero precursor and postcursor samples within the finite response,
including uncancelled tails and either sign of cursor. Bounds for arbitrary
binary neighboring symbols are then:

```text
low_max  = sum_(j != 0) max(r[j], 0)
high_min = a[0] + sum_(j != 0) min(r[j], 0)
Vt       = (low_max + high_min) / 2
         = (a[0] + sum_(j != 0) r[j]) / 2
gap      = high_min - low_max
         = a[0] - sum_(j != 0) abs(r[j])
```

The midpoint maximizes the smaller vertical margin to these two bounds. A
positive gap provides a separating threshold for this finite model with correct
past decisions. A nonpositive gap means this bound cannot establish a safe
threshold. The signed sum matters: using half the main cursor alone would
ignore the residual ISI offset. No tap normalization or additional matched
source factor is applied to the feedback voltages.

### Question 2: select a good sampling point without manual tuning

Search one complete UI of discrete phases near the main SBR peak. At each
candidate L, recompute its cursors, feedback taps, threshold, and bound gap.
Choose the phase maximizing a[0] minus the absolute sum of residual ISI.
This objective optimizes the predicted worst-case opening rather than the SBR
peak alone. The phase search uses the discrete response of the same finite
ramps and FIR as the random-data simulation, preserving amplitude and timing.

The bound uses correct past decisions. Actual decisions are evaluated afterward
by the recursive slicer loop, so the design calculation cannot hide error
propagation by substituting transmitted bits into the feedback path.

<!-- pagebreak -->

## Channel response and automatic phase search

![Matched transfer magnitude and sampled single-bit response for the stated source and channel. Its causal trapezoid uses the same finite-ramp definition as the random waveform.](../results/channel_and_sbr.png)

![Predicted residual-ISI gap across the discrete sampling-phase search. Each candidate has its own taps and midpoint threshold; selected phases maximize this finite-model objective.](../results/phase_search.png)

The continuous-pulse frequency-domain response is a separate sampling-convergence
reference. It is not substituted for the sampled finite-FIR response when
selecting taps or checking random-data decisions.

<!-- pagebreak -->

## Selected designs and numerical results

{markdown_table(design)}

{markdown_table(timing)}

Clock phase is shown relative to a 25 ps UI. Latency is the complete source-bit
to received-sample offset. The post-DFE eye is the analog summer output, not the
slicer's ideal 0/1 waveform. The eye plots use individually selected voltage
ranges and therefore different voltage-bin widths. Shared probability-per-bin
color limits do not imply the same probability density. Compare the recorded
voltage gaps rather than apparent eye size or color alone. Each new decision updates the first feedback
interval starting half a UI later; subsequent taps follow one UI apart.

The bound gap assumes correct feedback and every binary neighboring-symbol
combination within the finite model. The observed gap is the minimum sampled
high voltage minus the maximum sampled low voltage in this finite payload,
classified by transmitted truth. They need not be equal because a finite random
record need not realize the worst pattern. Neither quantity is a measured BER,
a noise-inclusive eye height, or a guarantee for an unmodeled physical link.

{bound_note}{markdown_table(decisions)}

The reported counts exclude both guard regions and compare each decision with
its latency-aligned source bit. The decisions come from compensated summer
samples. No gain or vertical offset is fitted between configurations.

{''.join(figure_pages)}

<!-- pagebreak -->

## Question 3: are decisions always correct?

{decision_text}

<!-- pagebreak -->

## Independent verification and numerical consistency

{evidence}

<!-- pagebreak -->

## Model sensitivity, limits, and provenance

{sensitivity}

### Model limits

The reported run is finite and noiseless. Sampling resolution, measured channel
bandwidth, Kaiser tapering, and causal FIR truncation limit accuracy. The
residual-ISI bound covers all cursors of the retained discrete response, not an
infinite physical channel. A finite run without observed errors does not prove
zero physical-link BER. A DFE also cannot cancel precursor ISI using past
symbols, and an incorrect decision can perturb later summer samples.

No measured waveform or independent circuit-simulator reference is supplied for
this case. These checks establish numerical consistency of the stated model.
The optional ADS ChannelSim exercise is not performed.

### Source provenance

The supplied assignment, instructions/HW9_DFE.pdf, requires the three tap counts,
post-DFE eyes, 1 ns waveform comparisons, and the three questions answered here.
The incomplete sample specifies a 0/1 slicer and fixed-rate decision timing.
The supplied channel is the same 50 cm stripline used in the preceding homework.
Input byte hashes are recorded by the reference manifest.

Course lecture "Decision Feedback Equalizer", PDF pages 4-6, supports the summer
node, 0/1 slicer, and postcursor cancellation convention; page 8 discusses DC
level changes; pages 9-10 discuss incorrect decisions and error propagation.
The lecture's page 7 numerical example uses 20 Gb/s and is not copied as this
40 Gb/s design. The residual-bound threshold and phase search above are derived
for this report. Lecture material is consulted as a reference, not redistributed
as an original report or attributed to an individual student.

### Reproduce from the repository root

```sh
python -m pip install -r requirements.txt
python homework/hw09/verify_inputs.py
python homework/hw09/analyze.py --bits 65536 --seed 20261025
python homework/hw09/verify_results.py
python -m unittest discover -s tests -v
python -m pip install -r homework/hw09/report/requirements.txt
python homework/hw09/report/build_report.py --refresh
```

The report is generated from the recorded numerical results and contains no
student identity or machine-specific location. Its PDF author field is empty.
"""


def verification_narrative(summary):
    v = summary["verification"]
    c = [case["verification"] for case in summary["cases"]]
    rows = [["Recorded numerical check", "Value"]]
    rows += [
        ["193-bit sectional versus single-FFT convolution, max error (V)",
         f"{v['sectional_vs_whole_fft_max_abs_v']:.6g}"],
        ["Waveform versus scalar DFE samples, max error over cases (V)",
         f"{max(x['waveform_vs_scalar_max_abs_v'] for x in c):.6g}"],
        ["Waveform versus scalar decision mismatches, all cases",
         str(sum(x['waveform_vs_scalar_decision_mismatches'] for x in c))],
        ["Source-bit cursor reconstruction, max error over cases (V)",
         f"{max(x['source_bits_cursor_reconstruction_max_abs_v'] for x in c):.6g}"],
        ["Eye probability-column sum error, max over cases",
         f"{max(x['column_normalization_max_error'] for x in c):.6g}"],
        ["Analytic versus discrete SBR, max recorded error (V)",
         f"{v['analytic_vs_discrete_sbr_peak_region_max_abs_v']:.6g}"],
        ["20 ns FIR DC gain", f"{v['finite_fir_dc_gain']:.9f}"],
        ["Frequency-domain source-to-load DC gain", f"{v['frequency_dc_gain']:.9f}"],
        ["Finite-FIR minus frequency-domain DC gain", f"{v['finite_fir_dc_difference']:+.6g}"],
    ]
    if (v["sectional_vs_whole_fft_max_abs_v"] >= 1e-10
            or any(x["waveform_vs_scalar_max_abs_v"] >= 1e-10
                   or x["waveform_vs_scalar_decision_mismatches"] != 0
                   or x["source_bits_cursor_reconstruction_max_abs_v"] >= 1e-10
                   or x["column_normalization_max_error"] >= 1e-12 for x in c)):
        raise ValueError("Recorded numerical consistency checks did not pass")
    evidence = json.loads((RESULTS / "test_results.json").read_text(encoding="utf-8"))
    if not evidence["successful"] or evidence.get("failures", 0) or evidence.get("errors", 0):
        raise ValueError("Recorded test suite did not pass")
    independent = json.loads((RESULTS / "verification.json").read_text(encoding="utf-8"))
    if not independent["passed"] or not independent["pre_payload_history_verified"]:
        raise ValueError("Independent result verification did not pass")
    if independent["source_sha256"] != v["input_data_sha256"]:
        raise ValueError("Independent verification used different channel data")
    if independent["payload_bits_per_case"] != summary["model"]["payload_bits"]:
        raise ValueError("Independent verification does not cover this payload")
    checks = [["Taps", "Decisions checked", "Phases checked", "Sample error (V)",
               "Snapshot error (V)"]]
    for item, case in zip(independent["cases"], summary["cases"]):
        n = case["design"]["tap_count"]
        if (item["tap_count"] != n or item["phase_candidates_verified"] != 50
                or item["decisions_verified"] != summary["model"]["payload_bits"]
                or item["snapshot"]["samples_verified"] != 2001
                or item["errors_verified"] != case["metrics"]["errors"]):
            raise ValueError("Independent evidence does not match this case")
        sample_error = max(item["post_channel_samples_max_abs_v"],
                           item["post_dfe_samples_max_abs_v"])
        snapshot_error = max(item["snapshot"]["post_channel_max_abs_v"],
                             item["snapshot"]["post_dfe_max_abs_v"])
        checks.append([str(n), f"{item['decisions_verified']:,}", "50",
                       f"{sample_error:.3g}", f"{snapshot_error:.3g}"])
    text = ("The independent verifier imports no production solver or channel "
            "module. It rebuilds the channel with NumPy's inverse FFT, convolves "
            "the sampled pulse directly, and reconstructs source-bit responses "
            "and scalar decision history including pre-payload clocks. It checks "
            "every persisted payload decision and all 2,001 samples of each "
            "required waveform CSV. Maximum differences are shown below.\n\n" +
            markdown_table(checks) + "\n\n" +
            f"The verifier also reproduces all {independent['stress']['errors_verified']:,} "
            "mistuning errors, its complete 1 ns snapshot, and all 12 persisted "
            "injected-error rows. Eye-file existence and complete sample counts "
            "are checked; image pixels and every histogram bin are not independently "
            "regenerated. Evidence is retained in results/verification.json.\n\n" +
            f"The recorded repository regression run passed {evidence['tests_run']} "
            f"tests, including {evidence['hw09_tests']} HW9 tests, with "
            f"{evidence['failures']} failures and {evidence['errors']} errors. "
            "Outcomes are retained in results/test_results.json.\n\n" +
            "### Additional solver consistency checks\n\n" +
            markdown_table(rows) + "\n\n" +
            "The scalar DFE uses past-decision dot products rather than held "
            "waveform blocks. Source-bit reconstruction separately checks full "
            "latency and waveform association. The analytic and discrete SBRs "
            "use continuous-time versus sampled-ramp approximations; their "
            "measured difference is recorded rather than forced to zero.")
    return text


def sensitivity_narrative(summary):
    rows = [["Taps", "Threshold change (mV)", "Bound-gap change (mV)",
             "Max tap change (mV)"]]
    for n in (1, 2, 3):
        item = summary["verification"]["truncation_20_to_40_ns"][str(n)]
        rows.append([str(n),
                     f"{item['fixed_latency_threshold_change_v'] * 1e3:+.6f}",
                     f"{item['fixed_latency_eye_bound_change_v'] * 1e3:+.6f}",
                     f"{item['max_tap_change_v'] * 1e3:.6f}"])
    return ("### Sensitivity to extending the retained FIR\n\n" +
            markdown_table(rows) + "\n\n" +
            "This independently checked diagnostic extends the FIR from 20 ns "
            "to 40 ns at each selected latency. It recomputes thresholds, taps "
            "and correct-history bounds; it does not repeat the random-data "
            "simulation or re-optimize phase. Baseline eyes and errors remain "
            "the 20 ns model. Values shown as zero are rounded, not asserted "
            "to vanish exactly.")


def decision_narrative(summary):
    cases = summary["cases"]
    bits = summary["model"]["payload_bits"]
    counts = ", ".join(f"{c['design']['tap_count']} tap(s): {c['metrics']['errors']:,}"
                       for c in cases)
    text = (f"Across {bits:,} scored payload bits per baseline configuration, "
            f"the observed wrong-decision counts are {counts}. These counts "
            "compare recursive slicer decisions with the source truth. They "
            "characterize this finite noiseless run, rather than every possible "
            "channel, noise realization or bit sequence.\n\n")
    stress = summary["stress"]
    error = stress["first_error"]
    text += ("### Separately labeled threshold-mistuning example\n\n"
             "Keep the three-tap coefficients and sampling phase, but raise the "
             "threshold by 0.6 times the main cursor. This deliberately changes "
             "the operating point and is not a baseline design result. ")
    text += (f"The mistuned threshold is {stress['threshold_v']:.6f} V; "
             f"the run produces {stress['metrics']['errors']:,} errors in "
             f"{bits:,} evaluated bits.\n\n")
    if error is None:
        text += "No wrong decision was observed in this diagnostic record.\n\n"
    else:
        rows = [["First observed stress error", "Recorded value"],
                ["Payload bit index (zero based)", str(error["payload_bit"])],
                ["Transmitted bit / slicer decision",
                 f"{error['true_bit']} / {error['decision']}"],
                ["Post-channel sample (V)", f"{error['post_channel_v']:.6f}"],
                ["Post-DFE summer sample (V)", f"{error['summer_v']:.6f}"],
                ["Mistuned threshold (V)", f"{error['threshold_v']:.6f}"],
                ["Source context / first context index",
                 "".join(map(str, error["source_context"])) +
                 f" / {error['source_context_start']}"]]
        text += markdown_table(rows) + "\n\n"
    text += (f"![Deliberately mistuned three-tap DFE: an exactly 1 ns waveform "
             f"snapshot. Wrong decisions are marked separately. This stress case "
             f"is distinct from the automatically designed baseline.]"
             f"(../results/{stress['snapshot']['file']})\n\n")
    injection = summary["injected_error"]
    shifts = ", ".join(f"{v * 1e3:+.3f}" for v in
                       injection["next_three_voltage_changes_v"])
    text += ("### Controlled error-propagation diagnostic\n\n"
             f"A second diagnostic forces only payload decision "
             f"{injection['payload_bit']} from {injection['forced_from']} to "
             f"{injection['forced_to']}, using the baseline threshold. "
             f"The next three summer-sample changes are {shifts} mV. "
             f"The remaining payload has "
             f"{injection['extra_decision_changes_after_injection']:,} additional "
             "decision changes relative to baseline. This is an artificial "
             "single-error injection, not a naturally occurring baseline error. "
             "Its voltage changes show how an incorrect past decision enters "
             "later feedback; a voltage perturbation need not cause another "
             "wrong decision. The sample-by-sample comparison is preserved "
             "in results/injected_error.csv.")
    return text


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
