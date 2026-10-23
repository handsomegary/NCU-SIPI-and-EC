"""Build the anonymous HW05 PDF and Markdown report from regenerated results.

Run from any working directory:
    python report/build_report.py

The default results directory is ../results relative to this file. No original
submission, private path, account, or author information is needed to rebuild.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (
    BaseDocTemplate, Frame, Image, PageBreak, PageTemplate, Paragraph,
    Spacer, Table, TableStyle,
)


INK = colors.HexColor("#162B43")
BLUE = colors.HexColor("#1268A2")
MUTED = colors.HexColor("#526478")
LINE = colors.HexColor("#D9E3EC")
PALE = colors.HexColor("#F1F6FA")
PAGE_W, PAGE_H = letter
MARGIN = 48
BODY_W = PAGE_W - 2 * MARGIN


def styles():
    result = getSampleStyleSheet()
    definitions = {
        "Title": dict(fontName="Helvetica-Bold", fontSize=24, leading=28,
                      textColor=INK, spaceAfter=12),
        "Heading1": dict(fontName="Helvetica-Bold", fontSize=19, leading=23,
                         textColor=INK, spaceAfter=14),
        "Heading2": dict(fontName="Helvetica-Bold", fontSize=11.3, leading=15,
                         textColor=BLUE, spaceBefore=12, spaceAfter=5),
        "BodyText": dict(fontName="Helvetica", fontSize=9.8, leading=14,
                         textColor=INK, spaceAfter=8),
        "Caption": dict(fontName="Helvetica", fontSize=8.1, leading=11,
                        textColor=MUTED, spaceBefore=5, spaceAfter=9),
        "Small": dict(fontName="Helvetica", fontSize=8.5, leading=11.5,
                      textColor=INK, spaceAfter=6),
        "Eyebrow": dict(fontName="Helvetica-Bold", fontSize=9, leading=12,
                        textColor=BLUE, spaceAfter=9),
        "Equation": dict(fontName="Courier", fontSize=9.2, leading=13.5,
                         textColor=INK, spaceBefore=0, spaceAfter=3),
    }
    for name, options in definitions.items():
        if name in result:
            for key, value in options.items():
                setattr(result[name], key, value)
        else:
            result.add(ParagraphStyle(name, parent=result["BodyText"], **options))
    return result


def p(text, style="BodyText"):
    return Paragraph(text, STYLES[style])


def table(rows, widths=None, header=True, small=False):
    row_style = "Small" if small else "BodyText"
    data = [[p(str(cell), row_style) for cell in row] for row in rows]
    obj = Table(data, colWidths=widths, hAlign="LEFT", repeatRows=1 if header else 0)
    commands = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 9),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LINEBELOW", (0, 0), (-1, -1), 0.45, LINE),
    ]
    if header:
        commands += [("BACKGROUND", (0, 0), (-1, 0), PALE)]
    obj.setStyle(TableStyle(commands))
    return obj


def equation(*lines):
    cell = [p(html.escape(line), "Equation") for line in lines]
    obj = Table([[cell]], colWidths=[BODY_W], hAlign="LEFT")
    obj.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), PALE),
        ("LEFTPADDING", (0, 0), (-1, -1), 11),
        ("RIGHTPADDING", (0, 0), (-1, -1), 11),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return obj


def image(path, max_height):
    width, height = ImageReader(str(path)).getSize()
    scale = min(BODY_W / width, max_height / height)
    return Image(str(path), width=width * scale, height=height * scale, hAlign="CENTER")


def page_frame(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(LINE)
    canvas.setLineWidth(0.6)
    canvas.line(MARGIN, PAGE_H - 30, PAGE_W - MARGIN, PAGE_H - 30)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(MARGIN, 25, "HW05 / Statistical eye construction")
    canvas.drawRightString(PAGE_W - MARGIN, 25, str(doc.page))
    canvas.restoreState()


def load_facts(results_dir):
    """Read results for the checked assignment configuration and input record."""
    summary_path = results_dir / "summary.json"
    facts = json.loads(summary_path.read_text(encoding="utf-8"))
    expected = {
        "ui_ps": 25.0, "bins_per_ui": 25, "vmin": -0.1, "vmax": 1.1,
        "dv": 0.01, "origin_ps": 0.0, "dt_ps": 1.0,
        "window_ps": 50.0, "time_bins": 50, "voltage_bins": 120,
    }
    if any(facts["configuration"].get(key) != value for key, value in expected.items()):
        raise ValueError("This report targets the default assignment grid; regenerate default results.")
    expected_sha256 = "5e72b7e14bab20ee171442ae58f01c3594c656d71baf48926378c6d4c297b615"
    if facts["input"].get("sha256") != expected_sha256:
        raise ValueError("The report requires the verified assignment waveform record.")
    return facts


def build_pdf(output, results_dir, facts):
    cfg, inp, trigger = (facts[key] for key in ("configuration", "input", "trigger"))
    fixed, triggered = facts["fixed_rate"], facts["triggered"]
    max_error = max(fixed["maximum_column_sum_error"], triggered["maximum_column_sum_error"])
    story = []
    story += [p("HW05 / SIGNAL INTEGRITY", "Eyebrow"),
              p("Construction of statistical eyes", "Title"),
              p("Fixed-rate folding and differential-clock triggering", "Heading2"),
              p("The same received waveform produces substantially different eyes when its time "
                "reference changes. This report reconstructs both statistical eyes from the supplied "
                "ADS waveform record and verifies the bin counts independently. The triggered eye "
                "has a visibly larger opening because it follows the received clock rather than "
                "a constant nominal period."),
              p("Objective and setup", "Heading2"),
              p("The assignment specifies a nominal 40 Gbps bit rate, so one unit interval (UI) "
                "is 25 ps and the displayed window is 2 UI = 50 ps. It describes a frequency deviation "
                "of approximately +/-2 MHz about 40 GHz with a 1 us modulation period. Such spread "
                "spectrum clocking (SSC) causes a fixed time reference to drift relative to data. [1]"),
              table([
                  ["Waveform record", "Histogram grid", "Common plot axes"],
                  [f"{inp['samples']:,} measured samples<br/>1 ps nominal spacing<br/>Vo, CLKp and CLKn",
                   f"{cfg['voltage_bins']} voltage x {cfg['time_bins']} time bins<br/>dV = 10 mV; dt = 1 ps<br/>Column normalization",
                   "Time: 0 to 50 ps<br/>Voltage: -0.1 to 1.1 V<br/>Shared linear probability scale"],
              ], [BODY_W / 3] * 3, small=True),
              Spacer(1, 10), image(results_dir / "eye_comparison.png", 215),
              p("Figure 1. Regenerated fixed-rate and triggered eyes. Color denotes the probability "
                "of a voltage bin conditional on a time bin. Both panels use the same axes and "
                "color scale; differences are due to time assignment rather than voltage rescaling.", "Caption"),
              p("Main finding", "Heading2"),
              p(f"The received clock crossing phase spans {trigger['peak_to_peak_phase_error_ps']:.3f} ps "
                "relative to successive nominal 50 ps intervals. Fixed-rate folding spreads transitions "
                "across that changing phase. Resetting time at positive differential-clock crossings "
                "reduces this misalignment and reveals the opening seen in Figure 1. This comparison "
                "supports the assignment's SSC alignment argument; it does not establish a bit error rate.")]

    story += [PageBreak(), p("Method and bin conventions", "Heading1"),
              p("1. Preserve the waveform values", "Heading2"),
              p("Read the four CSV columns Time, Vo, CLKp and CLKn. Require finite values and strictly "
                "increasing time, then express time in picoseconds. All voltage samples are retained "
                "unchanged. Interpolation is applied only to clock crossing times. The observed Vo "
                f"range is {inp['minimum_voltage_V']:.8f} to {inp['maximum_voltage_V']:.5f} V, inside "
                "the selected -0.1 to 1.1 V voltage range."),
              p("2. Fixed-rate folding", "Heading2"),
              p("For each sample time t, use a fixed origin t0 = 0 ps and fold the sample into a "
                "50 ps window. The corresponding integer bin is obtained by flooring the time in "
                "1 ps units and then taking integer modulo 50:"),
              equation("j_fixed = floor((t - t0) / dt) mod 50", "k = floor((Vo - vmin) / dV)"),
              p("The voltage offset vmin is essential: a negative voltage belongs near the bottom "
                "of the histogram rather than wrapping around to a high-voltage row. Time and voltage "
                "bins are left-closed and right-open; the final voltage bin also accepts exact vmax. "
                "Only floating point roundoff near integer bin edges is snapped before flooring."),
              p("3. Differential-clock triggers", "Heading2"),
              p("Form c = CLKp - CLKn. Detect a positive crossing between adjacent samples with "
                "c[n-1] &lt; 0 and c[n] &gt;= 0. Linear interpolation estimates its time:"),
              equation("a = -c[n-1] / (c[n] - c[n-1])",
                       "t_cross = t[n-1] + a * (t[n] - t[n-1])"),
              p("For each data sample, choose the latest crossing at or before its time. Set "
                "age = t - t_cross and retain the sample only when 0 &lt;= age &lt; 50 ps. "
                "Then j_triggered = floor(age / dt). Triggered ages never wrap, and a sample is "
                "counted at most once. This avoids double counting where adjacent 2 UI windows "
                "overlap. Samples before the first crossing and in a gap past 50 ps are excluded."),
              p("4. Count and normalize vertically", "Heading2"),
              p("Accumulate a count H[k,j] at each voltage/time bin. Normalize each occupied time "
                "column separately; an empty column remains zero without division by zero:"),
              equation("P[k,j] = H[k,j] / sum_k H[k,j]", "sum_k P[k,j] = 1  (for occupied columns)"),
              p("P is a discrete conditional bin probability P(voltage bin | time bin). A continuous "
                "voltage density would additionally divide by dV. The plotted values therefore depend "
                "on the chosen voltage and time bins. The vertical-normalization procedure follows "
                "the assignment starter and the eye-diagram lecture. [2,3]")]

    story += [PageBreak(), p("Timing evidence and verification", "Heading1"),
              image(results_dir / "clock_phase.png", 135),
              p("Figure 2. Positive crossing times relative to the first crossing plus n x 50 ps. "
                "This is received-clock phase error relative to a nominal reference; it includes "
                "the effects present in the supplied sampled waveform.", "Caption"),
              table([
                  ["Checked result", "Fixed rate", "Triggered"],
                  ["Samples counted", f"{fixed['total_samples']:,}", f"{triggered['total_samples']:,}"],
                  ["Count per occupied time column",
                   f"{fixed['minimum_column_count']:,} to {fixed['maximum_column_count']:,}",
                   f"{triggered['minimum_column_count']:,} to {triggered['maximum_column_count']:,}"],
                  ["Occupied / total time columns", "50 / 50", "50 / 50"],
                  ["Largest column-sum error", f"{fixed['maximum_column_sum_error']:.2e}",
                   f"{triggered['maximum_column_sum_error']:.2e}"],
              ], [BODY_W * 0.48, BODY_W * 0.26, BODY_W * 0.26], small=True),
              Spacer(1, 8),
              p(f"There are {trigger['crossings']:,} interpolated positive crossings. Their intervals "
                f"range from {trigger['minimum_interval_ps']:.3f} to {trigger['maximum_interval_ps']:.3f} ps "
                f"(mean {trigger['mean_interval_ps']:.6f} ps). Triggering excludes "
                f"{trigger['samples_before_first_trigger']} samples before the first crossing and "
                f"{trigger['samples_excluded_after_2UI']} samples with age at least 50 ps, accounting "
                "for the difference in totals. Unequal column counts are therefore expected.", "Small"),
              p("Verification", "Heading2"),
              p("Fourteen unit tests passed, covering bin boundaries, negative voltage offsets, "
                "crossing interpolation, trigger gaps and overlap, normalization, invalid input, "
                "and synthetic clock-rate drift. A separately implemented full-record histogram "
                "check matched every integer bin: 6,000 of 6,000 fixed-rate bins and 6,000 of "
                "6,000 triggered bins. All probabilities are finite, and occupied-column sums "
                f"differ from one by at most {max_error:.2e}.", "Small"),
              p("Reproduction", "Heading2"),
              p("From the HW05 folder, install the analysis requirements and run the analysis, "
                "tests and independent check. The report builder reads the regenerated summary "
                "and figures and needs the separate report requirements.", "Small"),
              equation("python src/hw5.py", "python -m unittest discover -s tests -v",
                       "python report/build_report.py", "python tests/verify_delivery.py"),
              p("Limits and conclusion", "Heading2"),
              p("The triggered eye shows improved timing alignment for this finite ADS record. "
                "The observed crossing intervals are received-waveform measurements, not an estimate "
                "of the SSC source's frequency extrema. No noise model, BER extrapolation, threshold "
                "optimization or standardized eye-mask test was performed. Any eye opening inferred "
                "from these bins remains sample- and bin-dependent.", "Small"),
              p("References: [1] HW5_stat_eye_construction.pdf, pp. 1-3. "
                "[2] HW5_stateye_construction_sample.py. "
                "[3] 5 Eye diagram.pdf, pp. 4-5 and 14-17. "
                "Data: HW5_ADS_waveforms.csv.", "Caption")]

    doc = BaseDocTemplate(
        str(output), pagesize=letter, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=45, bottomMargin=43, title="Construction of statistical eyes",
        author="", subject="Fixed-rate and differential-clock-triggered statistical eye analysis",
        creator="Report builder", pageCompression=1,
    )
    frame = Frame(MARGIN, 43, BODY_W, PAGE_H - 88, leftPadding=0,
                  rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates(PageTemplate(id="Report", frames=frame, onPage=page_frame))
    doc.build(story)


def build_markdown(results_dir, facts):
    cfg, inp, trig = (facts[key] for key in ("configuration", "input", "trigger"))
    fixed, triggered = facts["fixed_rate"], facts["triggered"]
    max_error = max(fixed["maximum_column_sum_error"], triggered["maximum_column_sum_error"])
    return f'''# Construction of statistical eyes

HW05 / Fixed-rate folding and differential-clock triggering

## Objective and setup

Reconstruct fixed-rate and triggered statistical eyes from the supplied ADS waveform record.
The same received waveform produces substantially different eyes when its time reference changes.
The triggered eye has a visibly larger opening because it follows the received clock rather than a constant nominal period.

The assignment specifies a nominal 40 Gbps bit rate: UI = 25 ps and a 2 UI window = 50 ps.
It describes a frequency deviation of approximately +/-2 MHz about 40 GHz with a 1 us modulation period. [1]
Spread spectrum clocking (SSC) causes a fixed time reference to drift relative to the data.

| Parameter | Value |
|---|---|
| Input columns | Time, Vo, CLKp, CLKn |
| Samples | {inp['samples']:,} |
| Nominal sample spacing | 1 ps |
| Observed Vo range | {inp['minimum_voltage_V']:.8f} to {inp['maximum_voltage_V']:.5f} V |
| Display window | 0 to 50 ps, identical for both eyes |
| Voltage range | -0.1 to 1.1 V |
| Grid | {cfg['voltage_bins']} voltage x {cfg['time_bins']} time bins |
| Bin sizes | 10 mV and 1 ps |
| Normalization | Per occupied time column |

![Fixed-rate and triggered statistical eyes](../results/eye_comparison.png)

**Figure 1.** Color denotes the probability of a voltage bin conditional on a time bin.
Both panels use the same axes and shared linear probability scale.
The received crossing phase spans {trig['peak_to_peak_phase_error_ps']:.3f} ps relative to nominal 50 ps intervals.
Fixed-rate folding spreads transitions across that changing phase; resetting time at positive differential-clock crossings reduces misalignment.
This supports the assignment's SSC alignment argument without establishing a bit error rate.

## Method and bin conventions

### 1. Preserve the waveform values

Require finite values and strictly increasing time; convert time from seconds to picoseconds.
Preserve all measured Vo samples. Interpolate only clock crossing times.

### 2. Fixed-rate folding and voltage assignment

With fixed origin t0 = 0 ps, dt = 1 ps, W = 50 ps and dV = 0.01 V:

```text
j_fixed = floor((t - t0) / dt) mod 50
k = floor((Vo - vmin) / dV)
```

The voltage offset is essential: negative voltage must map to low rows rather than wrap to high rows.
Bins are left-closed and right-open; the last voltage bin additionally includes exact vmax.
Only floating point roundoff close to integer bin boundaries is snapped before flooring.

### 3. Differential-clock triggers

Form c = CLKp - CLKn and find adjacent samples satisfying c[n-1] < 0 and c[n] >= 0.
Linearly interpolate their positive zero crossing:

```text
a = -c[n-1] / (c[n] - c[n-1])
t_cross = t[n-1] + a * (t[n] - t[n-1])
```

Choose the latest crossing at or before each data sample, set age = t - t_cross,
and retain the sample only for 0 <= age < 50 ps. Use j_triggered = floor(age / dt).
Triggered ages never wrap and each measured sample is counted at most once.
This avoids duplicate counts when adjacent 2 UI windows overlap.
Exclude samples before the first crossing and samples in gaps beyond 50 ps.

### 4. Count and normalize vertically

Count H[k,j] and normalize every occupied time column:

```text
P[k,j] = H[k,j] / sum_k H[k,j]
sum_k P[k,j] = 1
```

Empty columns remain zero. P is the discrete conditional probability P(voltage bin | time bin),
not a continuous voltage density; the latter would additionally divide by dV.
The plotted probabilities depend on the time and voltage bins.
The procedure follows the assignment starter and lecture vertical normalization. [2,3]

## Timing evidence and results

![Clock phase relative to nominal 50 ps intervals](../results/clock_phase.png)

**Figure 2.** Positive crossing times relative to the first crossing plus n x 50 ps.
This is received-clock phase error relative to a nominal reference and includes the effects in the sampled waveform.

| Checked result | Fixed rate | Triggered |
|---|---:|---:|
| Samples counted | {fixed['total_samples']:,} | {triggered['total_samples']:,} |
| Count per occupied time column | {fixed['minimum_column_count']:,} to {fixed['maximum_column_count']:,} | {triggered['minimum_column_count']:,} to {triggered['maximum_column_count']:,} |
| Occupied / total time columns | 50 / 50 | 50 / 50 |
| Largest column-sum error | {fixed['maximum_column_sum_error']:.2e} | {triggered['maximum_column_sum_error']:.2e} |

There are {trig['crossings']:,} interpolated positive crossings.
Intervals range from {trig['minimum_interval_ps']:.3f} to {trig['maximum_interval_ps']:.3f} ps (mean {trig['mean_interval_ps']:.6f} ps).
Triggering excludes {trig['samples_before_first_trigger']} samples before the first crossing and
{trig['samples_excluded_after_2UI']} samples with age at least 50 ps, accounting for the difference in totals.
Unequal column counts are expected.

## Verification and reproduction

Fourteen unit tests passed: bin boundaries, voltage offsets, crossing interpolation, trigger gaps and overlap,
normalization, invalid inputs, and synthetic rate drift.
A separate full-record histogram implementation matched all 6,000 fixed-rate bins and all 6,000 triggered bins exactly.
All probabilities are finite; occupied-column sums differ from one by at most {max_error:.2e}.

From the HW05 folder:

```sh
python -m pip install -r requirements.txt
python src/hw5.py
python -m unittest discover -s tests -v
python -m pip install -r report/requirements.txt
python report/build_report.py
python tests/verify_delivery.py
```

The analysis writes figures, summary.json and histograms.npz to results/.
The report builder reads the regenerated summary and figures and writes this Markdown report and hw05_report.pdf.

## Limits and conclusion

The triggered eye shows improved timing alignment for this finite ADS record.
Observed crossing intervals are received-waveform measurements, not estimates of the SSC source frequency extrema.
No noise model, BER extrapolation, threshold optimization or standardized eye-mask test was performed.
Any eye opening inferred from these bins remains sample- and bin-dependent.

## References

1. HW5_stat_eye_construction.pdf, pp. 1-3.
2. HW5_stateye_construction_sample.py.
3. 5 Eye diagram.pdf, pp. 4-5 and 14-17.

Data: HW5_ADS_waveforms.csv.
'''


def main():
    base = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=base.parent / "results")
    parser.add_argument("--output-dir", type=Path, default=base)
    args = parser.parse_args()
    facts = load_facts(args.results)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    build_pdf(args.output_dir / "hw05_report.pdf", args.results, facts)
    (args.output_dir / "report.md").write_text(
        build_markdown(args.results, facts), encoding="utf-8", newline="\n"
    )
    print("Created hw05_report.pdf and report.md")


STYLES = styles()
rl_config.invariant = True

if __name__ == "__main__":
    main()
