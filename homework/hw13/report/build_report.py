"""Build the anonymous HW13 report from recorded numerical results.

Run without arguments to rebuild the PDF from the editable report.md.
Use --refresh to regenerate report.md from verified saved results first.
The builder does not rerun HW12 and does not modify source data or results.
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
MARGIN = 44
WIDTH = PAGE_WIDTH - 2 * MARGIN
INK = colors.HexColor("#193044")
ACCENT = colors.HexColor("#087B8A")
MUTED = colors.HexColor("#506473")
LINE = colors.HexColor("#D4E1E6")
PALE = colors.HexColor("#F0F6F7")


def styles():
    definitions = {
        "Body": dict(fontName="Helvetica", fontSize=9.5, leading=13,
                     textColor=INK, spaceAfter=8),
        "Title": dict(fontName="Helvetica-Bold", fontSize=25, leading=29,
                      textColor=INK, spaceAfter=13),
        "Heading": dict(fontName="Helvetica-Bold", fontSize=16, leading=21,
                        textColor=INK, spaceAfter=10, keepWithNext=True),
        "Subheading": dict(fontName="Helvetica-Bold", fontSize=10.7, leading=14,
                           textColor=ACCENT, spaceBefore=6, spaceAfter=7,
                           keepWithNext=True),
        "Caption": dict(fontName="Helvetica", fontSize=8.3, leading=11,
                        textColor=MUTED, spaceBefore=4, spaceAfter=10),
        "TableBody": dict(fontName="Helvetica", fontSize=8.1, leading=10.7,
                          textColor=INK),
        "TableHead": dict(fontName="Helvetica-Bold", fontSize=8.1, leading=10.7,
                          textColor=colors.white),
        "Kicker": dict(fontName="Helvetica-Bold", fontSize=9, leading=12,
                       textColor=ACCENT, spaceAfter=10),
        "Code": dict(fontName="Courier", fontSize=8.0, leading=10.7,
                     textColor=INK, backColor=PALE, borderPadding=7,
                     spaceBefore=4, spaceAfter=10),
    }
    return {key: ParagraphStyle(key, **value) for key, value in definitions.items()}


STYLE = styles()


def rich(text):
    value = html.escape(str(text), quote=False)
    value = re.sub(r"`([^`]+)`", r'<font name="Courier" size="8.2">\1</font>', value)
    value = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", value)
    return re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", value)


def paragraph(text, style="Body"):
    return Paragraph(rich(text), STYLE[style])


def table(rows):
    count = len(rows[0])
    if count == 2:
        widths = [WIDTH * .38, WIDTH * .62]
    elif count == 3:
        widths = [WIDTH * .38, WIDTH * .31, WIDTH * .31]
    elif count == 4:
        widths = [WIDTH * .28, WIDTH * .24, WIDTH * .24, WIDTH * .24]
    else:
        widths = [WIDTH / count] * count
    cells = [[paragraph(cell, "TableHead" if i == 0 else "TableBody")
              for cell in row] for i, row in enumerate(rows)]
    result = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LINEBELOW", (0, -1), (-1, -1), .6, LINE),
    ]))
    return result


def figure(relative_path, caption):
    path = (HERE / relative_path).resolve()
    if path.parent != RESULTS.resolve():
        raise ValueError("Report figures must be saved numerical result figures")
    width, height = ImageReader(str(path)).getSize()
    limits = {
        "primary_cdf.png": 380,
        "primary_bathtub.png": 380,
        "primary_pdf.png": 285,
        "primary_case3_fit.png": 215,
        "primary_case4_fit.png": 300,
        "sensitivity.png": 235,
    }
    scale = min(WIDTH / width, limits.get(path.name, 275) / height)
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
    canvas.drawString(MARGIN, PAGE_HEIGHT - 25, "SIGNAL INTEGRITY")
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawRightString(PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 25, "HW13 / BER EXTRAPOLATION")
    canvas.line(MARGIN, 33, PAGE_WIDTH - MARGIN, 33)
    canvas.drawString(MARGIN, 21, "Anonymous numerical report")
    canvas.drawRightString(PAGE_WIDTH - MARGIN, 21, str(document.page))
    canvas.restoreState()


def build(source, output):
    rl_config.invariant = True
    doc = BaseDocTemplate(
        str(output), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=51, bottomMargin=45,
        title="HW13 - BER Extrapolation", author="",
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
    verification = json.loads((RESULTS / "verification.json").read_text(encoding="utf-8"))
    tests = json.loads((RESULTS / "test_results.json").read_text(encoding="utf-8"))
    if verification.get("passed") is not True or tests.get("passed") is not True:
        raise ValueError("Recorded independent verification and tests must pass")
    if verification.get("production_module_imported") is not False:
        raise ValueError("Saved-artifact verification must remain independent")
    return verification, tests


def report_markdown(summary):
    verification, tests = load_evidence()
    empirical_check = next(item for item in verification["checks"]
                           if item["check"] == "all primary and corrected legacy empirical exports")
    model_check = next(item for item in verification["checks"]
                       if item["check"] == "all saved model minima, stationarity, roots and widths with math.erfc")
    curve_check = next(item for item in verification["checks"]
                       if item["check"] == "every saved bathtub ordinate")
    curve_ordinates = (curve_check["sampled_curves"] + curve_check["model_curves"]) * curve_check["time_points"]
    if summary.get("schema_version") != 1:
        raise ValueError("Review report narrative for the changed result schema")
    params = summary["parameters"]
    for key, expected in {"ui_ps": 25, "dt_ps": .125, "rho_transition": .5,
                          "tail_min": 3e-4, "tail_max": 3e-2}.items():
        if not math.isclose(params[key], expected, rel_tol=1e-12):
            raise ValueError(f"Review report narrative for changed {key}")
    if params["ideal_boundaries_ps"] != [12.5, 37.5]:
        raise ValueError("Review report narrative for changed ideal boundaries")
    if summary["source"]["bits_per_case"] < 2**20:
        raise ValueError("The published report requires at least 2**20 evaluated bits")
    c3 = summary["cases"]["case3_rj"]
    c4 = summary["cases"]["case4_sj_rj"]
    for case in (c3, c4):
        if case["at_bound"] or case["converged_starts"] != case["starts"]:
            raise ValueError("Review convergence before publishing tail parameters")
    def width(case, target, unit="ps"):
        item = case["metrics"]["eye_widths"][target]
        return fnum(item[f"width_{unit}"]) if item["exists"] else "No opening"
    c3eye = c3["metrics"]["eye_widths"]["1e-12"]
    parameter_rows = [["Tail-model parameter", "Case 3 Rx", "Case 4 Rx"]]
    for label, key in [("mu1 (ps)", "mu1_ps"), ("sigma1 (ps)", "sigma1_ps"),
                       ("mu2 (ps)", "mu2_ps"), ("sigma2 (ps)", "sigma2_ps"),
                       ("tau1", "weight1")]:
        parameter_rows.append([label, fnum(c3["parameters"][key]), fnum(c4["parameters"][key])])
    parameter_rows.append(["tau2 = 1 - tau1", fnum(c3["weight2"]), fnum(c4["weight2"])])
    stats_rows = [["Distribution", "Slice observations", "Zero bins", "Smallest nonzero PMF"]]
    case_keys = ("case1_no_tx_jitter", "case2_sj", "case3_rj", "case4_sj_rj")
    for number, key in enumerate(case_keys, 1):
        for side in ("tx", "rx"):
            record = summary["sample_statistics"][f"{key}_{side}"]
            stats_rows.append([f"Case {number} {side.upper()}", f"{record['observations']:,.0f}",
                               record["zero_bins"], fnum(record["smallest_nonzero_pmf"])])
    eye_rows = [["Target BER", "Case 3 EW (ps)", "Case 3 EW (UI)", "Case 4 EW (ps)"]]
    for target in ("1e-06", "1e-09", "1e-12", "1e-15"):
        eye_rows.append([target, width(c3, target), width(c3, target, "ui"), width(c4, target)])
    fit_rows = [["Fit diagnostic", "Case 3 Rx", "Case 4 Rx"],
                ["Selected tail locations", c3["selected_bins"], c4["selected_bins"]],
                ["CDF/SF log10 RMSE", fnum(c3["tail_log10_rmse"]), fnum(c4["tail_log10_rmse"])],
                ["PDF-tail log10 RMSE", fnum(c3["pdf_tail_log10_rmse"]), fnum(c4["pdf_tail_log10_rmse"])],
                ["Converged starts", f"{c3['converged_starts']} / {c3['starts']}",
                 f"{c4['converged_starts']} / {c4['starts']}"]]
    sensitivity_rows = [["Tail min", "Tail max", "Case 3 min BER", "Case 3 EW (ps)", "Case 4 min BER"]]
    sensitivity_index = {(row["case"], row["tail_min"], row["tail_max"]): row
                         for row in summary["sensitivity"]}
    for low in (1e-4, 3e-4, 1e-3):
        for high in (.01, .03, .1):
            r3 = sensitivity_index[("case3_rj", low, high)]
            r4 = sensitivity_index[("case4_sj_rj", low, high)]
            sensitivity_rows.append([fnum(low), fnum(high), fnum(r3["metrics"]["minimum_ber"]),
                                     width(r3, "1e-12"), fnum(r4["metrics"]["minimum_ber"])])
    case3_windows = [row for row in summary["sensitivity"] if row["case"] == "case3_rj"]
    case3_min_range = [min(row["metrics"]["minimum_ber"] for row in case3_windows),
                       max(row["metrics"]["minimum_ber"] for row in case3_windows)]
    case3_max_width = max(row["metrics"]["eye_widths"]["1e-12"]["width_ps"] for row in case3_windows)
    legacy = summary["legacy"]
    legacy_rows = [["Comparison / Case 3", "Minimum BER", "EW at 1e-12 (ps)"]]
    legacy_rows.append(["Verified HW12 input, tail fit", fnum(c3["metrics"]["minimum_ber"]), width(c3, "1e-12")])
    legacy_refit = legacy["refitted_cases"]["case3_rj"]
    legacy_rows.append(["Corrected legacy CSV, tail fit", fnum(legacy_refit["metrics"]["minimum_ber"]), width(legacy_refit, "1e-12")])
    for key, label in [("archive_case3", "Archive variant, supplied parameters"),
                       ("standalone_case3", "Adjacent variant, supplied parameters")]:
        record = legacy["supplied_parameter_predictions"][key]
        legacy_rows.append([label, fnum(record["metrics"]["minimum_ber"]), width(record, "1e-12")])
    return f"""# Jitter fitting and BER extrapolation

> HW13 - Signal integrity numerical analysis

The verified full-length HW12 data give a two-Gaussian tail-model minimum BER of
**{fnum(c3['metrics']['minimum_ber'])} for Case 3 Rx** and
**{fnum(c4['metrics']['minimum_ber'])} for Case 4 Rx**. At BER = 1e-12, the nominal
Case 3 eye width is **{width(c3, '1e-12')} ps ({width(c3, '1e-12', 'ui')} UI)**;
Case 4 has no opening. These are extrapolated model predictions, not measured
error rates. The sensitivity analysis on the final page shows why the narrow
Case 3 opening must be interpreted cautiously.

## Scope and provenance

All eight sampled CDFs and bathtubs are calculated for four cases at Tx and Rx.
The two-Gaussian extrapolation then focuses on Case 3 and Case 4 Rx, as requested
in the assignment [1]. The primary input is the saved, independently verified
HW12 horizontal eye-slice distribution from **1,048,576 evaluated bits per case**.
That full channel simulation was reused; **HW12 was not rerun for this analysis**.

| Physical / numerical quantity | Setting |
| --- | --- |
| Signaling / UI | 40 Gb/s NRZ / 25 ps |
| Cases 1-4 | No Tx jitter; SJ; RJ; SJ + RJ |
| SJ / RJ | 3 ps peak at 1 MHz / zero-mean Gaussian, sigma = 1 ps |
| Source / load | 20 ohm / 85 ohm, each with 0.3 pF shunt capacitance |
| Rise and fall duration | 5 ps |
| Statistical time grid | 200 bins per UI; dt = 0.125 ps |
| Transition density | rhoT = 0.5 |
| Adjacent ideal boundaries | t_i = 12.5 ps; t_(i+1) = 37.5 ps |

### Probability definitions and timing convention

For normalized bin probability p[k], the PMF is dimensionless and sums to one.
The density f[k] = p[k]/dt has units 1/ps and integrates to one. The empirical
CDF at a bin center uses half that cell's mass; interpolation assumes uniform
density inside each centered 0.125 ps cell. This fixes the time origin and keeps
the CDF and density consistent, rather than treating endpoint-inclusive samples
as a different physical grid.

```text
F(u) = P(jitter <= u);  S(u) = P(jitter > u)
BER(t) = rhoT * [ S(t - 12.5 ps) + F(t - 37.5 ps) ]
f2G(u) = tau1 * phi((u-mu1)/sigma1)/sigma1
       + tau2 * phi((u-mu2)/sigma2)/sigma2
F2G(u) = tau1 * Phi((u-mu1)/sigma1)
       + tau2 * Phi((u-mu2)/sigma2);  tau2 = 1 - tau1
```

Here phi and Phi denote the standard-normal PDF and CDF. The timing-error
coordinate u is centered at an ideal transition; absolute sampling time t lies
between the two boundaries. Gaussian SF and CDF tails are evaluated directly
with stable log-domain normal-CDF/SF functions and log-sum addition, including the SF instead
of subtracting a nearly unit CDF from one [2].

### Reproduce from the repository root

```sh
python homework/hw13/verify_inputs.py
python homework/hw13/src/hw13.py
python -m unittest discover -s homework/hw13/tests -v
python homework/hw13/tests/verify_results.py
python -m pip install -r homework/hw13/report/requirements.txt
python homework/hw13/report/build_report.py --refresh
```

<!-- pagebreak -->

## 1. Sampled distributions: all eight cases

![Figure 1. Sampled CDFs for every Tx/Rx distribution. Color identifies the case; dashed and solid lines distinguish Tx and Rx. Zero CDF values are masked only for logarithmic display. No low-probability floor is added to the analysis.](../results/primary_cdf.png)

![Figure 2. The same eight distributions as PDFs in 1/ps. These are normalized eye-slice densities, not counts or probabilities per bin. The narrow ideal-Tx peak is a sampled, quantized distribution.](../results/primary_pdf.png)

The CDF is an accumulated PMF, not an accumulated PDF without a dt factor.
Case 1 Rx spreads despite its ideal Tx because channel ISI causes DDJ. Adding
SJ widens the distributions; RJ creates sampled tails. The central multimodal
Rx structure need not have a Gaussian shape and is excluded from the tail fit.
The eye slice also weights voltage occupancy and edge slope; it is not an
equal-weight histogram of independently detected crossing events [2, pp. 2-4].

<!-- pagebreak -->

## 2. Sampled bathtub curves and finite evidence

![Figure 3. All eight non-extrapolated bathtub curves, formed from the empirical CDF and survival function at neighboring transitions. Curves cover the 2-UI plotting interval; eye-width metrics use the eye between 12.5 and 37.5 ps.](../results/primary_bathtub.png)

The late left transition contributes S(t - t_i); the early right transition
contributes F(t - t_(i+1)). Moving the sampler toward either boundary increases its
corresponding tail. This is the two-transition construction in the lecture [2,
pp. 10-13]. The left PDF tail governs the right side of the bathtub, and the
right PDF tail governs its left side [1, p. 2].

{markdown_table(stats_rows)}

All eight PMFs sum to one and all eight PDFs integrate to one on the recorded
grid. The receiver slice counts differ from the bit count: Case 3 and Case 4
contain {summary['sample_statistics']['case3_rj_rx']['observations']:,.0f} and
{summary['sample_statistics']['case4_sj_rj_rx']['observations']:,.0f} threshold-voltage
occupancies. Their smallest nonzero PMF increments are approximately 7.4e-6.
Thus a zero empirical bathtub value is absence of observed slice support,
**not evidence of zero physical BER or direct observation of BER = 1e-12**.
The zeros remain exact in numerical arrays and are omitted only when plotting
on logarithmic axes.

<!-- pagebreak -->

## 3. Tail model and Case 3 Rx result

The two-Gaussian model has separate means, standard deviations and a single
independent mixing weight. Only tail locations with one-sided empirical
CDF or SF in **[3e-4, 3e-2]** are selected. Joint nonlinear least squares minimizes
unweighted residuals in log10(CDF) and log10(SF); it does not fit the central DJ
peak. Six starting guesses converge, and neither reported solution reaches a
parameter bound. This is a tail-fit procedure, **not maximum likelihood**.

The fit evaluates a selected bin's CDF at its upper edge, u[k] + dt/2, and its
SF at the lower edge, u[k] - dt/2, so each includes the full selected cell mass.
This half-bin correction is distinct from the center-interpolated plotted CDF.
The bounds are mu1 in [-8, 0] ps, mu2 in [0, 8] ps, both sigmas in [0.2, 5] ps,
and tau1 in [0.01, 0.99].

{markdown_table(parameter_rows)}

![Figure 4. Case 3 Rx: sampled PDF versus the two-Gaussian model and sampled bathtub versus extrapolation. Highlighted PDF points identify the selected tail bins; the model extends below the empirical evidence scale.](../results/primary_case3_fit.png)

{markdown_table(fit_rows)}

Case 3 fits left-center coordinates from {fnum(c3['left_center_range_ps'][0])} to
{fnum(c3['left_center_range_ps'][1])} ps and right-center coordinates from
{fnum(c3['right_center_range_ps'][0])} to {fnum(c3['right_center_range_ps'][1])} ps.
The selected probability threshold corresponds to approximately
{c3['minimum_tail_observations']:.1f} cumulative slice observations at its low end.
The nominal minimum is **{fnum(c3['metrics']['minimum_ber'])}** at
**{fnum(c3['metrics']['optimal_time_ps'])} ps**. At 1e-12, the crossing times are
{fnum(c3eye['left_ps'])} and {fnum(c3eye['right_ps'])} ps, giving
**{width(c3, '1e-12')} ps** of opening. The lower-than-empirical values follow the
Gaussian tail assumption, not new simulated errors.

<!-- pagebreak -->

## 4. Case 4 Rx and low-BER eye widths

![Figure 5. Case 4 Rx PDF-tail fit and extrapolated bathtub. The predicted minimum remains above 1e-9, so there is no modeled eye opening at 1e-9, 1e-12 or 1e-15.](../results/primary_case4_fit.png)

Case 4 fits left-center coordinates from {fnum(c4['left_center_range_ps'][0])} to
{fnum(c4['left_center_range_ps'][1])} ps and right-center coordinates from
{fnum(c4['right_center_range_ps'][0])} to {fnum(c4['right_center_range_ps'][1])} ps.
The nominal minimum is **{fnum(c4['metrics']['minimum_ber'])}** at
**{fnum(c4['metrics']['optimal_time_ps'])} ps**. An eye width exists only when the
target BER is at least the curve minimum; it is the distance between the two
target-BER crossings inside this eye.

{markdown_table(eye_rows)}

### Corrections and comparison with supplied code

The executable analysis preserves zero probabilities, keeps dt in ps throughout,
uses the physical 0.125 ps grid, and computes extreme Gaussian SF/CDF tails
directly. It resolves data files relative to the project, uses valid positive
logarithmic display limits, and calculates numerical minima and crossing roots.
The incomplete starter is retained as instructions. Both supplied legacy
variants contain an inherited 1e-15 display offset and inconsistent endpoint
coordinates; those artifacts are not accepted as fresh tail observations.

{markdown_table(legacy_rows)}

The legacy comparison subtracts the known display offset and restores the
physical grid before refitting. Its original observation counts and exact
generating run are unavailable. The two supplied parameter sets also predict
different Case 3 conclusions. They are references, not a unique numerical answer
and not the primary input used above.

<!-- pagebreak -->

## 5. Sensitivity and limitations

![Figure 6. Predicted minimum BER for nine tail-window choices per receiver case. The table also gives Case 3 eye widths at 1e-12. These model responses do not represent confidence intervals.](../results/sensitivity.png)

{markdown_table(sensitivity_rows)}

Across all nine windows, Case 3 minimum BER ranges from
{fnum(case3_min_range[0])} to {fnum(case3_min_range[1])}, and its 1e-12 eye width
ranges from no opening to {fnum(case3_max_width)} ps. Including
probabilities down to 1e-4 leaves only about 14 cumulative slice observations
and can remove the opening; the narrow high-tail window can predict a larger
opening. **These are fit-window sensitivity results, not confidence intervals.**
The cumulative points are correlated and their RMSE is descriptive. Case 4 has
no 1e-12 opening in all nine tested windows.

The Gaussian-tail and two-component assumptions, finite time/voltage bins,
one-UI folding of extreme jitter, and slope-weighted eye slices limit physical
accuracy. This report contains no low-BER measurement, hardware validation or
new ADS ChannelSim result. Saved-artifact and numerical checks validate the
stated calculation, not the model assumptions. All **{tests['tests_run']} tests passed**.
An independent verifier, without importing the production module, checked
{empirical_check['bins']:,} empirical bins, {model_check['models']} model sets and
{curve_ordinates:,} bathtub ordinates using direct math.erfc tails, cell-mass
integration and separate minimum/root solvers. Maximum relative model-curve
disagreement was {fnum(curve_check['maximum_model_relative_error'])}.

### References

[1] Homework: Jitter fitting and BER extrapolation, supplied assignment,
pp. 1-2 (required tasks, finite-pattern note, zero-floor rule and tail mapping).
[2] Jitter Basics (2), course handout 13, pp. 2-4 (eye slice and CDF),
pp. 10-13 (two-boundary BER and eye width), pp. 15-17 (tail extrapolation and
two-Gaussian mixture). The optional ADS example is not presented as a new result.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="regenerate report.md from verified saved results")
    args = parser.parse_args()
    load_evidence()
    source = HERE / "report.md"
    if args.refresh or not source.exists():
        summary = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
        source.write_text(report_markdown(summary), encoding="utf-8")
    build(source, HERE / "report.pdf")


if __name__ == "__main__":
    main()
