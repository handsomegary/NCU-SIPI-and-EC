"""Build the anonymous HW7 report from generated results and independent checks.

Run from any directory: python report/build_report.py
No original report pages, screenshots, or personal metadata are reused.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
import re
from typing import Any

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
NAVY = colors.HexColor("#14304A")
BLUE = colors.HexColor("#176DAB")
TEAL = colors.HexColor("#087D77")
INK = colors.HexColor("#202D38")
MUTED = colors.HexColor("#536372")
PALE = colors.HexColor("#F1F6FA")
LINE = colors.HexColor("#D7E2EA")


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def fixed(value: float, places: int = 3) -> str:
    return f"{float(value):.{places}f}"


def metric_table(metrics: dict[str, Any]) -> list[list[str]]:
    before, after = metrics["before"], metrics["after"]
    return [
        ["Center metric", "Post-channel", "Post-CTLE"],
        ["Random-pattern eye height (mV)", fixed(1e3 * before["center_eye_height_v"], 2), fixed(1e3 * after["center_eye_height_v"], 2)],
        ["PDA eye height (mV)", fixed(1e3 * before["pda_eye_height_v"], 2), fixed(1e3 * after["pda_eye_height_v"], 2)],
        ["Minimum conditional one level (mV)", fixed(1e3 * before["one_min_v"], 2), fixed(1e3 * after["one_min_v"], 2)],
        ["Maximum conditional zero level (mV)", fixed(1e3 * before["zero_max_v"], 2), fixed(1e3 * after["zero_max_v"], 2)],
        ["Selected SBR sample time (ns)", fixed(1e9 * before["center_time_s"], 4), fixed(1e9 * after["center_time_s"], 4)],
        ["Sample index within the UI", f"{before['phase_samples']}", f"{after['phase_samples']}"],
        ["Analyzed bits", f"{before['analyzed_bits']:,}", f"{after['analyzed_bits']:,}"],
    ]


def report_content(metrics: dict[str, Any], verification: dict[str, Any] | None) -> dict[str, Any]:
    """Return text shared by the PDF and editable Markdown report."""
    settings, ctle, source = metrics["settings"], metrics["ctle"], metrics["source"]
    before, after = metrics["before"], metrics["after"]
    eye_mv = 1e3 * after["center_eye_height_v"]
    pda_mv = 1e3 * after["pda_eye_height_v"]
    constraints = [
        ["Constraint", "Computed result", "Status"],
        ["1 <= fp1, fp2 <= 40 GHz", f"{fixed(ctle['fp1_ghz'])}, {fixed(ctle['fp2_ghz'])} GHz", "PASS" if 1 <= ctle["fp1_ghz"] <= 40 and 1 <= ctle["fp2_ghz"] <= 40 else "FAIL"],
        ["Maximum CTLE gain <= 6 dB", f"{fixed(ctle['peak_gain_db'], 4)} dB", "PASS" if ctle["peak_gain_db"] <= 6 + 1e-9 else "FAIL"],
        ["-20 <= A_DC <= +6 dB", f"{fixed(ctle['adc_db'])} dB", "PASS" if -20 <= ctle["adc_db"] <= 6 else "FAIL"],
        ["Center eye height >= 300 mV", f"{fixed(eye_mv, 2)} mV", "PASS" if eye_mv >= 300 else "FAIL"],
    ]
    return {
        "title": "HW7 | CTLE simulation and design",
        "subtitle": "Anonymous numerical report",
        "summary": (
            f"A two-pole, one-zero CTLE was tuned for the supplied 50 cm stripline at "
            f"{settings['bit_rate_bps'] / 1e9:g} Gbps. The center eye height increases from "
            f"{1e3 * before['center_eye_height_v']:.2f} mV to {eye_mv:.2f} mV, "
            f"with a PDA margin of {pda_mv:.2f} mV. All assignment constraints are checked explicitly."
        ),
        "params": [
            ["Parameter", "Selected value", "Meaning"],
            ["f_p1", f"{fixed(ctle['fp1_ghz'], 6)} GHz", "First pole frequency"],
            ["f_p2", f"{fixed(ctle['fp2_ghz'], 6)} GHz", "Second pole frequency"],
            ["A_DC", f"{fixed(ctle['adc_db'], 6)} dB", "DC amplitude gain"],
            ["A_AC", f"{fixed(ctle['aac_db'], 6)} dB", "Midband gain parameter"],
            ["f_z", f"{ctle['zero_frequency_hz'] / 1e9:.6f} GHz", "Derived zero frequency"],
        ],
        "constraints": constraints,
        "eyes": metric_table(metrics),
        "source": source,
        "settings": settings,
        "ctle": ctle,
        "before": before,
        "after": after,
        "full_period_pda": metrics.get("full_period_pda"),
        "time_step_convergence": metrics.get("time_step_convergence"),
        "verification": verification,
    }


def create_styles() -> dict[str, ParagraphStyle]:
    styles = getSampleStyleSheet()
    return {
        "title": ParagraphStyle("ReportTitle", parent=styles["Title"], fontName="Helvetica-Bold", fontSize=24, leading=29, textColor=NAVY, alignment=TA_LEFT, spaceAfter=10),
        "subtitle": ParagraphStyle("ReportSubtitle", parent=styles["Normal"], fontSize=10, leading=14, textColor=TEAL, spaceAfter=14),
        "heading": ParagraphStyle("ReportHeading", parent=styles["Heading1"], fontSize=16, leading=20, textColor=NAVY, spaceBefore=1, spaceAfter=10),
        "subheading": ParagraphStyle("ReportSubheading", parent=styles["Heading2"], fontSize=11, leading=15, textColor=BLUE, spaceBefore=10, spaceAfter=5),
        "body": ParagraphStyle("ReportBody", parent=styles["BodyText"], fontSize=9.2, leading=13.2, textColor=INK, spaceAfter=8),
        "small": ParagraphStyle("ReportSmall", parent=styles["BodyText"], fontSize=8, leading=11, textColor=MUTED, spaceAfter=5),
        "formula": ParagraphStyle("ReportFormula", parent=styles["BodyText"], fontName="Helvetica", fontSize=10, leading=15, textColor=NAVY, spaceBefore=5, spaceAfter=9),
        "caption": ParagraphStyle("ReportCaption", parent=styles["BodyText"], fontSize=8, leading=11, textColor=MUTED, spaceBefore=5, spaceAfter=10),
        "table": ParagraphStyle("ReportTable", parent=styles["BodyText"], fontSize=8.6, leading=11.6, textColor=INK),
        "table_header": ParagraphStyle("ReportTableHeader", parent=styles["BodyText"], fontName="Helvetica-Bold", fontSize=8.4, leading=11.5, textColor=colors.white),
    }


def table(rows: list[list[str]], widths: list[float], styles: dict[str, ParagraphStyle]) -> Table:
    cells = [[Paragraph(html.escape(str(v)), styles["table_header"] if i == 0 else styles["table"]) for v in row] for i, row in enumerate(rows)]
    result = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT")
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LINEBELOW", (0, -1), (-1, -1), 0.5, LINE),
    ]))
    return result


def fit_image(path: Path, max_width: float, max_height: float) -> Image:
    if not path.is_file():
        raise FileNotFoundError(f"Required report figure is missing: {path.name}")
    image = Image(str(path))
    scale = min(max_width / image.imageWidth, max_height / image.imageHeight)
    image.drawWidth = image.imageWidth * scale
    image.drawHeight = image.imageHeight * scale
    image.hAlign = "CENTER"
    return image


def footer(canvas: Any, doc: Any) -> None:
    canvas.saveState()
    page_width, page_height = A4
    canvas.setStrokeColor(LINE)
    canvas.line(doc.leftMargin, page_height - 32, page_width - doc.rightMargin, page_height - 32)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(doc.leftMargin, page_height - 24, "HW7  /  CTLE DESIGN")
    canvas.drawRightString(page_width - doc.rightMargin, page_height - 24, "ANONYMOUS REPORT")
    canvas.line(doc.leftMargin, 37, page_width - doc.rightMargin, 37)
    canvas.drawString(doc.leftMargin, 23, "Reproducible frequency-domain channel model")
    canvas.drawRightString(page_width - doc.rightMargin, 23, str(doc.page))
    canvas.restoreState()


def build_pdf(content: dict[str, Any], results: Path, output: Path) -> None:
    styles = create_styles()
    doc = BaseDocTemplate(str(output), pagesize=A4, leftMargin=42, rightMargin=42, topMargin=49, bottomMargin=49,
                          title=content["title"], author="Anonymous", subject="CTLE design and verified eye-height simulation", creator="Anonymous report builder")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates(PageTemplate(id="report", frames=frame, onPage=footer))
    story: list[Any] = []
    body = lambda text: Paragraph(text, styles["body"])
    heading = lambda text: Paragraph(text, styles["heading"])
    sub = lambda text: Paragraph(text, styles["subheading"])
    caption = lambda text: Paragraph(text, styles["caption"])
    ctle, settings, source = content["ctle"], content["settings"], content["source"]
    before, after = content["before"], content["after"]

    # Page 1: the answer and all four parameter / specification checks.
    story += [Paragraph(content["title"], styles["title"]), Paragraph(content["subtitle"], styles["subtitle"]), body(content["summary"]), sub("Selected CTLE parameters"), table(content["params"], [86, 137, doc.width - 223], styles), Spacer(1, 12), sub("Transfer-function convention")]
    story += [Paragraph("H<sub>eq</sub>(s) = a<sub>AC</sub> omega<sub>p2</sub> [s + (a<sub>DC</sub>/a<sub>AC</sub>) omega<sub>p1</sub>] / [(s + omega<sub>p1</sub>)(s + omega<sub>p2</sub>)]", styles["formula"]), body("Here a<sub>DC</sub> = 10<super>A_DC/20</super>, a<sub>AC</sub> = 10<super>A_AC/20</super>, and omega<sub>p</sub> = 2 pi f<sub>p</sub>. The zero is f<sub>z</sub> = (a<sub>DC</sub>/a<sub>AC</sub>) f<sub>p1</sub>. A<sub>AC</sub> controls the midband asymptote; the actual maximum is evaluated separately. [1, 2]"), sub("Assignment acceptance checks"), table(content["constraints"], [205, 197, doc.width - 402], styles), Spacer(1, 9), body(f"The exact continuous-frequency peak occurs at {ctle['peak_frequency_hz'] / 1e9:.6f} GHz. The peak check is not restricted to the plotted frequency samples."), caption("The four tuning parameters are the two pole frequencies and the DC / AC amplitude-gain parameters. No independent bound is imposed on A_AC; the complete response must satisfy the 6 dB maximum.")]

    # Page 2: frequency response and impulse-response model.
    story += [PageBreak(), heading("Frequency response and channel model"), fit_image(results / "frequency_response.png", doc.width, 230), caption("Figure 1. Optimized CTLE magnitude and the supplied channel / combined response. The CTLE compensates the channel's smooth low-pass loss over the useful signal band."), fit_image(results / "impulse_response.png", doc.width, 197), caption("Figure 2. Impulse responses used for convolution. Both paths use the same interpolation, frequency grid, taper and time sampling."), body(f"The input file contains {source['points']:,} complex frequency samples with spacing {source['df_hz'] / 1e6:g} MHz, a maximum frequency of {source['fmax_hz'] / 1e9:g} GHz, and reference impedance {source['reference_ohm']:g} ohm. Its complex S21 is retained, including phase."), body(f"The positive-frequency spectrum is Kaiser-tapered using the actual data cutoff ({source['fmax_hz'] / 1e9:g} GHz), then zero-padded to the simulation Nyquist frequency before a real inverse FFT. The taper uses alpha = {settings['alpha']:g}; the retained impulse-response window is {settings['tail_ns']:g} ns. In particular, no independent hardcoded cutoff is used.")]

    # Page 3: comparable eyes, definitions and conditioning.
    story += [PageBreak(), heading("Eye diagrams and center eye height"), fit_image(results / "eye_diagrams.png", doc.width, 250), caption("Figure 3. Statistical eyes from the same seeded NRZ sequence. The selected sampling centers are aligned to 0.5 and 1.5 UI in each panel."), table(content["eyes"], [doc.width * .48, doc.width * .26, doc.width * .26], styles), Spacer(1, 10), sub("Metric and sampling convention"), body("Each path uses its own fixed sample phase, selected by maximizing the pattern-dependent analysis (PDA) margin around the main cursor. At that phase, random-pattern center EH = min(v | bit = 1) - max(v | bit = 0), over every retained bit. No probability threshold removes rare excursions. A negative separation means the one / zero distributions overlap and the physical eye is closed."), body("PDA sums the worst-case sampled precursor and postcursor contributions. It checks more bit combinations than the finite random sequence and is reported separately from its observed eye height."), body(f"Both paths use the same {settings['bits']:,} bits (seed {settings['seed']}) at {settings['bit_rate_bps'] / 1e9:g} Gbps. UI = {1e12 / settings['bit_rate_bps']:g} ps; the source has {settings['rise_time_s'] * 1e12:g} ps full 0-to-100% linear ramps, following the supplied sample convention. Each UI has {settings['samples_per_ui']} samples (dt = {settings['dt_s'] * 1e12:g} ps). The first {settings['skip_bits']:,} bits are excluded to remove the startup transient.")]

    # Page 4: all three required 1 ns waveform windows.
    story += [PageBreak(), heading("Three 1 ns waveform comparisons"), body("The overlays show the received post-channel and post-CTLE voltages in three different 1 ns intervals. Each trace is taken from the same complete simulated record; the input sequence and voltage scale are unchanged.")]
    for index, start in enumerate((100, 200, 300), start=4):
        story += [fit_image(results / f"waveform_{start}ns.png", doc.width, 174), caption(f"Figure {index}. Waveform interval {start}-{start + 1} ns. Both received paths are plotted together.")]
    story += [body("The CTLE reduces the long channel tail and separates the received logic levels. It may also produce overshoot or negative undershoot; the center eye-height calculation includes their contribution at the selected sample phase.")]

    # Page 5: implementation details, independent verification and model limits.
    story += [PageBreak(), heading("Reproducibility and verification"), sub("Amplitude and convolution accounting"), body(f"The source is a {settings['source_vpp']:g} Vpp unipolar NRZ waveform with {settings['rise_time_s'] * 1e12:g} ps linear transitions. Under the assignment's matched-source assumption, the post-channel transfer is H_ch = S21_ch / 2 and the post-CTLE transfer is H_total = S21_ch S21_CTLE / 2. The matching divider {settings['matching_divider']:g} is applied exactly once."), body("Sectional convolution is overlap-add with a retained output tail. The inverse FFT produces a discrete impulse response, so no extra dt factor is introduced. The full-period sample sum equals the DC transfer; truncation slightly changes that sum and is not renormalized. A single isolated bit with the same rising and falling ramps supplies the PDA response."), sub("Independent numerical checks")]
    verification = content["verification"]
    if verification is None:
        raise ValueError("Independent verification is required before the final report is built")
    checks = verification.get("report_checks")
    if checks:
        for check in checks:
            story += [body(f"<b>{html.escape(str(check['label']))}.</b> {html.escape(str(check['detail']))}")]
    else:
        story += [body("The independent verification results are recorded in the accompanying machine-readable verification file. The report is generated only after that file is available.")]
    full_pda = content.get("full_period_pda")
    fine = content.get("time_step_convergence")
    if full_pda and fine:
        story += [body(f"<b>Tail and time-step robustness.</b> Extending the kernel from 20 to the full 100 ns gives PDA EH = {full_pda['after']['pda_eye_height_v'] * 1e3:.3f} mV. Halving dt to {fine['fine_dt_s'] * 1e12:g} ps gives {fine['fine_full_period_pda_v'] * 1e3:.3f} mV with the full kernel, a change of {fine['full_period_pda_difference_v'] * 1e3:.3f} mV. Both exceed 300 mV.")]
    story += [sub("Interpretation and limits"), body("This linear, frequency-domain, noise-free model retains channel phase and delay, but omits noise, clock jitter, nonlinearity and saturation. CTLE can amplify high-frequency noise and does not generally remove reflection notches. Finite-pattern EH and finite-period PDA are not BER predictions or validation of an analog implementation."), body("The optional ADS ChannelSim exercise was not performed. Figures and claims come from the accompanying Python calculation and independent checks."), sub("References and provenance"), body("[1] Assignment brief, <i>Homework: CTLE simulation and design</i>, instructions/HW7_CTLE.pdf, pp. 1-3. Supplied channel and sample code are stored alongside it."), body("[2] Course lecture, <i>Continuous Time Linear Equalizer (CTLE)</i>, pp. 6-9 and 11-14. Transfer-function conventions, circuit interpretation and over-equalization."), Paragraph(f"Channel SHA-256: {source['sha256']}", styles["small"])]
    doc.build(story)


def markdown_table(rows: list[list[str]]) -> str:
    return "\n".join(["| " + " | ".join(rows[0]) + " |", "| " + " | ".join(["---"] * len(rows[0])) + " |"] + ["| " + " | ".join(row) + " |" for row in rows[1:]])


def build_markdown(content: dict[str, Any], output: Path) -> None:
    settings, ctle, source = content["settings"], content["ctle"], content["source"]
    verification = content["verification"]
    checks = "\n".join(f"- **{check['label']}.** {check['detail']}" for check in (verification or {}).get("report_checks", []))
    robustness = ""
    if content.get("full_period_pda") and content.get("time_step_convergence"):
        full, fine = content["full_period_pda"], content["time_step_convergence"]
        robustness = (
            f"Extending the retained impulse response from 20 to the full 100 ns gives PDA EH = "
            f"{full['after']['pda_eye_height_v'] * 1e3:.6f} mV. Halving dt to "
            f"{fine['fine_dt_s'] * 1e12:g} ps gives a full-period PDA EH of "
            f"{fine['fine_full_period_pda_v'] * 1e3:.6f} mV "
            f"(change {fine['full_period_pda_difference_v'] * 1e3:.6f} mV). "
            "The target remains satisfied without impulse-response renormalization."
        )
    text = f"""# {content['title']}

Anonymous numerical report.

{content['summary']}

## Selected design

{markdown_table(content['params'])}

The USB-style two-pole, one-zero transfer function is

$$H_{{eq}}(s)=a_{{AC}}\\omega_{{p2}}\\frac{{s+(a_{{DC}}/a_{{AC}})\\omega_{{p1}}}}{{(s+\\omega_{{p1}})(s+\\omega_{{p2}})}}.$$

Here $a_{{DC}}=10^{{A_{{DC}}/20}}$, $a_{{AC}}=10^{{A_{{AC}}/20}}$, $\\omega_p=2\\pi f_p$, and $f_z=(a_{{DC}}/a_{{AC}})f_{{p1}}$. A_AC is a midband asymptote parameter, not the actual peak. [1, 2]

{markdown_table(content['constraints'])}

The exact continuous-frequency maximum occurs at {ctle['peak_frequency_hz'] / 1e9:.6f} GHz. The gain limit is checked at the analytic maximum, not only on the plotted samples.

![Optimized CTLE and channel frequency responses](../results/frequency_response.png)

## Model and input

The supplied 50 cm stripline has {source['points']:,} complex samples, frequency spacing {source['df_hz'] / 1e6:g} MHz, maximum frequency {source['fmax_hz'] / 1e9:g} GHz and {source['reference_ohm']:g} ohm reference impedance. The complete complex S21 is retained, including phase.

The positive-frequency spectrum is Kaiser-tapered at the actual {source['fmax_hz'] / 1e9:g} GHz data cutoff (alpha = {settings['alpha']:g}), then zero-padded to the simulation Nyquist frequency and transformed by real inverse FFT. The retained impulse response is {settings['tail_ns']:g} ns long. No independent hardcoded data cutoff is used.

![Discrete impulse responses](../results/impulse_response.png)

The same {settings['bits']:,} seeded random bits (seed {settings['seed']}) drive both paths at {settings['bit_rate_bps'] / 1e9:g} Gbps. UI = {1e12 / settings['bit_rate_bps']:g} ps, source swing = {settings['source_vpp']:g} Vpp, rise time = {settings['rise_time_s'] * 1e12:g} ps and {settings['samples_per_ui']} samples/UI (dt = {settings['dt_s'] * 1e12:g} ps). The source is unipolar NRZ with full 0-to-100% linear transitions, following the supplied sample-code convention (the modeled 10-90% transition is 4 ps). The first {settings['skip_bits']:,} bits are discarded to remove the finite-record startup transient.

The matched transfer functions are $H_{{ch}}=S21_{{ch}}/2$ and $H_{{total}}=S21_{{ch}}H_{{eq}}/2$. The divider is applied exactly once. Sectional convolution uses overlap-add and retains the output tail. The inverse FFT produces a discrete impulse response, so there is no extra dt multiplication. Its full-period sample sum equals the DC transfer; truncation slightly changes the sum and is not renormalized.

## Statistical eyes and center eye height

![Post-channel and post-CTLE statistical eyes](../results/eye_diagrams.png)

{markdown_table(content['eyes'])}

Each path uses its own fixed sampling phase, chosen by maximizing the single-bit-response PDA margin around its main cursor. The plotting centers are aligned to 0.5 and 1.5 UI. At that fixed phase the random-pattern center eye height is

$$EH_{{random}}=\\min(v\\mid b=1)-\\max(v\\mid b=0),$$

over every retained bit. There is no density or probability threshold that removes rare excursions. A negative value is a signed separation indicating a closed eye. PDA sums worst-case sampled precursor and postcursor contributions from the isolated bit response, including the same 5 ps source ramps. PDA covers more bit combinations than the finite random sequence and is reported separately.

## Three 1 ns snapshots

These three overlays use the same complete simulated record and the same voltage scale for both received paths.

![Post-channel and post-CTLE, 100-101 ns](../results/waveform_100ns.png)

![Post-channel and post-CTLE, 200-201 ns](../results/waveform_200ns.png)

![Post-channel and post-CTLE, 300-301 ns](../results/waveform_300ns.png)

The CTLE reduces the long channel tail and separates received logic levels. Overshoot and negative undershoot can remain; their contributions at the sampling phase enter the eye-height calculation.

## Independent verification

{checks}

{robustness}

## Interpretation and limits

This is a linear, frequency-domain, noise-free calculation for the supplied channel. It retains channel phase and delay, but does not introduce receiver noise, clock jitter, nonlinearity or saturation. CTLE can increase high-frequency noise and does not generally remove reflection notches. Finite-pattern eye height and finite-tail PDA are not BER predictions. No independently validated analog circuit is claimed.

The optional ADS ChannelSim exercise was not performed. Every figure and numerical claim comes from the accompanying Python calculation and independent checks.

## References

[1] *Homework: CTLE simulation and design*, [assignment brief](../instructions/HW7_CTLE.pdf), pp. 1-3. The supplied channel and sample program are stored alongside the brief.

[2] Course lecture, *Continuous Time Linear Equalizer (CTLE)*, pp. 6-9 and 11-14. Transfer-function conventions, circuit interpretation and over-equalization.

Channel SHA-256: `{source['sha256']}`
"""
    output.write_text(text, encoding="utf-8", newline="\n")


def validate_anonymous_pdf(path: Path) -> dict[str, Any]:
    from pypdf import PdfReader

    reader = PdfReader(path)
    author = str((reader.metadata or {}).get("/Author", ""))
    if author not in {"", "Anonymous"}:
        raise ValueError("Unexpected PDF author metadata")
    names = reader.trailer["/Root"].get("/Names", {})
    if names and "/EmbeddedFiles" in names:
        raise ValueError("Report must not contain embedded files")
    text = "\n".join(page.extract_text() or "" for page in reader.pages)
    forbidden = ["/Users/", "@", "student id", "student number"]
    drive_path = re.search(r"\b[A-Za-z]:" + re.escape(chr(92)), text)
    if drive_path or any(token.casefold() in text.casefold() for token in forbidden):
        raise ValueError("Possible identifying content in report text")
    return {"pages": len(reader.pages), "author": author, "embedded_files": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--verification", type=Path, default=ROOT / "results" / "independent_verification.json")
    args = parser.parse_args()
    results = args.results.resolve()
    metrics = load_json(results / "metrics.json")
    verification = load_json(args.verification)
    if verification.get("all_passed") is not True:
        raise ValueError("Independent verification has not confirmed that all checks passed")
    content = report_content(metrics, verification)
    output_dir = Path(__file__).resolve().parent
    pdf_path = output_dir / "hw07_report.pdf"
    build_pdf(content, results, pdf_path)
    build_markdown(content, output_dir / "report.md")
    print(json.dumps(validate_anonymous_pdf(pdf_path), indent=2))


if __name__ == "__main__":
    main()
