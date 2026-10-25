"""Build the anonymous HW10 report from regenerated simulation results.

The report contains no copied original report pages or personal metadata.
Run from any directory; the default results directory is ../results.
"""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

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
PAGE_W, PAGE_H = A4
MARGIN = 44
WIDTH = PAGE_W - 2 * MARGIN
INK = colors.HexColor("#173349")
TEAL = colors.HexColor("#007C83")
MUTED = colors.HexColor("#536472")
LINE = colors.HexColor("#D5E0E7")
PALE = colors.HexColor("#F0F5F8")

STYLES = {
    "title": ParagraphStyle("Title", fontName="Helvetica-Bold", fontSize=28,
                            leading=32, textColor=INK, spaceAfter=15),
    "heading": ParagraphStyle("Heading", fontName="Helvetica-Bold", fontSize=19,
                              leading=24, textColor=INK, spaceAfter=13),
    "subheading": ParagraphStyle("Subheading", fontName="Helvetica-Bold", fontSize=11,
                                 leading=15, textColor=TEAL, spaceBefore=8,
                                 spaceAfter=7, keepWithNext=True),
    "body": ParagraphStyle("Body", fontName="Helvetica", fontSize=10,
                           leading=14, textColor=INK, spaceAfter=8),
    "small": ParagraphStyle("Small", fontName="Helvetica", fontSize=8.5,
                            leading=11.8, textColor=MUTED, spaceAfter=7),
    "caption": ParagraphStyle("Caption", fontName="Helvetica", fontSize=8.5,
                              leading=11.5, textColor=MUTED, spaceBefore=3,
                              spaceAfter=12),
    "cell": ParagraphStyle("Cell", fontName="Helvetica", fontSize=8.5,
                           leading=11.5, textColor=INK),
    "headcell": ParagraphStyle("HeadCell", fontName="Helvetica-Bold", fontSize=8.5,
                               leading=11.5, textColor=colors.white),
    "code": ParagraphStyle("Code", fontName="Courier", fontSize=8.2,
                           leading=11.8, textColor=INK, backColor=PALE,
                           borderPadding=8, spaceBefore=4, spaceAfter=12),
    "kicker": ParagraphStyle("Kicker", fontName="Helvetica-Bold", fontSize=9,
                             leading=13, textColor=TEAL, spaceAfter=7),
}


def p(text: str, style: str = "body") -> Paragraph:
    return Paragraph(text, STYLES[style])


def plain(value: Any) -> str:
    return html.escape(str(value))


def num(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def table(rows: list[list[Any]], widths: list[float] | None = None) -> Table:
    if widths is None:
        widths = [WIDTH / len(rows[0])] * len(rows[0])
    data = [[p(plain(cell), "headcell" if i == 0 else "cell")
             for cell in row] for i, row in enumerate(rows)]
    item = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    item.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), TEAL),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LINEBELOW", (0, -1), (-1, -1), .5, LINE),
    ]))
    return item


def figure(path: Path, caption: str, height: float = 270) -> KeepTogether:
    width_px, height_px = ImageReader(str(path)).getSize()
    scale = min(WIDTH / width_px, height / height_px)
    image = Image(str(path), width=width_px * scale, height=height_px * scale)
    image.hAlign = "CENTER"
    return KeepTogether([image, p(caption, "caption")])


def code(text: str) -> Preformatted:
    return Preformatted(text, STYLES["code"])


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(LINE)
    canvas.line(MARGIN, 37, PAGE_W - MARGIN, 37)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(MARGIN, 25, "HW10 | Differential channel equalization")
    canvas.drawRightString(PAGE_W - MARGIN, 25, str(doc.page))
    canvas.restoreState()


def metadata(canvas, doc):
    canvas.setTitle("HW10 - Differential channel simulation and equalizer design")
    canvas.setAuthor("")
    canvas.setSubject("Reproducible differential channel simulation and equalizer verification")
    canvas.setCreator("HW10 report builder")
    canvas.setKeywords("differential signaling, CTLE, FFE, DFE, signal integrity")
    footer(canvas, doc)


def build_story(metrics: dict[str, Any], results: Path, audit_override: Path | None = None) -> list:
    """Create content using the solver's saved metrics; never invent results."""
    parameters = metrics["parameters"]
    equalizers = metrics["equalizers"]
    ctle = equalizers["ctle"]
    ffe = equalizers["ffe"]
    dfe = equalizers["dfe"]
    stages = metrics["stages"]
    final = stages["dfe"]
    validation = metrics.get("validation", {})
    passed = final["eye_height_v"] > .120 and final["eye_width_ps"] > 15
    stage_names = {"channel": "Channel", "ctle": "Channel + CTLE",
                   "ffe": "FFE + channel + CTLE", "dfe": "FFE + channel + CTLE + DFE"}
    story: list = []
    add = story.append
    add(p("HW10  /  SIGNAL INTEGRITY", "kicker"))
    add(p("Differential channel<br/>simulation and<br/>equalizer design", "title"))
    add(p("Rewritten implementation and numerical verification", "subheading"))
    add(p("A four-port differential channel is simulated with a transmitter feed-forward "
          "equalizer (FFE), a continuous-time linear equalizer (CTLE), and a decision-feedback "
          "equalizer (DFE). The design sequence follows the assignment: CTLE, then FFE, then DFE. "
          "The physical signal path places the FFE before the channel and the other two equalizers at the receiver."))
    add(Spacer(1, 7))
    status = "MEETS THE SPECIFIED FINITE-RECORD TARGET" if passed else "DOES NOT MEET THE SPECIFIED TARGET"
    add(p(status, "subheading"))
    add(p(f"The regenerated post-DFE eye has <b>{num(final['eye_height_v']*1000, 3)} mV</b> eye height "
          f"and <b>{num(final['eye_width_ps'], 3)} ps</b> eye width. The assignment requires strict "
          "inequalities: EH &gt; 120 mV and EW &gt; 15 ps. These measurements apply to the simulated "
          "finite bit record and the model and clock convention stated in this report."))
    add(table([
        ["Operating condition", "Value"],
        ["Bit rate / unit interval", f"{num(parameters['bit_rate_bps']/1e9, 1)} Gb/s / {num(parameters['ui_ps'], 1)} ps"],
        ["Driver rise time", f"{num(parameters['rise_ps'], 1)} ps"],
        ["Each single-ended ideal source swing", f"{num(parameters['driver_vpp_v'], 1)} V peak-to-peak"],
        ["Source and load per conductor", "50 ohm / 50 ohm"],
        ["Simulation record", f"{parameters['n_bits']:,} bits; seed {parameters['seed']}"],
        ["Time step / samples per UI", f"{num(parameters['dt_ps'], 3)} ps / {parameters['samples_per_ui']}"],
    ], [WIDTH*.43, WIDTH*.57]))
    add(Spacer(1, 14))
    add(table([["Signal path", "Signed EH (mV)", "EW (ps)"]]+[
        [stage_names[key], num(stages[key]["eye_height_v"]*1000, 3), num(stages[key]["eye_width_ps"], 3)]
        for key in stage_names], [WIDTH*.60, WIDTH*.20, WIDTH*.20]))
    add(p("Negative signed EH means the conditional one/zero populations overlap at the selected "
          "clock; the physical eye is closed. The signed value is retained to quantify that overlap.", "small"))

    add(PageBreak())
    add(p("Channel and voltage normalization", "heading"))
    add(p("The supplied Touchstone file is a four-port network with a 50-ohm reference impedance. "
          "Ports (1, 3) form the input pair and ports (2, 4) form the output pair. The positive "
          "conductors are ports 1 and 2. Mixed-mode normalization uses a 100-ohm differential "
          "reference and a 25-ohm common-mode reference."))
    add(code("Sdd21(f) = [S21(f) + S43(f) - S23(f) - S41(f)] / 2\nHchannel(f) = Sdd21(f) / 2"))
    add(p("The second factor of one half is the source-to-port impedance division. It is included "
          "once in the channel transfer function, so the pulse response, random-bit waveform, "
          "and DFE coefficients share the same voltage scale. Sdd21 itself is a wave transmission "
          "coefficient; the S-parameter comparison below does not include this source divider."))
    add(figure(results / "sparameters.png", "Figure 1. Single-ended S21 and differential Sdd21 of the supplied channel.", 275))
    add(p("Causal transmitter pulse", "subheading"))
    add(p("Data symbols are b[n] in {-1, +1}, with pure differential excitation. The differential source has levels +/-Vpp before "
          "loading, where Vpp is the peak-to-peak swing of each single-ended source. A finite-rise "
          "single-bit pulse is formed by smoothing a rectangular UI pulse with a causal rise-time rectangle."))
    add(code("P(f) = Vpp Tb sinc(f Tb) sinc(f Tr) exp[-j pi f (Tb + Tr)]\nSBR(t) = IFT{ P(f) Hchannel(f) W(f) }\ny(t) = sum_n b[n] SBR(t - n Tb)"))
    add(p(f"Here sinc(x) = sin(pi x)/(pi x). The frequency-domain Kaiser parameter is "
          f"{num(parameters['kaiser_alpha'], 2)}. The time-domain response duration is "
          f"{num(parameters['impulse_duration_ns'], 2)} ns. Numerically, the causal trapezoid is sampled "
          "and convolved with the IDFT channel kernel. The expression above states the corresponding "
          "continuous-time pulse convention. Frequency sampling, windowing, truncation, and time "
          "discretization are part of the numerical model.", "small"))

    add(PageBreak())
    add(p("CTLE design and gain limit", "heading"))
    add(p("The CTLE compensates the high-frequency loss of the differential channel. Its "
          "transfer function is expressed in a DC-normalized pole-zero form, with pole and zero "
          "frequencies converted to angular frequencies before evaluation."))
    add(code("Hctle(s) = Gdc product_z(1 + s/wz) / product_p(1 + s/wp)\nGdc = 10^(DC_gain_dB/20),   wp = 2 pi fp,   wz = 2 pi fz"))
    add(table([
        ["CTLE quantity", "Selected value", "Assignment limit"],
        ["Poles (GHz)", ", ".join(num(v, 5) for v in ctle["poles_ghz"]), "At most 3; highest <= 40 GHz"],
        ["Zeros (GHz)", ", ".join(num(v, 5) for v in ctle["zeros_ghz"]), "No separate zero constraint"],
        ["DC gain", f"{num(ctle['dc_gain_db'], 5)} dB", "Included in full gain check"],
        ["Numerator gain parameter", f"{num(ctle['numerator_gain_db'], 5)} dB", "Not the realized peak gain"],
        ["Maximum gain", f"{num(ctle['peak_gain_db'], 6)} dB", "<= 6 dB"],
        ["Peak-gain frequency", f"{num(ctle['peak_frequency_ghz'], 5)} GHz", "Numerically located peak"],
    ], [WIDTH*.25, WIDTH*.37, WIDTH*.38]))
    add(Spacer(1, 13))
    add(figure(results / "frequency_response.png", "Figure 2. CTLE and combined transfer functions; gain is evaluated over frequency.", 310))
    add(p("A numerator coefficient described as an 'AC gain' in a parameterization does not "
          "necessarily equal the maximum magnitude of the realized transfer function. The constraint "
          "is checked against the actual peak of |Hctle|, rather than the name or numerical value of "
          "a coefficient. The selected equalizer is evaluated together with the channel and transmitter response."))

    add(PageBreak())
    add(p("FFE and DFE coefficients", "heading"))
    add(p("The FFE acts on the differential symbol sequence. Taps are listed in causal order: "
          "tap k multiplies b[n-k]. Its absolute coefficient sum fixes the transmitter swing. "
          "The DFE subtracts the listed voltage coefficient times the corresponding previous decision; "
          "negative coefficients are permitted."))
    validation_seeds = [parameters["seed"]] + [v["seed"] for v in validation.get("seed_checks", [])]
    add(p("The coefficient search generated 16,384 symbols with seed 70 and retained 14,336 after "
          "1,024-symbol guards at each end. The final settings are checked on separate "
          f"{parameters['n_bits']:,}-bit records with seeds {plain(', '.join(map(str, validation_seeds)))}; "
          "the validation records are not the coefficient-fitting record.", "small"))
    add(code("x[n] = sum_k c[k] b[n-k],       sum_k |c[k]| = 1\nz[n] = y[n] - sum_i a[i] d[n-i],       d[n] in {-1, +1}"))
    rows = [["FFE tap k", "Delay (UI)", "Coefficient c[k]"]]
    rows += [[str(i), str(i), num(v, 9)] for i, v in enumerate(ffe["taps"])]
    add(table(rows, [WIDTH*.27, WIDTH*.27, WIDTH*.46]))
    add(p(f"FFE tap count: {len(ffe['taps'])}; absolute sum: {num(ffe['l1_norm'], 10)}. "
          "The largest-magnitude coefficient occurs at index "
          f"{max(range(len(ffe['taps'])), key=lambda i: abs(ffe['taps'][i]))}.", "small"))
    rows = [["DFE feedback i", "Previous decision", "Subtracted coefficient (V)"]]
    rows += [[str(i+1), f"d[n-{i+1}]", num(v, 9)] for i, v in enumerate(dfe["taps_v"])]
    add(table(rows, [WIDTH*.27, WIDTH*.27, WIDTH*.46]))
    add(p(f"DFE threshold: {num(dfe['threshold_v'], 3)} V. Decision phase sample: "
          f"{dfe['decision_phase_sample']} of {parameters['samples_per_ui']} samples/UI. "
          "Each feedback value is held from half a UI before its decision instant up to, but not "
          "including, half a UI afterward. The right endpoint uses the next symbol's feedback value. "
          "The update convention is identical in waveform synthesis and eye extraction.", "small"))
    add(figure(results / "sbr.png", "Figure 3. Single-bit responses used to select the linear equalizers and DFE feedback coefficients.", 215))

    add(PageBreak())
    add(p("Eye diagrams: channel and CTLE", "heading"))
    add(p("The same deterministic bit sequence is used across stages. Eye diagrams show the "
          "voltage axis in the correct increasing direction. Numerical eye measurements use the "
          "underlying trace extrema, rather than reading the image or its histogram bins.", "small"))
    for key, number in [("channel", 4), ("ctle", 5)]:
        stage = stages[key]
        add(figure(results / f"eye_{key}.png", f"Figure {number}. {stage_names[key]}. "
                   f"EH = {num(stage['eye_height_v']*1000, 3)} mV; EW = {num(stage['eye_width_ps'], 3)} ps.", 275))

    add(PageBreak())
    add(p("Eye diagrams: FFE and DFE", "heading"))
    add(p("DFE uses its own preceding decisions in the reported simulation. The eye therefore "
          "includes any decision-feedback errors that occur in the selected finite record; it is "
          "not generated by replacing past decisions with the known transmitted symbols.", "small"))
    for key, number in [("ffe", 6), ("dfe", 7)]:
        stage = stages[key]
        add(figure(results / f"eye_{key}.png", f"Figure {number}. {stage_names[key]}. "
                   f"EH = {num(stage['eye_height_v']*1000, 3)} mV; EW = {num(stage['eye_width_ps'], 3)} ps.", 275))

    add(PageBreak())
    add(p("Output waveforms and eye measurement", "heading"))
    add(figure(results / "waveforms.png", "Figure 8. Output waveforms for all four requested signal paths, with their timing convention shown on the plot.", 300))
    add(p("At each phase within one UI, traces are grouped by the transmitted symbol. Eye height "
          "is the minimum logic-one voltage minus the maximum logic-zero voltage at the selected "
          "decision clock. No percentile clipping is applied. Eye width is the contiguous interval "
          "around that clock where every logic-one sample is above 0 V and every logic-zero sample "
          "is below 0 V. Its two edges are interpolated using the more restrictive threshold margin."))
    add(p(f"The analysis retains {parameters['n_bits']:,} symbols after a "
          f"{parameters.get('warmup_symbols', 1024):,}-symbol warmup and adds "
          f"{parameters.get('trailing_guard_symbols', 1024):,} trailing guard symbols. "
          "The eye plots are histograms; the reported extrema are computed from every retained trace.", "small"))
    add(table([
        ["Constraint", "Measured or selected value", "Limit"],
        ["Post-DFE eye height", f"{num(final['eye_height_v']*1000, 3)} mV", "> 120 mV"],
        ["Post-DFE eye width", f"{num(final['eye_width_ps'], 3)} ps", "> 15 ps"],
        ["CTLE pole count", str(len(ctle["poles_ghz"])), "<= 3"],
        ["Highest CTLE pole", f"{num(max(ctle['poles_ghz']), 5)} GHz", "<= 40 GHz"],
        ["Peak CTLE gain", f"{num(ctle['peak_gain_db'], 6)} dB", "<= 6 dB"],
        ["FFE / DFE tap counts", f"{len(ffe['taps'])} / {len(dfe['taps_v'])}", "<= 5 / <= 5"],
        ["FFE absolute coefficient sum", num(ffe["l1_norm"], 10), "1"],
    ], [WIDTH*.36, WIDTH*.40, WIDTH*.24]))
    add(p("The final eye result is a simulation measurement. It is not a statistical BER guarantee "
          "and does not certify every possible bit pattern, noise process, or physical implementation.", "small"))

    add(PageBreak())
    add(p("Verification, corrections, and limits", "heading"))
    add(p("The legacy waveform applied the source divider after simulation; the rewrite carries "
          "it once in the transfer function. It also corrects the dimensionally invalid Tb Tr pulse "
          "phase to the causal Tb + Tr convention, uses a common Kaiser setting, aligns the DFE "
          "clock, reuses one bit record across stages, and preserves the eye's voltage orientation."))
    add(p("Reproducibility checks", "subheading"))
    add(p(f"The recorded post-DFE decision-error count is {final.get('decision_errors', 'not recorded')}. "
          "Random seeds, time-step studies, and numerical checks below are retained with the results "
          "so that the measurement can be repeated independently."))
    checks = [{"seed": parameters["seed"], "n_bits": parameters["n_bits"], **final}] + validation.get("seed_checks", [])
    if checks:
        rows = [["Seed", "Bits", "EH (mV)", "EW (ps)", "Errors"]]
        for item in checks:
            rows.append([item.get("seed", "-"), item.get("n_bits", "-"),
                         num(item["eye_height_v"]*1000, 3), num(item["eye_width_ps"], 3),
                         item.get("decision_errors", "-")])
        add(table(rows, [WIDTH*.12, WIDTH*.19, WIDTH*.23, WIDTH*.23, WIDTH*.23]))
        add(Spacer(1, 8))
    refinement = validation.get("refinement")
    if refinement:
        if isinstance(refinement, dict):
            refinement = [refinement]
        rows = [["Time step (ps)", "Bits", "EH (mV)", "EW (ps)"]]
        for item in refinement:
            rows.append([num(item["dt_ps"], 3), item.get("n_bits", "-"),
                         num(item["eye_height_v"]*1000, 3), num(item["eye_width_ps"], 3)])
        add(table(rows, [WIDTH*.25]*4))
        add(Spacer(1, 8))
    audit_path = audit_override or results / "independent_checks.json"
    if audit_path.exists():
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        audit_checks = audit.get("checks", [])
        if audit_checks:
            successes = sum(bool(item.get("pass")) for item in audit_checks)
            add(p(f"Independent verification: <b>{successes}/{len(audit_checks)} checks pass</b>. "
                  "The separate reconstruction applies FFE to the source levels, then the ramp and "
                  "time-domain channel convolution, followed by serial DFE and an independently "
                  "implemented eye measurement. Detailed comparisons are saved in independent_checks.json.", "small"))
            comparisons = audit.get("production_comparison", {})
            wave_errors = audit.get("waveform_max_abs_errors_V", {})
            if comparisons and wave_errors:
                eh_error = max(c["EH_abs_error_V"] for c in comparisons.values())
                ew_error = max(c["EW_abs_error_ps"] for c in comparisons.values())
                add(p(f"Across all four stages, the largest EH difference is {eh_error:.2e} V "
                      f"and the largest EW difference is {ew_error:.2e} ps. The maximum waveform "
                      f"difference over 2,000 checked samples per stage is {max(wave_errors.values()):.2e} V.", "small"))
            failures = [item["name"] for item in audit_checks if not item.get("pass")]
            if failures:
                add(p("Failed independent checks: " + plain(", ".join(failures)), "small"))
    add(p("Scope of validation", "subheading"))
    add(p("The model has finite frequency samples and response support, an ideal fixed-rate clock, "
          "and idealized equalizers. Finite-record extrema do not guarantee every bit pattern. "
          "Noise, jitter, nonlinearity, and adaptation are not modeled. The optional ADS ChannelSim "
          "comparison has not been run.", "small"))
    add(p("Sources and reproduction", "subheading"))
    add(p("[1] HW10_diff_signaling.pdf, supplied assignment, pp. 1 and 4.<br/>"
          "[2] Differential Signaling (1), course handout, pp. 14-17 and 20-31.<br/>"
          "[3] Supplied four-port channel and reference implementation. All plotted results "
          "are regenerated by the rewritten analysis.", "small"))
    add(code("python src/hw10.py --validation\npython -m unittest discover -s tests -v\npython tests/verify_reference.py\npython report/build_report.py"))
    return story


def build(results: Path, output: Path, audit: Path | None = None) -> None:
    metrics_path = results / "metrics.json"
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    story = build_story(metrics, results, audit)
    output.parent.mkdir(parents=True, exist_ok=True)
    rl_config.invariant = True
    document = BaseDocTemplate(str(output), pagesize=A4, leftMargin=MARGIN,
                               rightMargin=MARGIN, topMargin=41, bottomMargin=49,
                               title="HW10 - Differential channel simulation and equalizer design",
                               author="", pageCompression=1)
    frame = Frame(MARGIN, 49, WIDTH, PAGE_H - 90, leftPadding=0,
                  rightPadding=0, topPadding=0, bottomPadding=0)
    document.addPageTemplates(PageTemplate(id="Report", frames=[frame], onPage=metadata))
    document.build(story)
    from pypdf import PdfReader
    reader = PdfReader(output)
    if (reader.metadata or {}).get("/Author", ""):
        raise ValueError("The report Author metadata must remain empty")
    print(f"Created {output.name}: {len(reader.pages)} pages")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=HERE.parent / "results")
    parser.add_argument("--output", type=Path, default=HERE / "HW10_report.pdf")
    parser.add_argument("--audit", type=Path, help="Optional independent-check JSON; defaults to results/independent_checks.json")
    args = parser.parse_args()
    build(args.results.resolve(), args.output.resolve(), args.audit.resolve() if args.audit else None)


if __name__ == "__main__":
    main()
