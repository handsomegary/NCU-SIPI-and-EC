"""Build the anonymous HW14 report from actual simulation and verification outputs.

The report contains no contributor identity or machine-specific path.
Run after the analysis and independent verification; see the homework README.
"""

from __future__ import annotations

import argparse
import html
import hashlib
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
FIGURE_ROOT = RESULTS
PAGE_WIDTH, PAGE_HEIGHT = A4
MARGIN = 43
WIDTH = PAGE_WIDTH - 2 * MARGIN
INK = colors.HexColor("#193044")
ACCENT = colors.HexColor("#087B8A")
MUTED = colors.HexColor("#506473")
LINE = colors.HexColor("#D4E1E6")
PALE = colors.HexColor("#F0F6F7")


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
    if relative_path.startswith("../results/"):
        path = (FIGURE_ROOT / relative_path.removeprefix("../results/")).resolve()
    else:
        path = (HERE / relative_path).resolve()
    width, height = ImageReader(str(path)).getSize()
    limit = 330 if "eye_" in path.name else 240
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
    canvas.drawString(MARGIN, PAGE_HEIGHT - 25, "SIGNAL INTEGRITY / HW14")
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(MUTED)
    canvas.drawRightString(PAGE_WIDTH - MARGIN, PAGE_HEIGHT - 25, "NRZ AND PAM4")
    canvas.line(MARGIN, 33, PAGE_WIDTH - MARGIN, 33)
    canvas.drawString(MARGIN, 21, "Anonymous numerical report")
    canvas.drawRightString(PAGE_WIDTH - MARGIN, 21, str(document.page))
    canvas.restoreState()


def build(source, output):
    rl_config.invariant = True
    doc = BaseDocTemplate(
        str(output), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=51, bottomMargin=45,
        title="HW14 - NRZ versus PAM4", author="",
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


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def values(sequence, factor=1, digits=5):
    return ", ".join(f"{float(value) * factor:.{digits}g}" for value in sequence)


def report_markdown(summary, verification, test_results=None):
    """Compose the report from recorded evidence, never illustrative results."""
    if verification.get("successful") is not True:
        raise ValueError("Independent verification must pass before report generation")
    if test_results and test_results.get("successful") is not True:
        raise ValueError("The recorded unit-test suite did not pass")
    cases = summary["cases"]
    rates = sorted({float(case["bit_rate_gbps"]) for case in cases.values()})
    if len(rates) != 2 or len(cases) != 4:
        raise ValueError("This report requires exactly two rates and both NRZ and PAM4")
    ordered = []
    for rate in rates:
        for modulation in ("NRZ", "PAM4"):
            match = [(key, case) for key, case in cases.items()
                     if float(case["bit_rate_gbps"]) == rate
                     and case["modulation"].upper() == modulation]
            if len(match) != 1:
                raise ValueError("Missing or duplicate rate/modulation case")
            ordered.append(match[0])

    for case_key, case in ordered:
        eq = case["equalizers"]
        factor = 1 if case["modulation"].upper() == "NRZ" else 2
        if not math.isclose(case["symbol_rate_gbaud"] * factor, case["bit_rate_gbps"]):
            raise ValueError("Bit and symbol rates are inconsistent")
        if not math.isclose(case["ui_ps"] * case["symbol_rate_gbaud"], 1000):
            raise ValueError("The symbol UI is inconsistent with the symbol rate")
        if len(eq["ffe"]["taps"]) != 3 or len(eq["dfe"]["taps_v"]) != 5:
            raise ValueError("Update the report narrative for changed equalizer tap counts")
        if len(eq["ctle"]["poles_ghz"]) != 2 or eq["ctle"]["peak_gain_db"] > 18 + 1e-8:
            raise ValueError("The CTLE does not match the report's shared design constraints")
        if not math.isclose(case["design"]["design_clock_stress_ps"], 3.0):
            raise ValueError("Update the report narrative for changed phase-stress constraints")
        if not math.isclose(eq["ffe"]["l1_norm"], 1, abs_tol=1e-9):
            raise ValueError("The transmitter FFE swing constraint was not satisfied")
        expected_eyes = 1 if factor == 1 else 3
        if any(len(case["metrics"][key]) != expected_eyes for key in
               ("eye_heights_v", "eye_widths_ps", "eye_widths_ui")):
            raise ValueError("Incorrect number of measured eye openings")
        verified_case = verification['cases'][case_key]
        if (verified_case['measured_bits'] != case['measured_bits'] or
                verified_case['measured_symbols'] != case['measured_symbols']):
            raise ValueError("Verification refers to a different payload size")
        for metric_key in ('eye_heights_v', 'eye_widths_ps'):
            if not all(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)
                       for a, b in zip(verified_case[metric_key], case['metrics'][metric_key])):
                raise ValueError("Verification and reported simulation metrics do not match")

    params = summary["parameters"]
    for key, expected in {
        "source_vpp_v": 1.0, "rise_time_ps": 10.0, "sj_amplitude_ps": 1.5,
        "sj_frequency_mhz": 10.0, "rj_sigma_ps": 1.0,
        "ctle_maximum_absolute_gain_db": 18.0,
    }.items():
        if not math.isclose(params[key], expected, rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError(f"Update the report narrative for changed {key}")
    condition_rows = [
        ["Shared condition", "Value / convention"],
        ["Differential source", "1 Vpp open-circuit; levels centered on 0 V"],
        ["Reference impedances", "50 ohm at each single-ended source/load port"],
        ["Channel voltage transfer", "Sdd21 / 2 (matched Thevenin source convention)"],
        ["Transition duration", "10 ps, linear 0%-100% ramp for every level change"],
        ["Sinusoidal jitter", "1.5 ps peak; 10 MHz; evaluated in absolute time"],
        ["Random jitter", "Gaussian, standard deviation 1 ps"],
        ["Equalizer resources", "3 FFE taps; 5 DFE taps; CTLE with 1 zero and 2 poles"],
        ["CTLE gain limit", "18 dB maximum absolute frequency-response magnitude"],
        ["Transmitter FFE limit", "Sum of absolute tap coefficients is 1"],
        ["Design clock stress", "Nominal sampling phase and +/-3 ps; fixed DFE taps"],
        ["Payload / random seed", f"{params['bits']:,} common information bits / seed {params['seed']}"],
    ]
    comparison = [["Bit rate", "NRZ EH (mV)", "PAM4 EH (mV)", "NRZ EW (ps)", "PAM4 EW (ps)"]]
    findings = []
    for rate in rates:
        at_rate = {case["modulation"].upper(): case for _, case in ordered
                   if float(case["bit_rate_gbps"]) == rate}
        n, p = at_rate["NRZ"]["metrics"], at_rate["PAM4"]["metrics"]
        comparison.append([f"{rate:g} Gbps", f"{n['worst_eye_height_v']*1000:.3f}",
                           f"{p['worst_eye_height_v']*1000:.3f}",
                           f"{n['worst_eye_width_ps']:.3f}", f"{p['worst_eye_width_ps']:.3f}"])
        eh_winner = "NRZ" if n["worst_eye_height_v"] > p["worst_eye_height_v"] else "PAM4"
        ew_winner = "NRZ" if n["worst_eye_width_ps"] > p["worst_eye_width_ps"] else "PAM4"
        findings.append(f"At {rate:g} Gbps, {eh_winner} has the larger worst eye height, "
                        f"and {ew_winner} has the larger worst eye width in ps.")

    settings = [["Parameter"] + [f"{c['modulation'].upper()} {c['bit_rate_gbps']:g}G" for _, c in ordered]]
    settings.extend([
        ["Symbol rate (Gbaud)"] + [fnum(c["symbol_rate_gbaud"]) for _, c in ordered],
        ["Symbol UI (ps)"] + [fnum(c["ui_ps"]) for _, c in ordered],
        ["Time step (ps)"] + [fnum(c["dt_ps"]) for _, c in ordered],
        ["CTLE zero (GHz)"] + [fnum(c["equalizers"]["ctle"]["zero_ghz"]) for _, c in ordered],
        ["CTLE poles (GHz)"] + [values(c["equalizers"]["ctle"]["poles_ghz"]) for _, c in ordered],
        ["CTLE peak gain (dB)"] + [fnum(c["equalizers"]["ctle"]["peak_gain_db"]) for _, c in ordered],
        ["CTLE DC gain (dB)"] + [fnum(c["equalizers"]["ctle"]["dc_gain_db"]) for _, c in ordered],
        ["FFE taps"] + [values(c["equalizers"]["ffe"]["taps"]) for _, c in ordered],
        ["FFE L1 norm"] + [fnum(c["equalizers"]["ffe"]["l1_norm"]) for _, c in ordered],
        ["DFE taps (mV)"] + [values(c["equalizers"]["dfe"]["taps_v"], 1000) for _, c in ordered],
        ["DFE thresholds (mV)"] + [values(c["equalizers"]["dfe"]["thresholds_v"], 1000) for _, c in ordered],
        ["Main cursor (mV)"] + [fnum(c["equalizers"]["dfe"]["main_cursor_v"]*1000) for _, c in ordered],
    ])

    final_pages = []
    for key, case in ordered:
        m = case["metrics"]
        chunks = ["<!-- pagebreak -->", f"## Final {case['modulation'].upper()} eye at {case['bit_rate_gbps']:g} Gbps"]
        chunks.append(f"### {case['symbol_rate_gbaud']:g} Gbaud; symbol UI {case['ui_ps']:g} ps")
        chunks.append(f"![White dashed lines mark conditional extrema; dotted lines mark decision thresholds and the sampling clock. "
                      f"Voltages remain in physical receiver volts. Axis ranges differ between figures; compare labelled axes and measured values.](../results/eye_{key}.png)")
        chunks.append(markdown_table([
            ["Measured opening", "Eye height (mV)", "Eye width (ps)", "Eye width (UI)"],
            *[[str(i+1), f"{eh*1000:.3f}", f"{ew:.3f}", f"{eui:.4f}"]
              for i, (eh, ew, eui) in enumerate(zip(m["eye_heights_v"], m["eye_widths_ps"], m["eye_widths_ui"]))]
        ]))
        chunks.append(f"All {case['measured_symbols']:,} retained symbols contribute to the measured boundaries. "
                      f"The receiver recovered {case['measured_bits']:,} information bits with "
                      f"{m['bit_errors']} observed bit errors and {m['symbol_errors']} symbol errors. "
                      "The displayed eye is a histogram of the continuous held-feedback waveform; "
                      "its measured boundaries use the full retained record.")
        if case['modulation'].upper() == 'PAM4':
            chunks.append("Openings 1 through 3 run from the lowest to the highest threshold. "
                          "Every opening must remain usable; selecting only the best eye would overstate the link margin.")
        final_pages.append("\n\n".join(chunks))

    error_rows = [["Case", "Evaluated bits", "Evaluated symbols", "Bit errors", "Symbol errors"]]
    for _, case in ordered:
        m = case["metrics"]
        error_rows.append([f"{case['modulation'].upper()} / {case['bit_rate_gbps']:g}G",
                           f"{case['measured_bits']:,}", f"{case['measured_symbols']:,}",
                           str(m["bit_errors"]), str(m["symbol_errors"])])

    tests_sentence = "The saved-artifact independent verification completed successfully."
    if test_results:
        tests_count = test_results.get("tests_run", test_results.get("testsRun", "recorded"))
        tests_sentence += f" The separate unit-test suite passed {tests_count} tests."
    elif verification.get("unit_test_count"):
        tests_sentence += f" The independent unit-test suite passed {verification['unit_test_count']} tests."
    verification_section = verification_markdown(summary, verification)
    width_qualifications = []
    for rate in rates:
        keys = {case['modulation'].upper(): key for key, case in ordered
                if float(case['bit_rate_gbps']) == rate}
        nkey, pkey = keys['NRZ'], keys['PAM4']
        baseline_delta = cases[nkey]['metrics']['worst_eye_width_ps'] - cases[pkey]['metrics']['worst_eye_width_ps']
        studies = summary.get('validation', {})
        if nkey in studies and pkey in studies:
            ntrial = studies[nkey]['independent_seed']['metrics']['worst_eye_width_ps']
            ptrial = studies[pkey]['independent_seed']['metrics']['worst_eye_width_ps']
            if baseline_delta * (ntrial - ptrial) < 0:
                width_qualifications.append(
                    f"At {rate:g} Gbps the eye-width ordering changes with the independent seed "
                    f"(NRZ {ntrial:.3f} ps; PAM4 {ptrial:.3f} ps). The small baseline width "
                    "advantage is therefore seed-sensitive; it is not a robust format advantage.")
    caveat = ' '.join(width_qualifications)
    width_unit_notes = []
    for item in summary['comparison']:
        if item['larger_minimum_eye_width_ps'] != item['larger_minimum_eye_width_ui']:
            width_unit_notes.append(
                f"At {item['bit_rate_gbps']:g} Gbps, "
                f"{item['larger_minimum_eye_width_ps'].upper()} has the larger width in ps, "
                f"while {item['larger_minimum_eye_width_ui'].upper()} has the larger width in symbol UI. "
                "The distinction arises because PAM4's symbol interval is twice as long.")
    units_note = ' '.join(width_unit_notes)
    loss_rows = [["Bit rate", "NRZ Nyquist loss", "PAM4 Nyquist loss", "Difference"]]
    for item in summary['comparison']:
        loss_rows.append([f"{item['bit_rate_gbps']:g} Gbps",
                          f"{item['nrz_nyquist_loss_db']:.3f} dB @ {item['nrz_nyquist_ghz']:g} GHz",
                          f"{item['pam4_nyquist_loss_db']:.3f} dB @ {item['pam4_nyquist_ghz']:g} GHz",
                          f"{item['nyquist_loss_difference_db']:.3f} dB"])
    first_design = ordered[0][1]['design']

    return f"""> NUMERICAL EXPERIMENT / ANONYMOUS REPORT

# NRZ versus PAM4

## Same information rate, different channel demands

This experiment compares two-level NRZ and four-level PAM4 over the supplied
50 cm differential channel. The route contains 40 cm of differential stripline,
10 cm of differential microstrip, and distributed parasitic elements.
The four required final eye diagrams use FFE, the measured channel, CTLE, and
decision-feedback equalization at two information rates.

{markdown_table(comparison)}

{' '.join(findings)} These statements describe the finite recorded experiment
under the settings below. Each PAM4 entry is the minimum of its three openings.
Eye width is compared in absolute ps; normalized UI values are also reported.

{caveat}

{markdown_table(condition_rows)}

PAM4 carries two bits per symbol, so at a fixed bit rate its symbol interval is
twice the NRZ interval. With the same source Vpp, the adjacent PAM4 voltage
spacing is one third of the NRZ spacing. The experiment tests whether lower
channel loss at PAM4's reduced symbol rate compensates for that vertical penalty.

<!-- pagebreak -->

## Channel and voltage conventions

Input ports are 1 and 3; output ports are 2 and 4. S-parameters use the convention
Sij = wave leaving port i / wave incident on port j. The forward differential
term is `Sdd21 = (S21 - S23 - S41 + S43) / 2`.

The waveform generator describes the open-circuit differential source voltage.
Matched 50 ohm source and load ports therefore give the voltage transfer
`Hchannel = Sdd21 / 2`. An ideal matched through connection consequently produces
half the open-circuit source voltage at the load. Keeping this factor explicit
prevents an unnoticed factor-of-two error in eye height.

![Differential insertion loss with the NRZ and PAM4 Nyquist frequencies for each selected bit rate.](../results/channel.png)

{markdown_table(loss_rows)}

The source file spans 0 through 80 GHz on a 10 MHz grid. Its 100 ns inverse-transform
period limits the available unaliased delay record. A real Hermitian IDFT uses
the measured response directly and sets frequencies above 80 GHz to zero.
The causal kernel retains {params['impulse_ns']:g} ns, without DC renormalization.
Baseline waveforms use {params['samples_per_ui']} samples per symbol UI;
the numerical sensitivity to both choices is checked separately. No hardware
measurement or new circuit-simulator result is claimed.

The retained loaded-channel kernel sums to {first_design['channel_kernel_dc']:.9f},
compared with the full-record DC value {first_design['channel_expected_dc']:.9f}.
This small truncation difference is retained rather than hidden by rescaling.

For rate Rb, NRZ uses Rs = Rb and PAM4 uses Rs = Rb/2. The symbol Nyquist
frequency is Rs/2; this is distinct from the first rectangular-pulse spectral
null at Rs. The lecture's approximate 9.54 dB differential-loss crossover is a
screening heuristic, not a substitute for the final eye comparison.

<!-- pagebreak -->

## Equalization and clock model

Equalizer coefficients may differ among the four cases, but all cases use the
same tap counts, CTLE pole/zero count, transmitter swing constraint and CTLE
absolute gain cap. Twenty CTLE candidates are combined with a conservative
worst-pattern linear program at the nominal sampling phase and +/-3 ps offsets.
DFE coefficients and midpoint thresholds stay fixed across these design offsets.
This is a finite robustness objective; it is not a bound on Gaussian jitter.
The design uses the noiseless channel response before the independent payload
is evaluated. The result is a reproducible finite design search, not a claim
of globally optimal transceiver performance.

{markdown_table(settings)}

The transmitter first applies its symbol-spaced FFE to the data levels. All
resulting level transitions are then placed on the jittered transmitter clock.
The same absolute SJ amplitude and RJ standard deviation are used for both
formats. Using a common percentage of UI instead would give PAM4 twice the
absolute timing error and would violate the assignment.

The differential PAM4 Gray levels, from low to high, represent 00, 01, 11 and 10.
DFE feeds back previously recovered decisions and uses three voltage thresholds
for PAM4 or one for NRZ. It does not replace recovered decisions with the known
transmitted payload. The receiver clock and thresholds are fixed from the design.
The feedback value for symbol n is held from half a UI before its clock through
half a UI after it; the right endpoint belongs to the next hold. This remains
causal because every feedback term uses a decision from at least one symbol earlier.

{chr(10).join(final_pages)}

<!-- pagebreak -->

## Measurement definitions and finite-record evidence

Eye height is the finite-record separation between adjacent conditional voltage
envelopes at the receiver's design sampling phase. The upper extreme of the lower
symbol class is subtracted from the lower extreme of the next class. Negative
values indicate overlapping voltage distributions. No percentile trimming or
3-sigma opening is substituted for these extrema.

Eye width is the contiguous time interval around the sampling phase over which
the corresponding fixed decision threshold lies between those class envelopes.
Both ps and symbol UI are retained. PAM4's three eyes are checked separately;
the minimum opening is used for the format comparison.

{markdown_table(error_rows)}

The error counts come from recovered decisions matched to the transmitted
payload, including Gray decoding for PAM4. A zero count means no error was seen
in that record. It does not establish zero BER or demonstrate a 1e-12 BER target.
Gaussian timing jitter is unbounded, and finite-record extrema depend on record
length and seed. No receiver voltage noise, clock recovery, nonlinear device
effects or FEC is added to this model. The CTLE is an ideal linear gain stage;
receiver voltage rails and saturation are not modeled.

<!-- pagebreak -->

## Independent verification and sensitivity

{tests_sentence}

{verification_section}

<!-- pagebreak -->

## Interpretation and reproduction

{' '.join(findings)} The physical tradeoff is between reduced bandwidth demand
and reduced separation of adjacent voltage levels. Channel loss and equalizer
constraints determine which effect dominates at a particular bit rate.
The measured results, rather than the 9.54 dB heuristic alone, support the
conclusion for this specific channel and model.

{caveat}

{units_note}

CTLE magnitude is bounded by the same absolute cap across cases. The source FFE
normalization keeps the emitted symbol range within the shared Vpp budget.
PAM4 is assessed by its weakest of three openings. These controls prevent an
unfair comparison through larger transmitted swing, extra equalizer complexity,
selective choice of the best PAM4 eye, or per-figure amplitude normalization.

### Reproduce from the repository root

```sh
python -m pip install -r homework/hw14/requirements.txt
python homework/hw14/verify_inputs.py
python homework/hw14/src/hw14.py --validation
python -m unittest discover -s homework/hw14/tests -v
python homework/hw14/tests/verify_results.py
python -m pip install -r homework/hw14/report/requirements.txt
python homework/hw14/report/build_report.py
```

The simulator resolves its input paths relative to its own homework directory
and runs without an interactive plotting window. The result archive retains
the settings and numerical arrays required by the independent verifier.
This report is regenerated from the recorded JSON and figures rather than
manually entered measurement tables.

### Sources and provenance

[1] Supplied HW14 assignment, `instructions/HW14_NRZ_vs_PAM4.pdf`, pages 1-4.
[2] Supplied `instructions/HW14_diff_channel_50cm.s4p`.
[3] Course handout 15, PAM4, pages 2-5 and 8-11, for symbol rates, Gray encoding,
three eye openings, spectral comparison, the empirical loss rule and PAM4 DFE.
[4] N. Dikhaminjia et al., "High-speed serial link challenges using multi-level
signaling," IEEE EPEPS, 2015, as cited by the assignment.

The supplied assignment is preserved. Personally identifying export comments
were removed from the channel copy; its numerical data were preserved exactly
and checksums are recorded in `instructions/input_manifest.json`. This report,
its figures and its metadata do not identify a student or contributor.
"""


def verification_markdown(summary, verification):
    """Render compact independent-check and convergence evidence."""
    rows = [["Verification record", "Outcome"]]
    for key, value in verification.items():
        if key in ("successful", "passed") or isinstance(value, (dict, list)):
            continue
        if isinstance(value, (int, float, bool, str)):
            rows.append([key.replace("_", " "), fnum(value) if isinstance(value, float) else str(value)])
    body = markdown_table(rows) if len(rows) > 1 else "The complete independent check record is saved in `results/verification.json`."
    verified_cases = list(verification['cases'].values())
    tx_count = sum(case['end_to_end_source_samples'] for case in verified_cases)
    decision_count = sum(case['end_to_end_decisions_checked'] for case in verified_cases)
    cells = sum(case['end_to_end_histogram_cells_checked'] for case in verified_cases)
    decision_error = max(case['end_to_end_max_decision_voltage_error_v'] for case in verified_cases)
    envelope_error = max(case['end_to_end_max_conditional_envelope_error_v'] for case in verified_cases)
    body += (f"\n\nThe independent baseline reconstruction generated {tx_count:,} transmitter "
             f"samples from analytic ramps, used a full complex IDFT and overlap-add convolution, "
             f"and checked all {decision_count:,} receiver decision inputs, including guards. "
             f"Maximum decision-voltage disagreement was {decision_error:.3g} V. Conditional "
             f"boundaries at every saved phase agreed within {envelope_error:.3g} V, and all "
             f"{cells:,} histogram cells matched exactly. Separate checks cover odd-mode "
             "port excitation, source scaling, FFE shifts, recovered-decision feedback, "
             "Gray decoding and eye measurements.")
    validation = summary.get("validation")
    if validation:
        checks = [["Case / perturbation", "EH change (mV)", "EW change (ps)", "Bit errors"]]
        for key, study in validation.items():
            base = summary["cases"][key]
            for variant, trial in study.items():
                measured = trial["metrics"]
                if variant == "independent_seed":
                    label = f"Seed {trial['seed']}"
                elif variant == "time_refinement":
                    label = f"{trial['samples_per_ui']} samples/UI"
                elif variant == "impulse_extension":
                    label = f"{trial['impulse_ns']:g} ns impulse"
                else:
                    label = variant.replace("_", " ")
                checks.append([
                    f"{base['modulation'].upper()} {base['bit_rate_gbps']:g}G / {label}",
                    f"{(measured['worst_eye_height_v']-base['metrics']['worst_eye_height_v'])*1000:+.4f}",
                    f"{measured['worst_eye_width_ps']-base['metrics']['worst_eye_width_ps']:+.4f}",
                    str(measured["bit_errors"]),
                ])
        body += "\n\n" + markdown_table(checks)
        body += ("\n\nAll sensitivity runs freeze the equalizer coefficients and design clock; "
                 "the time-grid study changes the numerical grid, and the impulse study "
                 "changes the retained channel duration. The independent-seed study changes "
                 "the random payload and jitter realization. Changes are reported relative "
                 "to each baseline case. Seed sensitivity describes finite-record variability, "
                 "not a deterministic discretization error. The sensitivity runs retain "
                 "metrics rather than complete raw waveforms; the separate verifier checks "
                 "their consistency, while the baseline has the saved numerical archives.")
    else:
        raise ValueError("The published report requires the recorded --validation study")
    return body


def main():
    global FIGURE_ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument("--output", type=Path, default=HERE / "report.pdf")
    parser.add_argument("--refresh", action="store_true", help="accepted for consistency; results are always refreshed")
    args = parser.parse_args()
    FIGURE_ROOT = args.results.resolve()
    summary = load_json(args.results / "metrics.json")
    verification = load_json(args.results / "verification.json")
    for name, expected in verification.get('verified_file_sha256', {}).items():
        path = HERE.parent / name if name.startswith('src/') else args.results / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Verified file changed: {name}; rerun independent verification")
    tests_path = args.results / "test_results.json"
    tests = load_json(tests_path) if tests_path.exists() else None
    source = HERE / "report.md"
    source.write_text(report_markdown(summary, verification, tests), encoding="utf-8")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    build(source, args.output)


if __name__ == "__main__":
    main()
