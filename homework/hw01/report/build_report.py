"""Build an anonymous, reproducible HW01 PDF from the verified result files."""
from __future__ import annotations

import json
import os
from io import BytesIO
from pathlib import Path
from xml.sax.saxutils import escape

from pypdf import PdfReader, PdfWriter

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak,
    KeepTogether, Preformatted,
)

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
SUMMARY = RESULTS / "summary.json"
PDF = HERE / "report.pdf"
MARKDOWN = HERE / "report.md"

PAGE_W, PAGE_H = 595.276, 841.89
LEFT, RIGHT, TOP, BOTTOM = 52, 52, 57, 48
WIDTH = PAGE_W - LEFT - RIGHT
INK = colors.HexColor("#182C3B")
BLUE = colors.HexColor("#1D5F82")
MUTED = colors.HexColor("#516575")
RULE = colors.HexColor("#D4DEE5")
LIGHT = colors.HexColor("#F1F5F7")

# Embed a widely available Unicode font without binding the report to a host path.
try:
    windows_fonts = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    if (windows_fonts / "arial.ttf").is_file():
        font_files = [windows_fonts / name for name in ("arial.ttf", "arialbd.ttf", "consola.ttf")]
    else:
        import matplotlib
        font_dir = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
        font_files = [font_dir / name for name in ("DejaVuSans.ttf", "DejaVuSans-Bold.ttf", "DejaVuSansMono.ttf")]
    for family, path in zip(("ReportSans", "ReportSans-Bold", "ReportMono"), font_files):
        pdfmetrics.registerFont(TTFont(family, path))
    pdfmetrics.registerFontFamily("ReportSans", normal="ReportSans", bold="ReportSans-Bold")
    BODY, BOLD, MONO = "ReportSans", "ReportSans-Bold", "ReportMono"
except (ImportError, OSError):
    BODY, BOLD, MONO = "Helvetica", "Helvetica-Bold", "Courier"

STYLES = {
    "eyebrow": ParagraphStyle("eyebrow", fontName=BOLD, fontSize=9, leading=13, textColor=BLUE, spaceAfter=8),
    "title": ParagraphStyle("title", fontName=BOLD, fontSize=25, leading=30, textColor=INK, spaceAfter=10),
    "subtitle": ParagraphStyle("subtitle", fontName=BODY, fontSize=11, leading=16, textColor=MUTED, spaceAfter=18),
    "h1": ParagraphStyle("h1", fontName=BOLD, fontSize=17, leading=22, textColor=INK, spaceAfter=12),
    "h2": ParagraphStyle("h2", fontName=BOLD, fontSize=10.5, leading=15, textColor=BLUE, spaceBefore=10, spaceAfter=5),
    "body": ParagraphStyle("body", fontName=BODY, fontSize=9.25, leading=14.2, textColor=INK, spaceAfter=8),
    "small": ParagraphStyle("small", fontName=BODY, fontSize=8.2, leading=12, textColor=MUTED, spaceAfter=7),
    "caption": ParagraphStyle("caption", fontName=BODY, fontSize=8.1, leading=11.5, textColor=MUTED, spaceBefore=5, spaceAfter=11),
    "table": ParagraphStyle("table", fontName=BODY, fontSize=8.5, leading=12, textColor=INK),
    "tablehead": ParagraphStyle("tablehead", fontName=BOLD, fontSize=8.3, leading=11.5, textColor=INK),
    "equation": ParagraphStyle("equation", fontName=BODY, fontSize=10, leading=16, textColor=INK, spaceBefore=4, spaceAfter=9),
    "code": ParagraphStyle("code", fontName=MONO, fontSize=7.5, leading=10.4, textColor=INK, spaceBefore=5, spaceAfter=9),
}


def p(text, style="body"):
    return Paragraph(text, STYLES[style])


def fig(filename, max_height, caption):
    path = RESULTS / filename
    if not path.is_file():
        raise FileNotFoundError(path.name)
    width, height = ImageReader(str(path)).getSize()
    scale = min(WIDTH / width, max_height / height)
    image = Image(str(path), width=width * scale, height=height * scale)
    image.hAlign = "CENTER"
    return KeepTogether([image, p(caption, "caption")])


def table(headers, rows, widths):
    cells = [[p(str(v), "tablehead") for v in headers]]
    cells += [[p(str(v), "table") for v in row] for row in rows]
    result = Table(cells, colWidths=widths, hAlign="LEFT", repeatRows=1)
    result.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), LIGHT),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, RULE),
        ("LINEBELOW", (0, 1), (-1, -1), 0.35, RULE),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
        ("LEFTPADDING", (0, 0), (-1, -1), 9),
        ("RIGHTPADDING", (0, 0), (-1, -1), 9),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    return result


def page_header(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(RULE)
    canvas.setLineWidth(0.5)
    canvas.line(LEFT, PAGE_H - 35, PAGE_W - RIGHT, PAGE_H - 35)
    canvas.setFont(BOLD, 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(LEFT, PAGE_H - 26, "EE8091  /  HOMEWORK 01")
    canvas.setFont(BODY, 8)
    canvas.drawRightString(PAGE_W - RIGHT, PAGE_H - 26, "KAISER WINDOW SPECTRAL ANALYSIS")
    canvas.line(LEFT, 33, PAGE_W - RIGHT, 33)
    canvas.setFont(BODY, 7.5)
    canvas.drawString(LEFT, 21, "Reproducible numerical report")
    canvas.drawRightString(PAGE_W - RIGHT, 21, f"{doc.page:02d}")
    canvas.restoreState()


def number(value, digits=6):
    if abs(value) < 0.5 * 10**(-digits):
        value = 0.0
    return f"{value:.{digits}f}"


def build():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    sample = summary["sampling"]
    fit = summary["fit"]
    components = fit["components"]
    n = sample["n_samples"]
    fft_size = sample["fft_size"]
    rmse = fit["rmse_v"]
    max_residual = fit["max_abs_residual_v"]
    order_rows = fit.get("order_selection", [])
    fit_method = summary["analysis"]["fit_method"]

    story = []
    story += [p("SPECTRAL ESTIMATION", "eyebrow"), p("Spectrum analysis<br/>with a Kaiser window", "title")]
    story += [p("Homework 01 | Direct window implementation, amplitude calibration, and signal reconstruction", "subtitle")]
    story += [p("Objective", "h2"), p("Estimate the number of sinusoidal components Q, their amplitudes A<sub>i</sub>, and their frequencies f<sub>i</sub> from the supplied two-column timing waveform. The assignment requires a time-domain plot, spectra for at least four distinct Kaiser parameters, gain compensation, and Python source code [1].")]
    story += [fig("time_domain.png", 212, "Figure 1. The supplied time-domain waveform. The varying envelope indicates multiple frequencies; the waveform alone does not determine their number or amplitudes.")]
    story += [table(["Record property", "Value", "Interpretation"], [
        ["Samples", f"{n:,}", "Uniformly spaced voltage observations"],
        ["Sample interval / rate", f"{sample['dt_s']*1e12:.0f} ps / {sample['fs_hz']/1e9:.0f} GHz", "Nyquist frequency: 50 GHz"],
        ["Record length / endpoint span", f"{sample['record_duration_s']*1e9:.2f} ns / {sample['sample_span_s']*1e9:.2f} ns", "N Δt and (N - 1) Δt have different roles"],
        ["Native DFT spacing", f"{sample['unpadded_df_hz']/1e6:.0f} MHz", "1 / (N Δt)"],
    ], [140, 157, WIDTH-297])]
    story += [Spacer(1, 10), p(f"Result in brief: the combined window analysis and reconstruction support Q = {fit['q']}. Two approximately 1 V components occur near 525 and 602 MHz, and a much weaker component occurs near 888.5 MHz. The final time-domain fit has RMSE {rmse:.3e} V.", "body")]

    story += [PageBreak(), p("1. Windowing and amplitude calibration", "h1")]
    story += [p("Direct implementation", "h2"), p("For sample indices n = 0, ..., N - 1, the implemented window follows the assignment's modified-Bessel equation directly. No built-in Kaiser-window constructor is used.")]
    story += [p("w[n] = I<sub>0</sub>(π α √(1 - [2n/(N - 1) - 1]<super>2</super>)) / I<sub>0</sub>(π α)", "equation")]
    story += [p("This course uses α. Libraries that parameterize a Kaiser window by β use β = π α. Substituting α directly for β would create a different taper. The α = 0 case reduces to the rectangular window; increasing α reduces the endpoints and broadens the main lobe.")]
    story += [fig("kaiser_windows.png", 190, "Figure 2. Directly evaluated Kaiser windows. The plotted parameter labels use the assignment's α convention.")]
    story += [p("Window gain and a one-sided spectrum", "h2"), p("Let Y[k] be the M-point DFT of x[n]w[n], with zero padding when M &gt; N. The coherent gain is G<sub>c</sub> = (1/N) Σw[n]. For an isolated, bin-centered sinusoid away from DC and Nyquist, a calibrated single-sided amplitude is")]
    story += [p("A[k] = 2 |Y[k]| / Σw[n] = 2 |Y[k]| / (N G<sub>c</sub>).", "equation")]
    story += [p("DC is not doubled. For even M, the Nyquist bin is not doubled either; an odd-length real FFT has no Nyquist bin. The denominator uses the N real samples and their window sum, never the padded length M. Zero padding contributes no additional signal samples.")]
    story += [p("The plotted vertical scale is 20 log<sub>10</sub>(A / 1 V), labeled dBV. Thus 1 V corresponds to 0 dBV and 0.00316228 V corresponds to -50 dBV. The displayed floor is only a plotting limit; the fit uses the original voltage samples.")]
    story += [p("Calibration limits", "h2"), p("Coherent-gain correction removes the taper's amplitude loss. It does not remove off-grid scalloping or bias from overlapping lobes. These effects motivate refinement with a time-domain sinusoidal model instead of reading every amplitude directly from a peak height.")]

    story += [PageBreak(), p("2. The spectrum changes with α", "h1")]
    story += [p("The four individual spectrum figures on this page and the next use the same time record, frequency units, amplitude reference, and gain correction. Each figure labels its own α value, as required by the assignment.")]
    story += [fig("spectrum_alpha_0p1.png", 194, "Figure 3. α = 0.1. A narrow main lobe helps distinguish the nearby strong pair, while sidelobes mask the distant weak component.")]
    story += [fig("spectrum_alpha_1.png", 194, "Figure 4. α = 1. Sidelobe levels decrease, with a wider main lobe than for α = 0.1.")]
    story += [p("Reading the comparison", "h2"), p("Small α retains a narrow main lobe, which helps distinguish the two strong components separated by roughly 77 MHz. Large α lowers sidelobes so that leakage from the strong pair is less likely to mask a distant weak tone. The main lobe simultaneously widens, so no single window setting is best for both tasks.")]
    story += [p("Peak heights can change even after gain compensation when lobes overlap or a peak is sampled away from its maximum. The window comparison is therefore used to identify candidate structure; the final table comes from a joint fit.")]

    story += [PageBreak(), p("3. Recovering the weak component", "h1")]
    story += [fig("spectrum_alpha_4.png", 227, "Figure 5. α = 4. Strong-tone leakage is reduced, making the weak feature near 0.8885 GHz visible on the calibrated amplitude scale.")]
    story += [fig("spectrum_alpha_5.png", 227, "Figure 6. α = 5. The same weak feature persists under a different taper. This cross-check supports treating it as signal structure rather than relying on a single-window sidelobe.")]
    story += [p("Why Q = 2 is incomplete", "h2"), p("A small-α spectrum can emphasize the two large peaks and obscure the third component. A model containing only the strong pair leaves a coherent residual; adding the weak component removes that residual and improves the numerical fit. Window comparison and residual verification are both needed to establish Q.")]
    story += [p("Zero padding and resolving power", "h2"), p(f"The FFT uses M = {fft_size:,} points, giving a displayed grid spacing of {sample['plotted_df_hz']/1e6:.6g} MHz. The acquired record still lasts {sample['record_duration_s']*1e9:.0f} ns and its native spacing remains {sample['unpadded_df_hz']/1e6:.0f} MHz. Padding interpolates the sampled spectrum; it does not narrow a window's main lobe or add independent information. The joint fit uses a parametric model and the full time record.")]

    story += [PageBreak(), p("4. Estimated signal and reconstruction", "h1")]
    story += [p("Fitting procedure", "h2"), p(escape(str(fit_method)))]
    story += [p("At each trial set of frequencies, the sine coefficients, cosine coefficients, and constant offset are solved by linear least squares. Frequency refinement is nonlinear. The fitted amplitude is √(s<sub>i</sub><super>2</super> + c<sub>i</sub><super>2</super>), and the phase is atan2(c<sub>i</sub>, s<sub>i</sub>). This allows numerical checks of the zero-phase sinusoidal model in the prompt without forcing the phase in advance.")]
    story += [table(["Component", "Frequency (MHz)", "Amplitude (V)", "Amplitude (dBV)"], [
        [str(row["index"]), number(row["frequency_mhz"], 6), number(row["amplitude_v"], 10), number(row["amplitude_dbv"], 4)]
        for row in components
    ], [65, 137, 148, WIDTH-350])]
    story += [Spacer(1, 10), p(f"Selected model: Q = {fit['q']}. Fitted offset: {fit['offset_v']:.3e} V. Reconstruction RMSE: {rmse:.3e} V. Maximum absolute residual: {max_residual:.3e} V. Relative RMSE: {fit['relative_rmse']:.3e}.")]
    story += [fig("reconstruction.png", 265, "Figure 7. Reconstruction from all fitted components, together with the residual. The residual scale is reported explicitly so that an apparently perfect overlay is not mistaken for an independent accuracy check.")]
    story += [p("Precision and interpretation", "h2"), p("The table retains enough digits to reproduce the numerical result. These are parameter estimates for a supplied deterministic waveform; a small residual does not establish corresponding measurement accuracy or uncertainty for a physical experiment. The fitted phases and model-order diagnostics are retained in the machine-readable results.")]

    story += [PageBreak(), p("5. Reproducibility and verification", "h1")]
    story += [p("Numerical model-order evidence", "h2")]
    if order_rows:
        story += [table(["Q", "RMSE (V)", "Residual sum of squares", "BIC"], [
            [str(row["q"]), f"{row['rmse_v']:.3e}", f"{row['rss_v2']:.3e}", f"{row['bic']:.3f}"]
            for row in order_rows
        ], [40, 110, 190, WIDTH-340])]
        story += [Spacer(1, 9), p(f"The RSS column is in V<super>2</super>. For BIC, the code uses unit-RMS-normalized samples and floors normalized RSS at N × ({summary['analysis']['numerical_relative_tolerance']:.0e})<super>2</super>. It evaluates N ln(RSS<sub>norm,floored</sub>/N) + (3Q + 1) ln N. The normalization changes the common scale, not the model-order comparison. The search stops when the relative numerical residual tolerance is reached; the physical-voltage RMSE remains unfloored.")]
    story += [p("Core amplitude-calibration algorithm", "h2")]
    snippet = """n = np.arange(N, dtype=float)
z = 2*n/(N - 1) - 1
w = np.i0(np.pi*alpha*np.sqrt(np.maximum(0, 1 - z*z)))
w /= np.i0(np.pi*alpha)
Y = np.fft.rfft(x*w, n=M)
A = np.abs(Y) / np.sum(w)
if M % 2 == 0:
    A[1:-1] *= 2      # DC and Nyquist stay undoubled
else:
    A[1:] *= 2        # no Nyquist bin for odd M
f = np.fft.rfftfreq(M, d=dt)"""
    story += [Preformatted(snippet, STYLES["code"])]
    story += [p("Repository evidence", "h2"), p("The accompanying Python source performs input validation, computes every plot, estimates the model parameters, and writes summary.json and parameters.csv. Figures in this report are embedded directly from those generated result files. The repository README gives the command for reproducing the analysis and the checks; the complete source is supplied with the report, satisfying the code requirement.")]
    story += [table(["Check", "Evidence"], [
        ["Input integrity", "N = 2000; uniformly increasing time; finite samples"],
        ["Window definition", "Direct modified-Bessel evaluation; β = π α"],
        ["Amplitude handling", "Coherent gain; correct DC/Nyquist treatment; volts and dBV"],
        ["Signal completeness", "Multiple α values plus joint-fit residual and model-order results"],
    ], [130, WIDTH-130])]
    story += [p("References", "h2"), p("[1] Course-provided assignment, <i>Time-frequency transform (1) - Using Kaiser window for spectral analysis</i>. Source waveform: HW1_Kaiser.txt. The assignment and waveform are reference inputs; no personal identifiers are reproduced.", "small")]

    doc = SimpleDocTemplate(str(PDF), pagesize=(PAGE_W, PAGE_H), leftMargin=LEFT, rightMargin=RIGHT,
        topMargin=TOP, bottomMargin=BOTTOM, title="HW01 - Spectrum analysis with a Kaiser window",
        author="", subject="Anonymous reproducible spectral analysis report", creator="SIPI and EC report builder")
    doc.build(story, onFirstPage=page_header, onLaterPages=page_header)
    # Replace, rather than append to, metadata so author names and timestamps
    # cannot be inherited from an environment or an older document.
    reader = PdfReader(BytesIO(PDF.read_bytes()))
    writer = PdfWriter()
    writer.clone_document_from_reader(reader)
    writer.metadata = {
        "/Title": "HW01 - Spectrum analysis with a Kaiser window",
        "/Subject": "Anonymous reproducible spectral analysis report",
        "/Author": "",
        "/Creator": "SIPI and EC report builder",
    }
    with PDF.open("wb") as stream:
        writer.write(stream)
    write_markdown(summary)
    print(json.dumps({"pdf": PDF.name, "markdown": MARKDOWN.name, "rmse_v": rmse, "q": fit["q"]}, indent=2))


def write_markdown(summary):
    fit = summary["fit"]
    sample = summary["sampling"]
    rows = "\n".join(f"| {r['index']} | {number(r['frequency_mhz'], 6)} | {number(r['amplitude_v'], 10)} | {number(r['amplitude_dbv'], 4)} |" for r in fit["components"])
    text = f"""# HW01 - Spectrum analysis with a Kaiser window

EE8091. Anonymous reproducible numerical report.

## Problem and record

Estimate the number, amplitudes, and frequencies of the sinusoids in the supplied waveform. The record contains {sample['n_samples']} samples, with a {sample['dt_s']*1e12:.0f} ps sample interval, {sample['fs_hz']/1e9:.0f} GHz sampling rate, {sample['record_duration_s']*1e9:.0f} ns record duration, and {sample['sample_span_s']*1e9:.2f} ns endpoint span. Its native DFT spacing is {sample['unpadded_df_hz']/1e6:.0f} MHz.

![Time-domain waveform](../results/time_domain.png)

## Kaiser window and calibration

Evaluate the assignment equation directly using the modified Bessel function I0:

```text
w[n] = I0(pi*alpha*sqrt(1 - (2*n/(N-1) - 1)^2)) / I0(pi*alpha)
```

The course convention is beta = pi*alpha. No built-in Kaiser constructor is used. For the M-point real FFT Y of the windowed samples, divide abs(Y) by sum(w), then double the positive-frequency bins except Nyquist. DC is never doubled. An odd M has no Nyquist bin. The denominator refers to the original N samples rather than the zero-padded FFT length. Report 20*log10(A/1 V) in dBV.

![Direct Kaiser windows](../results/kaiser_windows.png)

## Window comparison

![Spectra for different Kaiser parameters](../results/spectra.png)

![Individual spectrum with alpha 0.1](../results/spectrum_alpha_0p1.png)

![Individual spectrum with alpha 1](../results/spectrum_alpha_1.png)

Small alpha helps separate the nearby strong tones through a narrow main lobe. Large alpha suppresses leakage and exposes the distant weak tone while broadening the main lobe. Coherent-gain compensation does not correct scalloping or overlapping lobes.

![Spectrum with alpha 4](../results/spectrum_alpha_4.png)

![Spectrum with alpha 5](../results/spectrum_alpha_5.png)

The weak feature persists for multiple tapers. A two-tone model leaves a structured residual. Combining window evidence with residual fitting supports Q = {fit['q']}.

The FFT length is {sample['fft_size']}; its grid spacing is {sample['plotted_df_hz']/1e6:.6g} MHz. Padding interpolates the spectrum and does not improve the acquired record's resolving power.

## Joint fit and reconstruction

{summary['analysis']['fit_method']}

For fixed frequencies, solve sine coefficients, cosine coefficients, and the offset by linear least squares; refine the frequencies nonlinearly. Amplitude is sqrt(s^2+c^2); phase is atan2(c,s). Retain model-order evidence in summary.json.

| Component | Frequency (MHz) | Amplitude (V) | Amplitude (dBV) |
| --- | ---: | ---: | ---: |
{rows}

Selected Q: {fit['q']}. Offset: {fit['offset_v']:.3e} V. RMSE: {fit['rmse_v']:.3e} V. Maximum absolute residual: {fit['max_abs_residual_v']:.3e} V. Relative RMSE: {fit['relative_rmse']:.3e}.

Model-order selection compares BIC using unit-RMS-normalized samples, with normalized RSS floored at N*(1e-10)^2. The penalty counts 3*Q+1 parameters (frequency, sine coefficient, cosine coefficient, and one overall offset). Physical-voltage RMSE is unfloored. The algorithm stops on its declared relative numerical residual tolerance.

![Reconstruction and residual](../results/reconstruction.png)

The displayed digits describe a deterministic numerical fit and do not establish experimental uncertainty. The complete Python source, input validation, plots, numerical outputs, and reproduction command accompany this report in the repository.

## Rebuild the report

First generate the analysis outputs using the command in the homework README. From the homework/hw01 directory:

```sh
python -m pip install -r report/requirements.txt
python report/build_report.py
```

The builder reads results/summary.json and the generated PNG figures, then rewrites report/report.pdf and report/report.md. It clears inherited author/date metadata and uses no personal information. ReportLab and pypdf are the report-building dependencies; Matplotlib is an optional font fallback on non-Windows systems.

## Reference

Course-provided assignment: *Time-frequency transform (1) - Using Kaiser window for spectral analysis*. Source waveform: HW1_Kaiser.txt. No personal identifiers are reproduced.
"""
    MARKDOWN.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    build()
