"""Build the anonymous HW4 report from the package's reproducible results.

Run ``python report/build_report.py`` after ``python src/hw4.py``.
Requires numpy, matplotlib, and reportlab. No external services are used.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak,
    KeepTogether,
)


HERE = Path(__file__).resolve().parent
BASE = HERE.parent
RESULTS = BASE / "results"
ASSETS = HERE / "assets"
PDF = HERE / "hw04_report.pdf"
MARKDOWN = HERE / "report.md"
NAVY = colors.HexColor("#15334f")
BLUE = colors.HexColor("#1558a6")
INK = colors.HexColor("#203040")
MUTED = colors.HexColor("#526579")
PALE = colors.HexColor("#eef4f9")
WIDTH, HEIGHT = A4
CONTENT = WIDTH - 104


def fmt(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}g}"


def load_inputs():
    summary = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
    audit = json.loads((RESULTS / "pattern_audit.json").read_text(encoding="utf-8"))
    spectra = np.load(RESULTS / "spectra.npz")
    return summary, audit, spectra


def make_assets(summary, spectra):
    ASSETS.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "savefig.dpi": 180})
    expressions = {
        "estimator": r"$\widehat S_1(f_j)=q_j\,\frac{\Delta t}{Z_0\sum_{n=0}^{N-1}w_n^2}\,\frac{1}{K}\sum_{k=1}^{K}|\mathrm{FFT}_M\{w_n v_{k,n}\}_j|^2$",
        "reference": r"$S_{\mathrm{NRZ},1}(f)=\frac{2A^2T_b}{Z_0}\,\mathrm{sinc}^2(fT_b)\,\mathrm{sinc}^2(fT_r),\quad f>0$",
        "coefficient": r"$c_m=\frac{A}{10}\,\mathrm{sinc}(f_mT_b)\,\mathrm{sinc}(f_mT_r)\,e^{-j\phi_m}\left[e^{-j2\pi m9/10}-e^{-j2\pi m8/10}\right]$",
        "hybrid": r"$S_{\mathrm{hybrid},1}(f)=0.8S_{\mathrm{NRZ},1}(f)+\sum_{m=1}^{\infty}\frac{2|c_m|^2}{Z_0}\,\delta(f-f_m),\quad f_m=m\,\mathrm{GHz}$",
        "linepower": r"$P_m=\frac{8A^2}{100Z_0}\,\mathrm{sinc}^2(m/10)\,\mathrm{sinc}^2(m/50)\,\sin^2(\pi m/10)$",
    }
    for name, expression in expressions.items():
        fig = plt.figure(figsize=(12, 0.8))
        fig.text(0.01, 0.52, expression, fontsize=17, va="center", color="#15334f")
        fig.savefig(ASSETS / f"{name}.png", bbox_inches="tight", pad_inches=0.09, facecolor="white")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
    for p, color in zip(summary["patterns"], ("#1558a6", "#c7473d", "#27816e")):
        axes[0].plot(np.arange(1, 11), p["statistics"]["ones_fraction_by_position_modulo_10"],
                     "o-", ms=4, label=f"Pattern {p['pattern']}", color=color)
        bits = np.array([int(x) for x in (BASE / "instructions" / p["input_file"]).read_text().strip() if x in "01"])
        running = np.r_[0, np.cumsum(2 * bits.astype(np.int64) - 1)]
        selected = np.unique(np.r_[np.arange(0, len(running), 40), len(running) - 1])
        axes[1].plot(selected / 1000, running[selected], lw=1, color=color, label=f"Pattern {p['pattern']}")
    axes[0].set(xlabel="Position in aligned 10-bit word (1-based)", ylabel="Observed fraction of ones", ylim=(-0.05, 1.05), xticks=np.arange(1, 11))
    axes[0].axhline(0.5, ls=":", lw=1, color="gray")
    axes[1].set(xlabel="Bits processed (thousands)", ylabel="Cumulative signed bit sum")
    for ax in axes:
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ASSETS / "pattern_diagnostics.png")
    plt.close(fig)

    m = np.arange(1, 50)
    a = summary["configuration"]["vpp"] / 2
    z = summary["configuration"]["z0"]
    power = 8 * a*a / (100*z) * np.sinc(m/10)**2 * np.sinc(m/50)**2 * np.sin(np.pi*m/10)**2
    fig, ax = plt.subplots(figsize=(10, 3.2))
    markerline, stemlines, baseline = ax.stem(m, power * 1e6, basefmt=" ")
    plt.setp(markerline, markersize=3.4, markerfacecolor="#1558a6", markeredgecolor="#1558a6")
    plt.setp(stemlines, color="#1558a6", linewidth=1.2)
    ax.set(xlabel="Positive harmonic frequency (GHz)", ylabel="Integrated line power (uW)", xlim=(0, 50), ylim=(0, float(np.max(power) * 1e6 * 1.13)))
    ax.grid(alpha=0.22)
    fig.tight_layout()
    fig.savefig(ASSETS / "hybrid_predicted_line_powers.png")
    plt.close(fig)

    f = spectra["frequency_hz"]
    selected = (f >= 0.1e9) & (f <= 9e9)
    theory = spectra["sampled_random_reference_w_per_hz"]
    agreement = {}
    for k in (10, 100):
        error = 10 * np.log10(spectra[f"pattern_1_K{k}_w_per_hz"][selected] / theory[selected])
        agreement[str(k)] = {
            "mean_db": float(np.mean(error)),
            "rms_db": float(np.sqrt(np.mean(error**2))),
            "mean_linear_ratio": float(np.mean(10**(error/10))),
        }
    return agreement


def make_text(summary, audit, agreement):
    cfg = summary["configuration"]
    der = summary["derived_parameters"]
    p = summary["patterns"]
    worst_parseval = max(e["parseval_relative_error"] for x in p for e in x["estimates"].values())
    raw_random = p[0]["estimates"]["100"]["mean_time_domain_power_w"] * 1000
    weighted_random = p[0]["estimates"]["100"]["mean_window_weighted_power_w"] * 1000
    theory_random = der["independent_random_sampled_expected_total_power_w"] * 1000
    with np.load(RESULTS / "spectra.npz") as spec:
        lowmask = (spec["frequency_hz"] >= 10e6) & (spec["frequency_hz"] <= 100e6)
        lowratio = float(10*np.log10(np.mean(spec["pattern_2_K100_w_per_hz"][lowmask]) / np.mean(spec["pattern_1_K100_w_per_hz"][lowmask])))
    verification = json.loads((RESULTS / "verification.json").read_text(encoding="utf-8"))
    return {
        "cfg": cfg, "der": der, "patterns": p, "worst_parseval": worst_parseval,
        "raw_random": raw_random, "weighted_random": weighted_random,
        "theory_random": theory_random, "lowratio": lowratio,
        "agreement": agreement, "verification": verification,
    }


def write_markdown(d):
    cfg, der, p = d["cfg"], d["der"], d["patterns"]
    rows = "\n".join(
        f"| {x['pattern']} | {k} | {x['estimates'][str(k)]['mean_window_weighted_power_w']*1000:.6f} | "
        f"{x['estimates'][str(k)]['integrated_full_band_psd_w']*1000:.6f} | {x['estimates'][str(k)]['parseval_relative_error']:.2e} |"
        for x in p for k in (10, 100)
    )
    md = f"""# HW4: Spectrum of digital signals

Anonymous technical report

## Scope and conclusions

Three supplied 200,000-bit records are converted to bipolar trapezoidal NRZ and analyzed with a correctly normalized, one-sided averaged periodogram. Six required PSD plots cover 0-50 GHz in dBm/Hz at K = 10 and K = 100. Pattern 1 is consistent with independent random NRZ. Pattern 2 is the assignment's intended 8b/10b-like data, with strong low-frequency suppression; its supplied bitstream nevertheless fails basic strict 8b/10b invariants. Pattern 3 has eight random-looking payload bits followed by an exact 01 suffix in every aligned 10-bit word, producing a continuous background plus a 1 GHz harmonic comb.

## Parameters and signal model

| Quantity | Value |
|---|---|
| Input length per file | 200,000 bits |
| Bit rate / UI | 10 Gbit/s / 100 ps |
| Full record duration | 20 us |
| Voltage / impedance | -0.5 V and +0.5 V / 50 ohm |
| Linear transition duration | 20 ps |
| Sample interval / Nyquist frequency | 2 ps / 250 GHz |
| Samples per bit / per transition | 50 / 10 |
| Segment length | 2,000 bits = 100,000 samples = 200 ns |
| Unwindowed resolution scale | 1/T = 5 MHz |
| Kaiser window | beta = pi, equivalently alpha = 1 when beta = pi alpha |
| Window equivalent noise bandwidth | {der['window_enbw_hz']/1e6:.6f} MHz |
| FFT length | 524,288 = 2^19 samples |
| Display frequency grid | {der['fft_grid_step_hz']/1e6:.9f} MHz |
| Averages | 10 and 100 consecutive, nonoverlapping segments |

The full record contains 200,000 bits, not 2,000,000. At 10 Gbit/s it lasts 20 us. K = 10 uses the first 20,000 bits (2 us); K = 100 uses all 200,000 bits (20 us). The two estimates share their first 10 segments and are not independent experiments. Segment duration and window remain fixed so that changes in K primarily alter averaging variance, not spectral resolution.

Let A = Vpp/2 = 0.5 V. Bit 0 maps to -A and bit 1 to +A. At a transition the preceding level is retained at the boundary, followed by a linear ramp reaching the new level after 20 ps. Successive segments use the actual preceding bit rather than restarting the waveform. Only the start of the complete record assumes preceding bit 0. This is causal rectangular smoothing of ideal NRZ; its delay changes spectral phase but not the PSD. No mean subtraction or detrending is applied.

## Estimator and normalization

For N samples per segment, M FFT points and q_j = 2 at interior positive bins, q_j = 1 at DC and at the even-length FFT's Nyquist bin,

$$\\widehat S_1(f_j)=q_j\\frac{{\\Delta t}}{{Z_0\\sum_n w_n^2}}\\frac1K\\sum_{{k=1}}^K|\\mathrm{{FFT}}_M\\{{w_nv_{{k,n}}\\}}_j|^2.$$

This is equivalent to dividing the squared continuous-time transform approximation |dt FFT|^2 by T Z0 U, where U = mean(w^2). Average squared magnitudes, never complex FFT coefficients. The stored spectra are W/Hz. Plot values are 10 log10(1000 S1) dBm/Hz. DC and Nyquist must not be doubled.

The window power gain is U = {der['window_power_gain_U']:.8f}, coherent gain Gc = mean(w) = {der['window_coherent_gain']:.8f}, and ENBW = U/Gc^2 times 1/T = {der['window_enbw_hz']/1e6:.6f} MHz. U corrects broadband PSD power. Coherent gain instead corrects the amplitude of a bin-centered deterministic tone; substituting Gc^2 for U in a PSD changes the noise-floor normalization by the ENBW factor. Zero padding provides a smoother frequency grid but does not improve the 200 ns observation's physical resolving power. Both 1/T and ENBW are below 10 MHz; the plotted grid is below 1 MHz.

For independent, equiprobable bits, the assignment's one-sided continuous-time reference is

$$S_{{NRZ,1}}(f)=\\frac{{2A^2T_b}}{{Z_0}}\\operatorname{{sinc}}^2(fT_b)\\operatorname{{sinc}}^2(fT_r),\\quad f>0,$$

with normalized sinc(x) = sin(pi x)/(pi x). It is a reference for random NRZ, not a prediction that the other two inputs should follow. Pulse-shape nulls occur at multiples of 10 GHz and the first transition-filter null is 50 GHz. At the single DC grid point the code uses the undoubled two-sided reference value, consistent with the rFFT endpoint convention. The finite-sample waveform's exact reference includes sinc(f dt)^(-4) relative to the continuous formula away from removable zeros; this is only {der['sampled_vs_continuous_reference_at_50_GHz_db']:.4f} dB by 50 GHz in its limiting ratio. Both models are stored for verification.

## Verification and uncertainty

The complete positive sampled band, 0-250 GHz, satisfies the DFT Parseval check sum(S1) df = mean(sum(w^2 v^2)/sum(w^2))/Z0. Its largest relative error across all six estimates is {d['worst_parseval']:.2e}. Full-band integration must use the DFT bin sum, not trapezoidal endpoint half weights.

| Pattern | K | Weighted time-domain power (mW) | Full-band PSD power (mW) | Relative error |
|---|---:|---:|---:|---:|
{rows}

For pattern 1 at K = 100, ordinary time-domain power is {d['raw_random']:.6f} mW and window-weighted power is {d['weighted_random']:.6f} mW. The independent-bit sampled-model expectation is {d['theory_random']:.6f} mW; the ideal instantaneous-transition value of 5 mW is reduced by finite ramps. In 0.1-9 GHz, the RMS pointwise log error against the sampled random reference decreases from {d['agreement']['10']['rms_db']:.4f} dB at K = 10 to {d['agreement']['100']['rms_db']:.4f} dB at K = 100. The corresponding mean linear ratios are {d['agreement']['10']['mean_linear_ratio']:.6f} and {d['agreement']['100']['mean_linear_ratio']:.6f}. Neighboring zero-padded bins are correlated, so these are descriptive diagnostics, not independent sample counts or confidence tests.

For independent periodograms of a continuous Gaussian-like random spectrum, relative standard deviation is approximately 1/sqrt(K): 31.6% at K = 10 and 10.0% at K = 100 (about 1.37 and 0.43 dB under a small-error approximation). This heuristic is not a line-amplitude confidence interval and does not apply unmodified at spectral nulls, DC, deterministic lines, or correlated coded data. Finite windows broaden lines and fill apparent nulls through leakage; report a line's integrated power, not its plotted peak density alone.

## Pattern 1: independent random NRZ

![Pattern 1, K = 10](../results/pattern_1_K10.png)

![Pattern 1, K = 100](../results/pattern_1_K100.png)

Figures 1 and 2. Pattern 1 follows the random NRZ envelope, with noticeably reduced fluctuations at K = 100. The 10 GHz pulse null spacing and the 50 GHz rise-time null are consistent with Tb = 100 ps and Tr = 20 ps. Exact null depths depend on finite-window leakage and the plotted grid.

## Pattern 2: the intended 8b/10b-like record

![Pattern 2, K = 10](../results/pattern_2_K10.png)

![Pattern 2, K = 100](../results/pattern_2_K100.png)

Figures 3 and 4. Pattern 2 strongly suppresses low-frequency power relative to independent random NRZ. The ratio of average PSDs over 10-100 MHz at K = 100 is {d['lowratio']:.2f} dB relative to pattern 1. Transition and disparity constraints create correlations between bits; replacing these correlations by independent bits would lose the low-frequency spectral shaping. This identification is qualified by the supplied-file audit below.

## Pattern 3: random payload plus periodic suffix

![Pattern 3, K = 10](../results/pattern_3_K10.png)

![Pattern 3, K = 100](../results/pattern_3_K100.png)

Figures 5 and 6. A comb of lines at integer GHz is superposed on a random background. The aligned 10-bit block lasts 1 ns. The fixed 01 suffix creates a nonzero, periodic phase-conditioned mean even though its average over a block is zero. Missing lines at 10, 20, 30, 40 and 50 GHz agree with the suffix difference and pulse-shape zeros. The comb comes from the repeated mean structure; the entire 200,000-bit sequence is not periodic.

## Classification from the bit files

| Check | Pattern 1 | Pattern 2 | Pattern 3 |
|---|---:|---:|---:|
| Bits | 200,000 | 200,000 | 200,000 |
| Ones | 99,981 | 100,294 | 100,007 |
| Maximum equal-bit run | 17 | 6 | 9 |
| Final signed bit sum | -38 | 588 | 14 |
| Aligned distinct 10-bit words | 1,024 | 440 | 256 |

![Phase and cumulative-sum diagnostics](assets/pattern_diagnostics.png)

Pattern 1 uses all 1,024 10-bit words and has no strong position-dependent bias; its statistics and PSD are consistent with random NRZ. A finite record cannot prove stochastic independence. Pattern 3 has exactly 20,000 aligned words; position 9 is always 0 and position 10 is always 1. Its first eight positions are close to balanced and all 256 8-bit payload values occur. The exact suffix is proven by the file; an independent, equiprobable payload is the useful statistical model, not a proven generator specification.

Genuine 8b/10b does allow repeated bits. Its relevant properties are controlled DC balance, bounded running disparity and run-length limits, which permit AC coupling and aid clock recovery. The original IBM code gives a maximum run length of five and digital sum variation of six [3]; AMD describes the character-end running-disparity state as +1 or -1 [4]. Thus low DC content should be explained by controlled imbalance and correlation, not by an absence of consecutive zeros or ones.

The supplied pattern 2 contains 42 runs of six bits. Its cumulative signed bit sum ranges from -2 to 590 and ends at 588, inconsistent with the bounded imbalance of a strict 8b/10b stream. At the file's initial 10-bit alignment, 119 words have disparity +4. That alignment alone is not decisive: shifting the word start by six bits removes all absolute word disparities above two, but the long-term cumulative drift persists. The report therefore identifies pattern 2 as the intended coded example, with 8b/10b-like spectral shaping, and does not claim a standards-compliant encoding or a successful full decoder check. Input data are preserved unchanged.

## Exact hybrid model and spectral lines

Write the bipolar symbols as a random zero-mean part plus a periodic mean. Under the independent, equiprobable payload model, the ten phase-conditioned means are [0,0,0,0,0,0,0,0,-A,+A]. Only eight out of ten positions have variance A^2, so the time-averaged continuous spectrum is 0.8 times the independent random NRZ reference, or -0.9691 dB. A periodic variance alone is not a deterministic line source; the nonzero periodic mean supplies the lines here.

For fm = m/(10 Tb) = m GHz, the Fourier-series coefficient has magnitude given by

$$c_m=\\frac{{A}}{{10}}\\operatorname{{sinc}}(f_mT_b)\\operatorname{{sinc}}(f_mT_r)e^{{-j\\phi_m}}[e^{{-j2\\pi m9/10}}-e^{{-j2\\pi m8/10}}].$$

The phase phi_m collects the common NRZ-pulse and smoothing delays; it does not affect power. The total one-sided model is

$$S_{{hybrid,1}}(f)=0.8S_{{NRZ,1}}(f)+\\sum_{{m=1}}^\\infty\\frac{{2|c_m|^2}}{{Z_0}}\\delta(f-f_m).$$

Consequently each positive-frequency line has integrated power

$$P_m=\\frac{{8A^2}}{{100Z_0}}\\operatorname{{sinc}}^2(m/10)\\operatorname{{sinc}}^2(m/50)\\sin^2(\\pi m/10).$$

![Predicted positive-frequency line powers](assets/hybrid_predicted_line_powers.png)

The fundamental is 1 GHz. The suffix difference supplies sin(pi m/10), which vanishes at m = 10, 20, ... and at DC; the NRZ pulse also vanishes at these 10 GHz multiples, and the transition filter first vanishes at 50 GHz. These zeros explain why not every integer GHz line is present. A finite-length window convolves each ideal delta line with its spectral response. Plot peak heights thus depend on T and the window even when physical line power is fixed. The model describes the ideal ensemble; finite payload realizations fluctuate about the 0.8 background.

## Reproduction and limitations

From this homework directory, install the packages in requirements.txt and report/requirements.txt, run `python src/hw4.py`, run `python -m unittest discover -s tests -v`, run `python tests/verify_reference.py` for the independent SciPy audit, run `python report/build_report.py`, then run `python tests/verify_delivery.py` to check delivery integrity. The analysis reads only instructions/HW4_bit_pattern_1.txt through HW4_bit_pattern_3.txt. The report builder reads package-local results and those bit files, creates its own diagnostic assets, and emits this PDF and Markdown. Figures and numerical spectra retain no personal identifiers. Generated PDF metadata uses a generic title and an empty author field.

Independent verification passed all {d['verification']['unit_test_count']} tests with zero failures and zero errors. A separate construction convolves the complete bipolar stream with a causal boxcar, then splits the waveform and applies scipy.signal.periodogram with density normalization and no detrending. It does not call the rewritten estimator. Across all six saved spectra (262,145 bins each), the largest relative L2 error is {max(x['relative_L2_error'] for x in d['verification']['spectral_comparisons']):.3e}; the largest relative full-band power discrepancy is {max(x['relative_full_band_power_error'] for x in d['verification']['spectral_comparisons']):.3e}. Full details are in results/verification.json.

`results/summary.json` records every signal and window parameter, per-file statistics, integrated powers and Parseval errors. `results/spectra.npz` retains the full one-sided frequency grid and all six linear PSD arrays, as well as both reference curves. `results/pattern_audit.json` supplies a separate file audit including SHA-256 checksums and every 10-bit starting alignment. The six requested plots use the continuous-time reference; the sampled reference quantifies discretization error. No claim of a strict 8b/10b implementation is made for pattern 2, and the random payload's source distribution cannot be recovered uniquely from a finite record.

## References

1. Course assignment, *Homework: Spectrum of digital signals*, pp. 1-2. [Included assignment](../instructions/HW4_spectrum_of_digital_signals.pdf).
2. Course lecture, *Spectrum of digital signals*, pp. 2-8, 15-16 and 40-43. PSD definitions, periodogram averaging, periodic/random spectra and window gains. Lecture material is a read-only teaching reference and is not reproduced in this report.
3. A. X. Widmer and P. A. Franaszek, *Transmission code for high-speed fibre-optic data networks*, 17 March 1983. [IBM Research](https://research.ibm.com/publications/transmission-code-for-high-speed-fibre-optic-data-networks).
4. AMD, *Versal Adaptive SoC GTY and GTYP Transceivers Architecture Manual (AM002): Running Disparity*. [Official documentation](https://docs.amd.com/r/en-US/am002-versal-gty-transceivers/Running-Disparity).
5. *8b/10b encoding*, Wikipedia. [Assignment-requested overview](https://en.wikipedia.org/wiki/8b/10b_encoding). Claims about run length and disparity are additionally checked against [3]-[4].
"""
    MARKDOWN.write_text(md, encoding="utf-8")


def build_pdf(d):
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle("ReportTitle", fontName="Helvetica-Bold", fontSize=27, leading=31, textColor=NAVY, spaceAfter=14))
    styles.add(ParagraphStyle("Heading", fontName="Helvetica-Bold", fontSize=17, leading=21, textColor=NAVY, spaceAfter=12))
    styles.add(ParagraphStyle("Subhead", fontName="Helvetica-Bold", fontSize=10.3, leading=14, textColor=BLUE, spaceBefore=10, spaceAfter=5))
    styles.add(ParagraphStyle("Text", fontName="Helvetica", fontSize=9.3, leading=13.2, textColor=INK, spaceAfter=8))
    styles.add(ParagraphStyle("SmallText", fontName="Helvetica", fontSize=8.2, leading=11.4, textColor=INK, spaceAfter=6))
    styles.add(ParagraphStyle("Caption", fontName="Helvetica", fontSize=8.0, leading=10.7, textColor=MUTED, spaceAfter=7))
    styles.add(ParagraphStyle("TableText", fontName="Helvetica", fontSize=8.4, leading=11.1, textColor=INK))
    styles.add(ParagraphStyle("Kicker", fontName="Helvetica-Bold", fontSize=9, leading=13, textColor=BLUE, spaceAfter=8))
    story = []
    def para(text, kind="Text"):
        story.append(Paragraph(text, styles[kind]))
    def title(text):
        para(text, "Heading")
    def image(name, width=CONTENT, maxheight=None, result=False):
        path = RESULTS / name if result else ASSETS / name
        from reportlab.lib.utils import ImageReader
        iw, ih = ImageReader(str(path)).getSize()
        height = width * ih/iw
        if maxheight and height > maxheight:
            width *= maxheight / height
            height = maxheight
        img = Image(str(path), width, height)
        img.hAlign = "CENTER"
        story.append(img)
    def equation(name, maxheight=36):
        image(name + ".png", maxheight=maxheight)
        story.append(Spacer(1, 6))
    def table(rows, widths, small=False):
        cellstyle = styles["SmallText"] if small else styles["TableText"]
        headerstyle = ParagraphStyle("HeaderCell", parent=cellstyle, textColor=colors.white, fontName="Helvetica-Bold")
        converted = [[Paragraph(str(c), headerstyle if i == 0 else cellstyle) for c in row] for i, row in enumerate(rows)]
        t = Table(converted, colWidths=widths, hAlign="LEFT", repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,0), NAVY), ("TEXTCOLOR", (0,0), (-1,0), colors.white),
            ("ROWBACKGROUNDS", (0,1), (-1,-1), [PALE, colors.white]),
            ("VALIGN", (0,0), (-1,-1), "TOP"),
            ("LEFTPADDING", (0,0), (-1,-1), 8), ("RIGHTPADDING", (0,0), (-1,-1), 8),
            ("TOPPADDING", (0,0), (-1,-1), 6), ("BOTTOMPADDING", (0,0), (-1,-1), 5),
            ("LINEBELOW", (0,0), (-1,0), 0.5, NAVY),
        ]))
        story.append(t)
        story.append(Spacer(1,8))
    def page():
        story.append(PageBreak())

    der, cfg, patterns = d["der"], d["cfg"], d["patterns"]
    para("ANONYMOUS TECHNICAL REPORT", "Kicker")
    para("HW4: Spectrum of<br/>digital signals", "ReportTitle")
    para("One-sided PSD estimation, bit-pattern identification and reproducible verification.")
    para("Three supplied 200,000-bit files are converted to trapezoidal NRZ and analyzed with the averaged periodogram. Six PSD plots cover 0-50 GHz at K = 10 and K = 100. The analysis checks the estimator's power normalization and audits the actual bit files before interpreting their spectra.")
    table([
        ["File", "Evidence-based classification"],
        ["Pattern 1", "Consistent with independent random NRZ."],
        ["Pattern 2", "The intended 8b/10b-like example; low-frequency shaping is clear, but the supplied stream fails basic strict 8b/10b invariants."],
        ["Pattern 3", "Eight random-looking payload bits plus an exact 01 suffix in each 10-bit word; continuous background plus a 1 GHz comb."],
    ], [77, CONTENT-77])
    para("Signal and analysis parameters", "Subhead")
    table([
        ["Quantity", "Value"],
        ["Input / full duration", "200,000 bits per file / 20 us at 10 Gbit/s"],
        ["Bit period / voltage", "100 ps / -0.5 V and +0.5 V (1 V peak-to-peak)"],
        ["Linear transition / impedance", "20 ps / 50 ohm"],
        ["Sampling", "dt = 2 ps; 50 samples per bit; 10 per transition; Nyquist = 250 GHz"],
        ["Segment", "2,000 bits = 100,000 samples = 200 ns"],
        ["Window", f"Kaiser beta = pi (alpha = 1); ENBW = {der['window_enbw_hz']/1e6:.6f} MHz"],
        ["FFT / grid", f"M = 524,288; df = {der['fft_grid_step_hz']/1e6:.9f} MHz"],
        ["Resolution / averages", "1/T = 5 MHz; K = 10 and 100; no segment overlap"],
    ], [147, CONTENT-147])
    para("K = 10 uses the first 20,000 bits (2 us); K = 100 uses all 200,000 bits (20 us). Both use the same segment length and window, and share the first 10 segments. The input is 200,000 bits, correcting the 2,000,000-bit statement in the earlier report.", "SmallText")
    page()

    title("Estimator, normalization and checks")
    para("Map each bit to bipolar level +/-A, A = 0.5 V. At a transition, start at the preceding level and reach the new level after 20 ps. The actual preceding bit is retained across segment boundaries; only the full record's start assumes preceding bit 0. No detrending or mean subtraction is applied.")
    equation("estimator", 57)
    para("Here N = 100,000 is the segment sample count, M is the zero-padded FFT length, and q<sub>j</sub> = 2 at interior positive frequencies. At DC and the even FFT's Nyquist bin q<sub>j</sub> = 1. Average squared magnitudes; averaging complex FFT coefficients would cancel random phases. Stored spectra use W/Hz; figures use 10 log<sub>10</sub>(1000 S) dBm/Hz.")
    para(f"The window power gain is U = mean(w<sup>2</sup>) = {der['window_power_gain_U']:.8f}; coherent gain G<sub>c</sub> = mean(w) = {der['window_coherent_gain']:.8f}. Broadband PSD uses U. G<sub>c</sub> instead corrects a bin-centered tone's amplitude. ENBW = (U/G<sub>c</sub><sup>2</sup>)/T = {der['window_enbw_hz']/1e6:.6f} MHz. Both ENBW and 1/T are below 10 MHz. Zero padding gives the required sub-1-MHz display grid; it does not improve physical resolving power.")
    equation("reference", 34)
    para("This is the assignment's continuous-time, one-sided independent-bit reference, with normalized sinc(x) = sin(pi x)/(pi x), for f &gt; 0. The single DC grid point uses the undoubled value. It is a reference envelope, not the theory for coded or periodic data. The sampled waveform's exact reference differs by sinc(f dt)<sup>-4</sup> away from removable zeros: a limiting ratio of only " + f"{der['sampled_vs_continuous_reference_at_50_GHz_db']:.4f} dB at 50 GHz.", "SmallText")
    para("Full-band power conservation", "Subhead")
    rows = [["Pattern / K", "Weighted time power", "PSD bin-sum power", "Relative error"]]
    for p in patterns:
        for k in (10,100):
            e = p["estimates"][str(k)]
            rows.append([f"{p['pattern']} / {k}", f"{e['mean_window_weighted_power_w']*1000:.6f} mW", f"{e['integrated_full_band_psd_w']*1000:.6f} mW", f"{e['parseval_relative_error']:.2e}"])
    table(rows,[84,137,137,CONTENT-358],small=True)
    para("Over 0-250 GHz, sum(S<sub>1</sub>) df equals the window-weighted time-domain power, including DC. Use the DFT bin sum, not trapezoidal endpoint half weights. The largest relative discrepancy above is " + f"{d['worst_parseval']:.2e}.", "SmallText")
    para(f"Pattern 1, K = 100: ordinary time power {d['raw_random']:.6f} mW; sampled independent-bit expectation {d['theory_random']:.6f} mW. Finite ramps reduce the ideal 5 mW. Over 0.1-9 GHz, RMS pointwise log error versus the sampled reference decreases from {d['agreement']['10']['rms_db']:.3f} to {d['agreement']['100']['rms_db']:.3f} dB. Zero-padded neighboring bins are correlated; these are descriptive checks.", "SmallText")
    page()

    figure_text = {
        1: ("Pattern 1: independent random NRZ", "The broad envelope follows the independent random NRZ reference. Larger K reduces pointwise fluctuations while preserving the 200 ns segment's resolving power. Null spacing reflects the 100 ps UI; the transition filter first vanishes at 50 GHz."),
        2: ("Pattern 2: intended coded example", f"Low-frequency power is suppressed; the average PSD over 10-100 MHz at K = 100 is {d['lowratio']:.2f} dB relative to pattern 1. Coding-related correlations explain the shaping. The bit audit on page 6 qualifies any strict 8b/10b claim."),
        3: ("Pattern 3: hybrid random / periodic", "The periodic 01 suffix adds integer-GHz lines to a random background. The 10-bit word lasts 1 ns, so its mean repeats at 1 GHz. Missing lines at 10 GHz multiples follow from the deterministic suffix and pulse zeros; the complete bitstream is not periodic."),
    }
    for p in (1,2,3):
        title(figure_text[p][0])
        para(figure_text[p][1], "SmallText")
        for j,k in enumerate((10,100)):
            image(f"pattern_{p}_K{k}.png", width=CONTENT-8, maxheight=273, result=True)
            para(f"Figure {2*p-1+j}. Pattern {p}, K = {k}; T = 200 ns. Blue: estimated one-sided PSD. Red: continuous-time independent random NRZ reference.", "Caption")
        page()

    title("Bit-file evidence and the 8b/10b caveat")
    rows = [["Audit (200,000 bits each)", "Pattern 1", "Pattern 2", "Pattern 3"],
            ["Ones", "99,981", "100,294", "100,007"],
            ["Maximum equal-bit run", "17", "6", "9"],
            ["Final signed bit sum", "-38", "588", "14"],
            ["Distinct aligned 10-bit words", "1,024", "440", "256"]]
    table(rows,[CONTENT*0.43, CONTENT*0.19, CONTENT*0.19, CONTENT*0.19])
    image("pattern_diagnostics.png", maxheight=176)
    para("Figure 7. Left: fraction of ones by aligned word position. Right: cumulative sum of +1 for a one and -1 for a zero (display downsampled). Exact extrema and counts come from the full record.", "Caption")
    para("Pattern 1 uses all 1,024 10-bit words and has no strong phase bias. Together with its PSD, this supports independent random NRZ; a finite sample cannot prove independence. Pattern 3 has 20,000 words with position 9 always 0 and position 10 always 1. All 256 possible 8-bit payload values occur. The suffix is exact; independent, equiprobable payload bits are a statistical model.")
    para("Why genuine 8b/10b suppresses low frequencies", "Subhead")
    para("8b/10b allows repeated bits. The relevant properties are controlled DC balance, bounded running disparity and constrained run length. These correlations reduce slow imbalance, permit AC coupling and aid clock recovery. The IBM code specifies maximum run length five and digital sum variation six [3]. AMD describes character-end running-disparity states of +1 or -1 [4]. The explanation is controlled imbalance, not an absence of consecutive zeros or ones.")
    para("What the actual supplied pattern 2 proves", "Subhead")
    para("Pattern 2 contains <b>42 runs of six bits</b>. Its cumulative signed bit sum spans <b>-2 to 590</b> and ends at <b>588</b>, inconsistent with a strict bounded-disparity 8b/10b stream. At the initial 10-bit alignment, 119 words have disparity +4. Alignment alone is not decisive: offset six bits removes all absolute word disparities above two, but the cumulative drift persists.")
    para("The evidence supports the assignment's intended coded example and 8b/10b-like spectral shaping. It does not validate a standards-compliant encoder or a full decoder check. The supplied bit files are preserved unchanged.", "SmallText")
    page()

    title("Hybrid model: background plus a comb")
    para("Under the independent, equiprobable payload model, the ten phase-conditioned bipolar means are [0, 0, 0, 0, 0, 0, 0, 0, -A, +A]. The average mean is zero, but the phase-dependent mean is nonzero and periodic. Eight positions out of ten have variance A<sup>2</sup>, so the time-averaged continuous spectrum is 0.8 times the random NRZ reference (-0.9691 dB). The periodic mean supplies the deterministic lines.")
    para("Let f<sub>m</sub> = m/(10T<sub>b</sub>) = m GHz. The mean waveform's Fourier-series coefficients are", "SmallText")
    equation("coefficient", 35)
    para("The common phase phi<sub>m</sub> contains the pulse and smoothing delays; it does not affect power. The one-sided model is", "SmallText")
    equation("hybrid", 35)
    para("Each positive-frequency line therefore has integrated power", "SmallText")
    equation("linepower", 35)
    image("hybrid_predicted_line_powers.png", maxheight=192)
    para("Figure 8. Analytical integrated positive-frequency line powers, not finite-window PSD peak heights. The random background is omitted from this diagnostic.", "Caption")
    para("The 1 GHz fundamental follows from a 1 ns mean period. The suffix difference gives sin(pi m/10), zero at DC and m = 10, 20, ...; the NRZ pulse also vanishes at those 10 GHz multiples. The transition filter first vanishes at 50 GHz. This explains both the comb spacing and missing harmonics.")
    para("A finite window broadens each ideal delta line. A line's displayed dBm/Hz peak depends on T and the window even when its physical integrated power is fixed. Finite payload realizations fluctuate about the 0.8 continuum; observing all payload values does not prove an independent random source.", "SmallText")
    page()

    title("Uncertainty, reproduction and references")
    para("Uncertainty and interpretation", "Subhead")
    para("For independent periodograms of a Gaussian-like continuous random spectrum, the relative standard deviation is approximately 1/sqrt(K): 31.6% at K = 10 and 10.0% at K = 100. The small-error dB approximations are 1.37 and 0.43 dB. This heuristic is not a line-amplitude confidence interval and does not apply unchanged to nulls, DC, deterministic lines or correlated coded data.")
    para("Window leakage limits apparent null depth and broadens narrow lines. Neighboring zero-padded bins are correlated. The continuous reference assumes an ideal pulse model; the sampled reference quantifies the small discretization difference. The bit audit identifies the supplied pattern 2's strict-code inconsistencies; no source data are corrected to hide them.")
    para("Reproduce this report", "Subhead")
    table([
        ["Step", "Command from the homework directory"],
        ["Install dependencies", "python -m pip install -r requirements.txt -r report/requirements.txt"],
        ["Regenerate analysis", "python src/hw4.py"],
        ["Run verification", "python -m unittest discover -s tests -v"],
        ["Independent SciPy audit", "python tests/verify_reference.py"],
        ["Rebuild PDF and Markdown", "python report/build_report.py"],
        ["Delivery integrity check", "python tests/verify_delivery.py"],
    ], [124, CONTENT-124])
    para(f"All {d['verification']['unit_test_count']} tests pass. A separate waveform construction and SciPy density periodogram reproduce all six saved spectra (262,145 bins each) without calling the rewritten estimator. Worst relative L2 error: {max(x['relative_L2_error'] for x in d['verification']['spectral_comparisons']):.3e}; worst full-band power discrepancy: {max(x['relative_full_band_power_error'] for x in d['verification']['spectral_comparisons']):.3e}. Detailed results: results/verification.json.", "SmallText")
    para("summary.json stores configuration, window gains, bit statistics and six power checks. spectra.npz stores the full frequency grid, all linear PSD estimates and both random-reference curves. pattern_audit.json stores a separate audit including SHA-256 checksums and all ten word alignments. All inputs and generated paths are package-relative. The PDF title is generic and the author metadata is empty.", "SmallText")
    para("References", "Subhead")
    references = [
        '[1] Course assignment, <i>Homework: Spectrum of digital signals</i>, pp. 1-2. Included as instructions/HW4_spectrum_of_digital_signals.pdf.',
        '[2] Course lecture, <i>Spectrum of digital signals</i>, pp. 2-8, 15-16 and 40-43. PSD definitions, averaging, periodic/random spectra and window gains. Used as a read-only teaching reference.',
        '[3] A. X. Widmer and P. A. Franaszek, <i>Transmission code for high-speed fibre-optic data networks</i>, 17 March 1983. <link href="https://research.ibm.com/publications/transmission-code-for-high-speed-fibre-optic-data-networks" color="#1558a6">IBM Research: original code reference</link>.',
        '[4] AMD, <i>Versal Adaptive SoC GTY and GTYP Transceivers Architecture Manual (AM002): Running Disparity</i>. <link href="https://docs.amd.com/r/en-US/am002-versal-gty-transceivers/Running-Disparity" color="#1558a6">Official running-disparity documentation</link>.',
        '[5] <i>8b/10b encoding</i>, Wikipedia. <link href="https://en.wikipedia.org/wiki/8b/10b_encoding" color="#1558a6">Assignment-requested overview</link>. Run-length and disparity claims are additionally checked against [3]-[4].',
    ]
    for ref in references:
        para(ref,"SmallText")

    class AnonymousCanvas(canvas.Canvas):
        def __init__(self,*args,**kwargs):
            kwargs["invariant"] = 1
            super().__init__(*args,**kwargs)
            self.setTitle("HW4: Spectrum of digital signals - anonymous report")
            self.setAuthor("")
            self.setSubject("One-sided averaged periodogram and bit-pattern verification")
            self.setCreator("Anonymous reproducible report builder")

    def decoration(c, doc):
        c.saveState()
        c.setFont("Helvetica-Bold",8)
        c.setFillColor(BLUE)
        c.drawString(52, HEIGHT-31,"HW04  |  SPECTRUM OF DIGITAL SIGNALS")
        c.setStrokeColor(colors.HexColor("#c9d8e5"))
        c.setLineWidth(0.6)
        c.line(52,HEIGHT-41,WIDTH-52,HEIGHT-41)
        c.line(52,39,WIDTH-52,39)
        c.setFillColor(MUTED)
        c.setFont("Helvetica",7.4)
        c.drawString(52,26,"Anonymous technical report")
        c.drawRightString(WIDTH-52,26,str(doc.page))
        c.restoreState()
    doc = SimpleDocTemplate(str(PDF), pagesize=A4, rightMargin=52, leftMargin=52, topMargin=58, bottomMargin=51,
                            title="HW4: Spectrum of digital signals - anonymous report", author="")
    doc.build(story,onFirstPage=decoration,onLaterPages=decoration,canvasmaker=AnonymousCanvas)


def main():
    summary, audit, spectra = load_inputs()
    agreement = make_assets(summary,spectra)
    data = make_text(summary,audit,agreement)
    write_markdown(data)
    build_pdf(data)
    checks = {"random_reference_diagnostics": agreement,
              "worst_full_band_parseval_relative_error": data["worst_parseval"],
              "pattern_2_vs_1_mean_psd_10_to_100_MHz_db": data["lowratio"]}
    (HERE / "report_metrics.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    print(f"Created {PDF.name} and {MARKDOWN.name}")


if __name__ == "__main__":
    main()
