"""Build the anonymous HW11 report from the verified solver's saved outputs.

Run from any working directory: python report/build_report.py
No original personal report or instructor example image is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer,
    Table, TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
NAVY = colors.HexColor("#16324A")
TEAL = colors.HexColor("#087F8C")
GRAY = colors.HexColor("#52616B")
PALE = colors.HexColor("#EEF5F8")
rl_config.invariant = 1


def num(value, precision=5):
    return f"{float(value):.{precision}g}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "report" / "HW11_Crosstalk_Anonymous_Report.pdf")
    args = parser.parse_args()
    data = json.loads((args.results / "metrics.json").read_text(encoding="utf-8"))
    validation = json.loads((ROOT / "tests" / "verification.json").read_text(encoding="utf-8"))
    if validation["status"] != "passed":
        raise ValueError("Independent verification must pass before building the final report")
    if validation["retained_bits_per_transmitter"] != data["configuration"]["bits"]:
        raise ValueError("Verification and analysis records have different bit counts")
    analysis_hash = hashlib.sha256((ROOT / "analyze.py").read_bytes()).hexdigest()
    if analysis_hash != data["reproducibility"]["analysis_source_sha256"]:
        raise ValueError("Analysis source changed after results were generated")
    if analysis_hash != validation["analysis_source_sha256"]:
        raise ValueError("Verification does not match the current analysis source")
    for artifact in ("metrics.json", "channel.npz", "simulation.npz"):
        digest = hashlib.sha256((args.results / artifact).read_bytes()).hexdigest()
        if digest != validation["audited_artifact_sha256"][artifact]:
            raise ValueError(f"Verification is stale for {artifact}; rerun tests/verify_results.py")
    cfg, inp = data["configuration"], data["input"]
    figures = args.results / "figures"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    work_pdf = args.output.with_suffix(".working.pdf")

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle("ReportTitle", fontName="Helvetica-Bold", fontSize=27,
                              leading=32, textColor=NAVY, spaceAfter=16))
    styles.add(ParagraphStyle("ReportSubtitle", fontName="Helvetica", fontSize=13,
                              leading=18, textColor=GRAY, spaceAfter=18))
    styles.add(ParagraphStyle("Section", fontName="Helvetica-Bold", fontSize=18,
                              leading=22, textColor=NAVY, spaceAfter=14))
    styles.add(ParagraphStyle("Subsection", fontName="Helvetica-Bold", fontSize=11,
                              leading=15, textColor=TEAL, spaceBefore=10, spaceAfter=6))
    styles.add(ParagraphStyle("BodyReport", fontName="Helvetica", fontSize=9.5,
                              leading=14, textColor=NAVY, spaceAfter=9))
    styles.add(ParagraphStyle("CaptionReport", fontName="Helvetica", fontSize=8.1,
                              leading=11, textColor=GRAY, spaceBefore=5, spaceAfter=12))
    styles.add(ParagraphStyle("SmallReport", fontName="Helvetica", fontSize=8,
                              leading=11, textColor=NAVY, spaceAfter=5))
    styles.add(ParagraphStyle("TableHeaderReport", parent=styles["SmallReport"],
                              textColor=colors.white, fontName="Helvetica-Bold"))
    styles.add(ParagraphStyle("EquationReport", fontName="Courier", fontSize=8.5,
                              leading=14, textColor=NAVY, leftIndent=10,
                              backColor=PALE, borderPadding=8, spaceBefore=6, spaceAfter=12))
    story = []
    width = A4[0] - 84

    def p(text, style="BodyReport"):
        if style == "EquationReport":
            text = text.replace(" ", "&nbsp;")
        return Paragraph(text, styles[style])

    def add(text, style="BodyReport"):
        story.append(p(text, style))

    def title(text):
        add(text, "Section")

    def table(rows, col_widths=None):
        cells = [[p(str(cell), "TableHeaderReport" if index == 0 else "SmallReport")
                  for cell in row] for index, row in enumerate(rows)]
        out = Table(cells, colWidths=col_widths, hAlign="LEFT", repeatRows=1)
        out.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), NAVY),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("LEFTPADDING", (0, 0), (-1, -1), 8),
            ("RIGHTPADDING", (0, 0), (-1, -1), 8),
            ("TOPPADDING", (0, 0), (-1, -1), 7),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
            ("LINEBELOW", (0, 0), (-1, 0), .5, TEAL),
        ]))
        story.append(out)
        story.append(Spacer(1, 9))

    def figure(name, caption, max_height=330):
        path = figures / name
        if not path.is_file():
            raise FileNotFoundError(f"Required regenerated figure missing: {name}")
        image = Image(str(path))
        scale = min(width / image.imageWidth, max_height / image.imageHeight)
        image.drawWidth, image.drawHeight = image.imageWidth * scale, image.imageHeight * scale
        image.hAlign = "CENTER"
        story.append(image)
        add(caption, "CaptionReport")

    def opening(case, rx):
        return data["cases"][case][f"rx{rx}"]["opening_v"]

    # Page 1: the conclusion is tied to actual saved results.
    add("HW11", "ReportSubtitle")
    add("Crosstalk simulation<br/>of eight coupled channels", "ReportTitle")
    add("Anonymous technical report | Reproducible frequency-to-time-domain analysis",
        "ReportSubtitle")
    add("Scope and result", "Subsection")
    add("The supplied 16-port network is terminated by the specified source and load RC networks. "
        "A completed Python solver computes all transmitter-to-receiver transfer functions, "
        "then combines eight independently generated NRZ signals using linear convolution. "
        "The results below use the regenerated run saved with this package, rather than the "
        "example images in the assignment.")
    table([
        ["Empirical voltage-cluster quantile gap", "Rx0", "Rx4"],
        ["Only corresponding transmitter on", f"{num(opening('only_tx0', 0))} V",
         f"{num(opening('only_tx4', 4))} V"],
        ["All eight transmitters on", f"{num(opening('all_on', 0))} V",
         f"{num(opening('all_on', 4))} V"],
        ["Decrease from isolated case", f"{num(opening('only_tx0', 0)-opening('all_on', 0))} V",
         f"{num(opening('only_tx4', 4)-opening('all_on', 4))} V"],
    ], [width*.48, width*.26, width*.26])
    add("This diagnostic summarizes empirical voltage-cluster separation over a finite "
        "deterministic record. It is not a bit-conditioned or worst-case eye opening, "
        "bit-error-rate estimate, or receiver compliance limit. The exact definition "
        "and numerical checks are given on the following pages.")
    add("Deliverables", "Subsection")
    add("The report includes all sixteen S<sub>1j</sub> curves, the isolated and all-on eyes at "
        "Rx0 and Rx4, the requested 4 ns Rx1-Rx4 waveform snapshot with only Tx4 active, "
        "and all eight all-on receiver eyes. The complete source, saved numerical summaries, "
        "and validation tests are supplied alongside this document.")
    add("Optional ADS comparison", "Subsection")
    add("ADS was not run for this report. The supplied assignment includes instructor example "
        "ADS images; no numerical agreement with a fresh ADS simulation is claimed here.")
    story.append(PageBreak())

    # Page 2: the specification is explicit and uses unambiguous indices.
    title("1 | Problem definition")
    add("Physical Touchstone ports are numbered 1-16. The program's channel labels are "
        "zero-based: Tx k drives physical port 2k+1 and Rx k observes physical port 2k+2, "
        "for k = 0,...,7. Thus Tx4 is physical port 9 and Rx4 is physical port 10.")
    table([
        ["Parameter", "Specified / simulated value"],
        ["Signaling", "10 Gbit/s single-ended NRZ; UI = 100 ps; source voltage levels 0 and 1 V"],
        ["Rise and fall transitions", "20 ps linear transitions"],
        ["Source termination", f"R<sub>S</sub> = {num(cfg['rs_ohm'])} ohm; C<sub>S</sub> = {num(cfg['cs_f']*1e12)} pF shunt at each Tx port"],
        ["Receiver termination", f"R<sub>L</sub> = {num(cfg['rl_ohm'])} ohm; C<sub>L</sub> = {num(cfg['cl_f']*1e12)} pF in parallel at each Rx port"],
        ["Initial source delays", "Tx0, Tx2, Tx4, Tx6: 0 ps; Tx1, Tx3, Tx5, Tx7: 50 ps (0.5 UI)"],
        ["Input data", f"16-port S network; {inp['frequency_count']:,} frequency samples; {num(inp['df_hz']/1e6)} MHz spacing; maximum {num(inp['fmax_hz']/1e9)} GHz"],
        ["Reference impedance", f"Common real Z<sub>0</sub> = {num(inp['z0_ohm'])} ohm"],
        ["Retained record", f"{cfg['bits']:,} bits per transmitter after the steady-state guard"],
    ], [width*.31, width*.69])
    add("Required-output map", "Subsection")
    table([
        ["Assignment item", "Location in this report"],
        ["1. S11, S12, ..., S1,16", "Section 4"],
        ["2-3. Rx0 with only Tx0 / all Tx on", "Section 5, top row"],
        ["4-5. Rx4 with only Tx4 / all Tx on", "Section 5, bottom row"],
        ["6. Rx1-Rx4 snapshot with only Tx4 on", "Section 6"],
        ["Compute the eyes at every Rx", "Appendix A, all eight all-on eyes"],
    ], [width*.60, width*.40])
    add("Inactive transmitters", "Subsection")
    add("Turning a transmitter off sets its source waveform to zero while retaining its "
        "specified source resistor and shunt capacitor. Removing inactive terminations would "
        "change the circuit and would not be the isolated-transmitter comparison used here.")
    story.append(PageBreak())

    # Page 3: model derivation, not an unexplained black-box calculation.
    title("2 | Multiport circuit model")
    add("For each frequency, voltage-wave scattering uses v<sup>-</sup> = S v<sup>+</sup>. "
        "With port currents defined into the S-parameter network and a common real reference "
        "impedance, v = v<sup>+</sup> + v<sup>-</sup> and Z<sub>0</sub>i = "
        "v<sup>+</sup> - v<sup>-</sup>. Eliminating the wave amplitudes gives:")
    add("(I - S) v - Z0 (I + S) i = 0", "EquationReport")
    add("A source v<sub>Tj</sub> in series with R<sub>S</sub> is replaced by a Norton "
        "current v<sub>Tj</sub>/R<sub>S</sub> in parallel with R<sub>S</sub>. "
        "The external port admittance is G + j omega C, with diagonal entries 1/R<sub>S</sub> "
        "and C<sub>S</sub> at Tx ports, and 1/R<sub>L</sub> and C<sub>L</sub> at Rx ports. "
        "KCL and the scattering relation form the 32-by-32 system:")
    add("[ G+j*omega*C       I             ] [v]   [B]\n<br/>"
        "[ I-S          -Z0*(I+S)         ] [i] = [0]", "EquationReport")
    add("The eight columns of B excite one transmitter each. Column j has 1/R<sub>S</sub> "
        "at Tx j's physical port and zero elsewhere, so the voltage part of the solution is "
        "already normalized to a 1 V Thevenin source. Selecting the eight Rx rows yields "
        "H<sub>rj</sub>(f) = V<sub>Rx r</sub>(f) / V<sub>Tj</sub>(f). "
        "All right-hand sides are solved together for each frequency; a matrix inverse is "
        "not explicitly formed.")
    add("y_r(t) = sum_j [ h_rj(t) * x_j(t) ]", "EquationReport")
    add("The network and external RC terminations are linear and time invariant. Therefore "
        "one may compute each transmitter's response at all receivers and add them. "
        "The stored transfer-function axes are frequency, receiver, transmitter; selecting "
        "a transmitter means selecting the third axis.")
    add("Why S/2 alone is insufficient here", "Subsection")
    add("With all ports resistively matched to Z<sub>0</sub> and no extra capacitors, the "
        "off-diagonal source-to-load transfer is S<sub>ij</sub>/2. This useful normalization "
        "check is not the actual HW11 circuit: the 30 ohm / 65 ohm terminations and shunt "
        "capacitors require the loaded multiport solve above.")
    add("Coupling interpretation", "Subsection")
    add("A diagonal channel H<sub>rr</sub> contains its loaded transmission response. "
        "Off-diagonal H<sub>rj</sub> terms are crosstalk paths. Nearby aggressors contribute "
        "the largest interference in this data set; their 0.5 UI source offset places many "
        "transitions near another channel's decision region. The channel's own intersymbol "
        "interference and the added aggressor terms must both be retained.")
    story.append(PageBreak())

    # Page 4: numerical choices and scope of numerical verification.
    title("3 | Frequency-to-time method")
    add("The supplied S-parameter matrix is read without assuming reciprocity and "
        "checked for uniform frequency spacing. The loaded H spectrum is multiplied by "
        "the Kaiser window below over the supplied 0-40 GHz band, then zero-padded to "
        "the 250 GHz sampling Nyquist frequency. DC is constrained to be real. "
        "Hermitian negative-frequency symmetry and a 50,000-point real inverse DFT "
        "produce dimensionless discrete-convolution weights.")
    add("W(f) = I0(pi*a*sqrt(1-(f/fmax)^2)) / I0(pi*a)", "EquationReport")
    table([
        ["Numerical setting", "Saved run"],
        ["Samples per UI", str(cfg['samples_per_ui'])],
        ["Time step", f"{num(cfg['dt_s']*1e12)} ps"],
        ["Kaiser frequency-window parameter", str(cfg['kaiser_a'])],
        ["Retained negative impulse time", f"{num(cfg['pre_time_s']*1e9)} ns"],
        ["Retained positive impulse time", f"{num(cfg['post_time_s']*1e9)} ns"],
        ["Steady-state guard / convolution section", f"{cfg['guard_ui']} UI / {cfg['block_size']:,} samples"],
        ["NRZ random streams", "Separate deterministic SeedSequence child streams, one per Tx"],
        ["Master seed", str(cfg['seed'])],
    ], [width*.52, width*.48])
    add("A finite frequency bandwidth produces pre-ringing and long impulse tails. The "
        "solver retains the specified negative- and positive-time support rather than "
        "discarding all negative lags or circularly wrapping the response into the data. "
        "Input bits are generated for a guard on both sides of the retained record. "
        "Finite 20 ps linear edges and the alternating 50 ps source delays are applied "
        "to each transmitter before convolution.")
    add("Sectional convolution", "Subsection")
    add("Linear convolution is evaluated with FFT overlap-add. Section responses are "
        "added at their correct absolute offsets; padding prevents within-section circular "
        "aliasing. Only the guarded central steady-state interval is used for eyes and "
        "metrics. Separate child streams avoid the accidental identical-pattern "
        "correlation that would result from resetting every transmitter to one seed.")
    add("Eye and voltage-cluster diagnostic", "Subsection")
    add("Each two-UI eye is a histogram of voltage samples folded by the fixed UI. "
        "Each phase column is normalized to observed bin mass; the color scale is "
        "log<sub>10</sub>(probability per voltage bin), not probability density or BER. "
        "The eye phase is absolute sample time modulo 100 ps; the data are not "
        "recentered or aligned to recovered bits.")
    add("The saved <font face='Courier'>opening_v</font> diagnostic splits voltage samples "
        "at half the main-path DC gain (about 0.340664 V). At each of 50 phases, it "
        "subtracts the lower cluster's 99.9% quantile from the upper cluster's "
        "0.1% quantile, requiring at least 20% of samples in each cluster. "
        "The largest gap is reported. The histogram has 3,000 bins from -0.25 to "
        "1.25 V, giving 0.5 mV voltage resolution. This voltage-conditioned "
        "diagnostic is not bit-conditioned and cannot establish a target-BER "
        "eye opening or detect all possible error trajectories.")
    add("The bandwidth, window, impulse support, sampling, and finite record define "
        "this approximation. No random jitter, receiver noise, clock recovery, or "
        "equalizer optimization is included.", "SmallReport")
    story.append(PageBreak())

    title("4 | S-parameter responses")
    figure("s_parameters.png", "Figure 1. Regenerated magnitude responses for all sixteen "
           "S<sub>1j</sub> elements. Indices here are the original one-based physical "
           "Touchstone port indices. These curves belong to the bare network, before "
           "the external source and load RC networks are added.", max_height=440)
    add("The S<sub>1j</sub> request is a row of the multiport matrix. It must not be "
        "confused with the eight loaded source-to-Rx responses. In particular, "
        "S<sub>12</sub> connects the two ends of channel 0, while S<sub>11</sub> is "
        "the reflection at physical port 1. The simulation uses the full matrix "
        "throughout, including all coupling and reflections.")
    add("The physical channel port pairs are (1,2), (3,4), (5,6), (7,8), (9,10), "
        "(11,12), (13,14), and (15,16). Pairing adjacent entries in a block of "
        "odd-numbered ports would therefore be incorrect.")
    story.append(PageBreak())

    title("5 | Required receiver-eye comparisons")
    image_rows = []
    for names in [("eye_rx0_only_tx0.png", "eye_rx0_all_on.png"),
                  ("eye_rx4_only_tx4.png", "eye_rx4_all_on.png")]:
        row = []
        for name in names:
            image = Image(str(figures / name))
            scale = min((width-10)/2/image.imageWidth, 240/image.imageHeight)
            image.drawWidth, image.drawHeight = image.imageWidth*scale, image.imageHeight*scale
            row.append(image)
        image_rows.append(row)
    eyes_table = Table(image_rows, colWidths=[width/2, width/2], hAlign="CENTER")
    eyes_table.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                                   ("LEFTPADDING", (0, 0), (-1, -1), 0),
                                   ("RIGHTPADDING", (0, 0), (-1, -1), 0),
                                   ("BOTTOMPADDING", (0, 0), (-1, -1), 12)]))
    story.append(eyes_table)
    add("Figure 2. Top: Rx0 with only Tx0 active (left) and all transmitters active "
        "(right). Bottom: Rx4 with only Tx4 active (left) and all transmitters active "
        "(right). All plots are regenerated from this run with the same numerical "
        "settings and transmitter data reused between each isolated/all-on comparison.", "CaptionReport")
    add("The source voltage levels remain 0 and 1 V; the receiver levels differ because "
        "of loaded transmission loss and reflections. Adding aggressors increases the "
        "range of trajectories at a given phase and reduces the voltage-cluster gap. "
        "Rx4 lies within the group of eight channels, so coupling arrives from both "
        "sides; Rx0 is an edge channel.")
    add("The display represents finite sample occupancy. Unoccupied histogram bins "
        "mean that the simulated record did not visit those bins. This does not "
        "establish a zero error probability.")
    story.append(PageBreak())

    title("6 | Tx4-only waveform snapshot")
    figure("tx4_only_snapshot.png", "Figure 3. Requested 4 ns snapshot at Rx1, Rx2, "
           "Rx3, and Rx4 with only Tx4 active. The plot explicitly identifies any "
           "5x display scaling of Rx1-Rx3 crosstalk; Rx4 is the directly driven "
           "channel. Display scaling is not part of the simulation or eye metrics.", max_height=430)
    add("The response at Rx4 carries the main NRZ pattern. The responses at the "
        "other receivers are caused by coupling from Tx4, with no independently "
        "driven data on those receivers' corresponding transmitters. They show "
        "transition-sensitive pulses and ringing rather than a separate NRZ data stream.")
    add("Receiver labels are zero-based. Rx1, Rx2, Rx3, and Rx4 correspond to physical "
        "ports 4, 6, 8, and 10, respectively; the active source Tx4 is physical port 9. "
        "The closest displayed coupled receiver is Rx3.")
    add("The snapshot is taken from the guarded steady-state output, not the startup "
        "transient. Independent random patterns are required for the all-on test; "
        "this isolated snapshot intentionally has only one nonzero source waveform.")
    story.append(PageBreak())

    title("Appendix A | All eight receiver eyes")
    figure("eyes_all_on.png", "Figure 4. All eight all-on eyes computed from the "
           "same retained record and the same absolute time-modulo-UI reference. "
           "The appendix satisfies the requirement to compute all receiver eyes, "
           "including those not highlighted in the main comparisons.", max_height=410)
    metric_rows = [["Receiver", "Physical Rx port", "Cluster gap (V)", "Best phase (UI)"]]
    for rx in range(8):
        met = data["cases"]["all_on"][f"rx{rx}"]
        metric_rows.append([f"Rx{rx}", str(2*rx+2), num(met["opening_v"]), num(met["best_phase_ui"])])
    table(metric_rows, [width*.17, width*.24, width*.31, width*.28])
    story.append(PageBreak())

    title("Appendix B | Numerical verification")
    add("The solver is accompanied by independent circuit and convolution checks. "
        "The saved verification summary records the numerical outcomes; tests "
        "exercise matched-load normalization, independent circuit equations, "
        "port ordering, source independence, linear superposition, and finite "
        "convolution alignment. Truncation diagnostics and a larger-support "
        "short-record comparison quantify numerical approximation rather than "
        "claiming an exact infinite-bandwidth waveform.")
    add("Recorded numerical diagnostics", "Subsection")
    diag, conv = data["impulse_diagnostics"], data["convergence"]
    sampling = data["sampling_convergence"]
    def maximum_abs(value):
        if isinstance(value, list):
            return max(maximum_abs(v) for v in value)
        return abs(float(value))
    table([
        ["Quantity", "Measured value"],
        ["DFT length / frequency-grid period", f"{diag['dft_length']:,} / {num(diag['dft_period_s']*1e9)} ns"],
        ["Maximum discarded impulse L1 norm, one path", num(maximum_abs(diag['discarded_l1_by_path']), 7)],
        ["Maximum discarded impulse L2 norm, one path", num(maximum_abs(diag['discarded_l2_by_path']), 7)],
        ["Maximum retained-tap DC-sum error", f"{num(diag['dc_sum_max_abs_error']*1e3, 7)} mV per 1 V source"],
        ["Negative-lag energy fraction", num(diag['negative_lag_energy_fraction'], 7)],
        ["Larger support for comparison", f"-{num(conv['reference_pre_time_s']*1e9)} ns to +{num(conv['reference_post_time_s']*1e9)} ns"],
        ["Comparison record, all eight Tx active", f"{conv['test_bits']:,} bits per transmitter"],
        ["Maximum waveform difference, larger support", f"{num(conv['maximum_abs_waveform_difference_v']*1e3, 7)} mV"],
        ["RMS waveform difference, larger support", f"{num(conv['rms_waveform_difference_v']*1e3, 7)} mV"],
        ["Sampling comparison, all eight Tx active", f"{sampling['coarse_samples_per_ui']} vs {sampling['fine_samples_per_ui']} samples/UI; {sampling['test_bits']} bits"],
        ["Maximum / RMS sampling difference", f"{num(sampling['maximum_abs_waveform_difference_v']*1e3, 7)} / {num(sampling['rms_waveform_difference_v']*1e3, 7)} mV"],
        ["Simulation arithmetic", "Float32 FFT and accumulation; float64 circuit/diagnostics"],
    ], [width*.62, width*.38])
    add("The larger-support comparison checks this particular short all-on sequence; "
        "it is not a universal worst-case bound. The residual DC-sum error explicitly "
        "records a small support-truncation bias. Kernel samples are dimensionless "
        "discrete-convolution weights, so no additional time-step multiplier is applied.", "SmallReport")
    add("Independent checks: passed", "Subsection")
    unit = validation["unit_tests"]
    add(f"All {unit['run']} independent unit tests passed with {unit['failures']} failures "
        f"and {unit['errors']} errors. Saved-artifact verification additionally checked "
        f"all {validation['loaded_transfer_checked_frequencies']:,} frequencies and "
        f"{validation['loaded_transfer_checked_paths']} loaded transfer paths.", "SmallReport")
    table([
        ["Independent comparison", "Maximum absolute error"],
        ["MNA vs wave-boundary circuit equations", num(validation['max_mna_vs_wave_boundary_abs_error'], 7)],
        ["Saved impulse vs direct trigonometric sum", num(validation['max_saved_impulse_vs_direct_trigonometric_sum_error'], 7)],
        ["Float64 overlap-add vs direct convolution", num(validation['max_float64_overlap_add_vs_direct_convolution_error'], 7)],
        ["Saved Tx4 snapshot vs direct convolution", f"{num(validation['max_tx4_snapshot_vs_direct_convolution_error_v']*1e6, 7)} microV"],
        ["Maximum off-diagonal source-bit correlation", num(validation['maximum_absolute_offdiagonal_source_bit_correlation'], 7)],
    ], [width*.65, width*.35])
    story.append(PageBreak())
    title("Appendix C | Reproduction and source record")
    add("Audit scope", "Subsection")
    add("Every all-on receiver histogram conserves 13,107,200 samples, or 262,144 "
        "samples per phase, with no samples outside the histogram voltage range. "
        "Ten desired-channel percentile gaps were independently recalculated. "
        "The entire full-record waveform is not stored or independently reconstructed; "
        "the direct waveform check covers the saved snapshot.")
    add("Reproduce the complete saved run", "Subsection")
    add("Use Python 3.12 or newer and install <font face='Courier'>requirements.txt</font>. "
        "In the HW11 directory, run <font face='Courier'>python verify_inputs.py</font>, "
        "<font face='Courier'>python analyze.py</font>, "
        "<font face='Courier'>python tests/verify_results.py</font>, then "
        "<font face='Courier'>python report/build_report.py</font>. "
        "The artifact audit includes the unit tests. The report builder requires matching "
        "hashes for the analysis source and all three saved numerical artifacts.")
    add("The explicit analysis command, run from the repository root, is shown "
        "wrapped below; enter it as one line:")
    add("python homework/hw11/analyze.py --bits 262144 --seed 2026<br/>"
        "  --samples-per-ui 50 --block-size 32768<br/>"
        "  --kaiser-a 1 --pre-ns 5 --post-ns 20", "EquationReport")
    rep = data["reproducibility"]
    table([["Recorded environment", "Version"],
           ["Python", rep["python_version"]], ["NumPy", rep["numpy_version"]],
           ["SciPy", rep["scipy_version"]], ["Matplotlib", rep["matplotlib_version"]]],
          [width*.5, width*.5])
    add("The master seed and separate child-stream construction are fixed. "
        "The software versions are recorded because exact pseudorandom and floating-point "
        "reproduction can depend on the numerical environment. "
        "The source hash below binds the saved results to the analysis code.")
    add("Source record", "Subsection")
    add("The authoritative assignment is <font face='Courier'>instructions/hw11_assignment.pdf</font>. "
        "The supplied template, clarification notes, network data, and legacy support files "
        "are retained as references. The manifest records original and staged hashes. "
        "The incomplete template and obsolete manual same-seed advice are not the executable "
        "workflow. The prior personal student report is excluded.")
    add("Dataset SHA-256: " + str(inp.get("sha256", inp.get("SHA256", "see source_manifest.json"))), "SmallReport")
    add("Analysis SHA-256: " + rep["analysis_source_sha256"], "SmallReport")

    def page(canvas, doc):
        canvas.saveState()
        w, h = A4
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(1)
        canvas.line(42, h-31, w-42, h-31)
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(GRAY)
        canvas.drawString(42, 25, "HW11 | Eight-channel crosstalk simulation | Anonymous")
        canvas.drawRightString(w-42, 25, f"{doc.page}")
        canvas.restoreState()

    document = SimpleDocTemplate(str(work_pdf), pagesize=A4, rightMargin=42,
                                 leftMargin=42, topMargin=48, bottomMargin=43,
                                 title="HW11 Crosstalk Simulation - Anonymous Technical Report",
                                 author="", creator="HW11 reproducible report builder",
                                 subject="Eight coupled channels with RC terminations")
    document.build(story, onFirstPage=page, onLaterPages=page)
    # Keep only deliberately neutral metadata, including no Author field.
    reader = PdfReader(work_pdf)
    writer = PdfWriter()
    for sheet in reader.pages:
        writer.add_page(sheet)
    writer.add_metadata({"/Title": "HW11 Crosstalk Simulation - Anonymous Technical Report",
                         "/Subject": "Eight coupled channels with RC terminations",
                         "/Creator": "HW11 reproducible report builder",
                         "/Producer": "ReportLab and pypdf"})
    with args.output.open("wb") as stream:
        writer.write(stream)
    work_pdf.unlink()
    print(f"Created {args.output.name}: {len(reader.pages)} pages")


if __name__ == "__main__":
    main()
