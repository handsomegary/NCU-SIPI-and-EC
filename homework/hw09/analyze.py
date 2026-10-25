"""HW9: verified, fixed-rate 0/1 decision-feedback equalization.

All waveforms, cursor values and taps use the same matched transfer S21/2.
The main simulation constructs finite-rise NRZ and uses linear overlap-add.
Tap/phase/threshold design uses the exact single-bit response of that discrete
model. A separate scalar decision recurrence verifies the held-block feedback.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.signal import fftconvolve

# Support both direct execution and importlib-based independent tests.
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from channel import (read_touchstone, frequency_impulse, get_sbr, trapezoid_pulse,
                     levels_to_waveform, sectional_convolve, eye_histogram)

TB = 25e-12
TR = 5e-12
DT = 0.5e-12
VPP = 2.0
ALPHA = 0.5
T2 = 20e-9
P = 50


def cursor_model(sbr, samples_per_ui, latency):
    """All finite-model cursor offsets: positive j is a postcursor."""
    sbr = np.asarray(sbr, dtype=float)
    p, latency = int(samples_per_ui), int(latency)
    if sbr.ndim != 1 or not len(sbr) or not np.all(np.isfinite(sbr)):
        raise ValueError("SBR must be a finite nonempty vector")
    if p < 1 or not 0 <= latency < len(sbr):
        raise ValueError("Invalid UI or full sampling latency")
    offsets = np.arange(-(latency // p), (len(sbr) - 1 - latency) // p + 1)
    return offsets, sbr[latency + offsets * p]


def design_at(sbr, samples_per_ui, latency, ntaps):
    """Correct-history 0/1 residual-ISI bounds, including all residual cursors."""
    if ntaps < 0 or int(ntaps) != ntaps:
        raise ValueError("Tap count must be a nonnegative integer")
    p, latency = int(samples_per_ui), int(latency)
    offsets, cursors = cursor_model(sbr, p, latency)
    post = np.arange(1, int(ntaps) + 1)
    taps = [-float(cursors[offsets == j][0]) if np.any(offsets == j) else 0.0
            for j in post]
    residual = cursors.copy()
    residual[np.isin(offsets, post)] = 0.0
    main = float(cursors[offsets == 0][0])
    residual = residual[offsets != 0]
    low = float(np.maximum(residual, 0).sum())
    high = float(main + np.minimum(residual, 0).sum())
    return {
        "tap_count": int(ntaps),
        "latency_samples": latency,
        "phase_samples": latency % p,
        "main_cursor_v": main,
        "taps_v": taps,
        "threshold_v": (low + high) / 2,
        "low_upper_v": low,
        "high_lower_v": high,
        "worst_case_eye_v": high - low,
        "residual_abs_sum_v": float(np.abs(residual).sum()),
        "residual_sum_v": float(residual.sum()),
        "cursor_offset_min": int(offsets[0]),
        "cursor_offset_max": int(offsets[-1]),
    }


def design_dfe(sbr, samples_per_ui, ntaps):
    """Search one UI around the SBR maximum; first maximum breaks ties."""
    p = int(samples_per_ui)
    if p < 1 or len(sbr) < p:
        raise ValueError("SBR must span at least one valid UI")
    start = min(max(int(np.argmax(sbr)) - p // 2, 0), len(sbr) - p)
    candidates = [design_at(sbr, p, lag, ntaps) for lag in range(start, start + p)]
    best = candidates[int(np.argmax([x["worst_case_eye_v"] for x in candidates]))]
    return dict(best), candidates


def _validate_dfe(y, samples_per_ui, taps, threshold):
    y, taps = np.asarray(y, float), np.asarray(taps, float)
    if y.ndim != 1 or taps.ndim != 1 or not np.all(np.isfinite(y)) or not np.all(np.isfinite(taps)):
        raise ValueError("Waveform and taps must be finite one-dimensional vectors")
    if int(samples_per_ui) != samples_per_ui or samples_per_ui < 2:
        raise ValueError("At least two integral samples per UI are required")
    if not np.isfinite(threshold):
        raise ValueError("Threshold must be finite")
    return y, taps, int(samples_per_ui)


def fixed_rate_dfe(y, samples_per_ui, taps, threshold, phase):
    """Return summer waveform, sample indices and 0/1 decisions.

    Tap c[0] from decision at k is held over [k+P//2, k+P//2+P).
    Thus the present sample is never modified by its own decision. Initial
    pre-record decision history is zero. Feedback always uses detected bits.
    """
    y, taps, p = _validate_dfe(y, samples_per_ui, taps, threshold)
    if int(phase) != phase:
        raise ValueError("Sampling phase must be an integer")
    indices = np.arange(int(phase) % p, len(y), p)
    decisions = np.zeros(len(indices), dtype=np.uint8)
    z = y.copy()
    for position, k in enumerate(indices):
        decision = z[k] >= threshold
        decisions[position] = decision
        if decision:
            for i, tap in enumerate(taps):
                start = k + p // 2 + i * p
                z[start:start + p] += tap
    return z, indices, decisions


def sample_reference(samples, taps, threshold, flip_decision_at=None):
    """Independent scalar, past-decision dot product; no waveform block writes."""
    samples, taps = np.asarray(samples, float), np.asarray(taps, float)
    if samples.ndim != 1 or taps.ndim != 1 or not np.all(np.isfinite(samples)) or not np.all(np.isfinite(taps)) or not np.isfinite(threshold):
        raise ValueError("Finite one-dimensional samples/taps and threshold required")
    corrected = np.empty(len(samples), dtype=float)
    detected = np.zeros(len(samples), dtype=np.uint8)
    for k, value in enumerate(samples):
        count = min(k, len(taps))
        feedback = float(np.dot(taps[:count], detected[k-count:k][::-1]))
        corrected[k] = value + feedback
        detected[k] = corrected[k] >= threshold
        if k == flip_decision_at:
            detected[k] = 1 - detected[k]
    return corrected, detected


def conditional_metrics(values, truth, threshold, decisions=None):
    values, truth = np.asarray(values, float), np.asarray(truth, np.uint8)
    predicted = (values >= threshold).astype(np.uint8) if decisions is None else np.asarray(decisions)
    low = float(values[truth == 0].max())
    high = float(values[truth == 1].min())
    wrong = np.flatnonzero(predicted != truth)
    return {"evaluated_bits": int(len(truth)), "errors": int(len(wrong)),
            "observed_error_fraction": float(len(wrong) / len(truth)),
            "observed_low_max_v": low, "observed_high_min_v": high,
            "observed_eye_v": high-low,
            "minimum_signed_margin_v": float(np.min(np.where(truth == 1, values-threshold, threshold-values))),
            "first_error_payload_bit": int(wrong[0]) if len(wrong) else None}


def save_csv(path, header, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _plot_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .22,
                         "figure.dpi": 130, "savefig.dpi": 170,
                         "axes.spines.top": False, "axes.spines.right": False})
    return plt


def _save_figure(plt, fig, output):
    fig.savefig(output, bbox_inches="tight", metadata={"Software": "HW9 reproducible analysis"})
    plt.close(fig)


def plot_eye(plt, z, start, bit_count, cfg, output):
    p = P
    # Start half a UI before the first scored sampling instant: sample at 0.5 UI.
    segment = z[start-p//2:start-p//2 + bit_count*p]
    prob, extent = eye_histogram(segment, p, 420)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    image = ax.imshow(np.log10(np.maximum(np.hstack([prob, prob]), 1e-6)),
                      origin="lower", extent=[0, 2, extent[2], extent[3]],
                      aspect="auto", cmap="viridis", vmin=-5, vmax=-1)
    ax.axhline(cfg["threshold_v"], color="white", ls="--", lw=1, label="Threshold")
    ax.axvline(.5, color="white", ls=":", lw=.9)
    ax.axvline(1.5, color="white", ls=":", lw=.9)
    ax.set(xlabel="Time (UI; 1 UI = 25 ps)", ylabel="Summer voltage (V)",
           title=f'{cfg["tap_count"]}-tap DFE: summer-output statistical eye')
    ax.grid(False)
    fig.colorbar(image, ax=ax, label="log10 probability per voltage bin")
    fig.tight_layout()
    _save_figure(plt, fig, output)
    return {"histogram_samples": int(len(segment)),
            "column_normalization_max_error": float(np.max(np.abs(prob.sum(axis=0)-1))),
            "voltage_bins": 420}


def plot_snapshot(plt, y, z, sample_indices, truth, decisions, threshold, first_sample, filename, out, title):
    start = int(first_sample - P//2)
    count = round(1e-9 / DT)
    ix = np.arange(start, start + count + 1)
    tps = (ix - start) * DT * 1e12
    selected = (sample_indices >= start) & (sample_indices <= ix[-1])
    sample_t = (sample_indices[selected] - start) * DT * 1e12
    vals = z[sample_indices[selected]]
    bad = decisions[selected] != truth[selected]
    fig, ax = plt.subplots(figsize=(9.2, 3.8))
    ax.plot(tps, y[ix], color="#738596", lw=1.15, label="Post-channel")
    ax.plot(tps, z[ix], color="#007F88", lw=1.3, label="Post-DFE summer")
    ax.axhline(threshold, color="#A96121", ls="--", lw=1, label="Threshold")
    ax.scatter(sample_t[~bad], vals[~bad], s=12, color="#007F88", label="Correct decisions", zorder=3)
    if np.any(bad):
        ax.scatter(sample_t[bad], vals[bad], marker="x", s=35, color="#B42A35", label="Wrong decisions", zorder=4)
    ax.set(xlim=(0, 1000), xlabel="Time from snapshot start (ps)", ylabel="Voltage (V)", title=title)
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    fig.tight_layout()
    _save_figure(plt, fig, out/f"{filename}.png")
    save_csv(out/f"{filename}.csv", ["time_ps", "post_channel_v", "post_dfe_v"],
             zip(tps, y[ix], z[ix]))
    return {"file": f"{filename}.png", "csv": f"{filename}.csv",
            "start_sample": start, "start_time_s": float(start*DT),
            "samples": int(len(ix)), "duration_ns": float((ix[-1]-ix[0])*DT*1e9),
            "scored_decisions_in_snapshot": int(selected.sum()),
            "errors_in_snapshot": int(bad.sum())}


def run(bits_count=65536, seed=20261025, output=None):
    if bits_count < 1024:
        raise ValueError("Use at least 1024 payload bits")
    out = Path(output) if output else HERE/"results"
    out.mkdir(parents=True, exist_ok=True)
    plt = _plot_setup()
    data = HERE/"data/HW7_stripline_50cm.s2p"
    f, s, z0 = read_touchstone(data)
    H = s[:, 1, 0] / 2
    t, h_period = frequency_impulse(f, H, DT, ALPHA)
    nh = round(T2/DT)
    h = h_period[:nh]
    pulse = trapezoid_pulse(TB, TR, DT, VPP)
    g = fftconvolve(pulse, h)
    ta, analytic = get_sbr(f, H, TB, TR, DT, VPP, ALPHA)
    configs, searches = {}, {}
    for n in (1, 2, 3):
        configs[n], searches[n] = design_dfe(g, P, n)
    guard = int(np.ceil(len(g)/P)) + 8
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, bits_count+2*guard, dtype=np.uint8)
    truth = bits[guard:guard+bits_count]
    x = levels_to_waveform(np.r_[VPP*bits.astype(float), 0.0], P, round(TR/DT))
    print(f"Channel convolution: {bits_count} payload bits, {guard} guard bits each side", flush=True)
    y = sectional_convolve(x, h, section_samples=16384)
    del x
    # Separate full-record convolution validates the block overlap/add placement.
    short_x = levels_to_waveform(np.r_[VPP*bits[:193].astype(float), 0.0], P, round(TR/DT))
    short_ola = sectional_convolve(short_x, h, section_samples=1031)
    short_fft = fftconvolve(short_x, h)
    verification = {
        "sectional_vs_whole_fft_max_abs_v": float(np.max(np.abs(short_ola-short_fft))),
        "analytic_vs_discrete_sbr_peak_region_max_abs_v": float(np.max(np.abs(g[:len(h)]-analytic[:len(h)]))),
        "finite_fir_dc_gain": float(h.sum()),
        "frequency_dc_gain": float(H[0].real),
        "finite_fir_dc_difference": float(h.sum()-H[0].real),
        "input_data_sha256": hashlib.sha256(data.read_bytes()).hexdigest(),
    }
    # Quantify model sensitivity to extending FIR truncation from 20 ns to 40 ns.
    g40 = fftconvolve(pulse, h_period[:round(40e-9/DT)])
    sensitivity = {}
    for n in (1,2,3):
        longer = design_at(g40, P, configs[n]["latency_samples"], n)
        sensitivity[str(n)] = {
            "fixed_latency_threshold_change_v": longer["threshold_v"]-configs[n]["threshold_v"],
            "fixed_latency_eye_bound_change_v": longer["worst_case_eye_v"]-configs[n]["worst_case_eye_v"],
            "max_tap_change_v": float(np.max(np.abs(np.array(longer["taps_v"])-configs[n]["taps_v"])))
        }
    verification["truncation_20_to_40_ns"] = sensitivity
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.7))
    axes[0].plot(f/1e9, 20*np.log10(np.maximum(np.abs(H),1e-15)), color="#007F88")
    axes[0].set(xlabel="Frequency (GHz)", ylabel="Matched transfer magnitude (dB)", title="H = S21 / 2")
    axes[1].plot(np.arange(len(g))*DT*1e9, g, color="#007F88", label="Discrete model")
    axes[1].plot(ta*1e9, analytic, ls="--", lw=.9, color="#A96121", label="Analytic pulse reference")
    peak_t = np.argmax(g)*DT*1e9
    axes[1].set(xlim=(peak_t-.12, peak_t+.4), xlabel="Time (ns)", ylabel="SBR (V)", title="Causal finite-rise single-bit response")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    _save_figure(plt, fig, out/"channel_and_sbr.png")
    fig, ax = plt.subplots(figsize=(7.5,3.7))
    for n in (1,2,3):
        candidates=searches[n]
        ax.plot([c["latency_samples"]*DT*1e12 for c in candidates],
                [c["worst_case_eye_v"] for c in candidates], label=f"{n} tap(s)")
        ax.scatter([configs[n]["latency_samples"]*DT*1e12], [configs[n]["worst_case_eye_v"]], s=25)
    ax.set(xlabel="Full source-to-sample latency (ps)", ylabel="Correct-history worst-case eye (V)", title="Automatic sampling search over one UI")
    ax.legend()
    fig.tight_layout()
    _save_figure(plt,fig,out/"phase_search.png")
    cases=[]
    for n in (1,2,3):
        cfg=configs[n]
        latency=cfg["latency_samples"]
        scored=(guard+np.arange(bits_count))*P+latency
        z, indices, decisions=fixed_rate_dfe(y,P,cfg["taps_v"],cfg["threshold_v"],cfg["phase_samples"])
        scalar, scalar_decisions=sample_reference(y[indices],cfg["taps_v"],cfg["threshold_v"])
        locations=(scored-cfg["phase_samples"])//P
        predicted=decisions[locations]
        comparison=float(np.max(np.abs(z[indices]-scalar)))
        mismatch=int(np.count_nonzero(decisions!=scalar_decisions))
        offsets, values=cursor_model(g,P,latency)
        reconstructed=np.convolve(bits.astype(float),values)[guard+np.arange(bits_count)-offsets[0]]
        reconstruction=float(np.max(np.abs(reconstructed-y[scored])))
        metrics=conditional_metrics(z[scored],truth,cfg["threshold_v"],predicted)
        rawdesign=design_at(g,P,latency,0)
        rawmetrics=conditional_metrics(y[scored],truth,rawdesign["threshold_v"])
        check={
            "waveform_vs_scalar_max_abs_v":comparison,
            "waveform_vs_scalar_decision_mismatches":mismatch,
            "source_bits_cursor_reconstruction_max_abs_v":reconstruction,
            "correct_history_bound_available":bool(cfg["worst_case_eye_v"]>0),
            "observed_gap_minus_bound_v":metrics["observed_eye_v"]-cfg["worst_case_eye_v"],
        }
        histcheck=plot_eye(plt,z,scored[0],bits_count,cfg,out/f"eye_{n}tap.png")
        check.update(histcheck)
        first_error=metrics["first_error_payload_bit"]
        snap_bit=min(max((first_error if first_error is not None else 512)-10,0),bits_count-41)
        snapshot=plot_snapshot(plt,y,z,scored,truth,predicted,cfg["threshold_v"],scored[snap_bit],
                               f"waveform_{n}tap",out,f"{n}-tap DFE: exactly 1 ns waveform comparison")
        save_csv(out/f"decisions_{n}tap.csv",
                 ["payload_bit","truth","decision","post_channel_v","post_dfe_v","threshold_v"],
                 zip(range(bits_count),truth,predicted,y[scored],z[scored],np.full(bits_count,cfg["threshold_v"])))
        save_csv(out/f"phase_search_{n}tap.csv",
                 ["latency_samples","phase_samples","threshold_v","worst_case_eye_v","main_cursor_v"],
                 ([c[k] for k in ("latency_samples","phase_samples","threshold_v","worst_case_eye_v","main_cursor_v")] for c in searches[n]))
        case={"design":cfg,"metrics":metrics,"raw_channel":rawmetrics,
              "raw_channel_threshold_v":rawdesign["threshold_v"],"verification":check,"snapshot":snapshot,
              "eye_file":f"eye_{n}tap.png"}
        cases.append(case)
        print(f'{n} tap(s): phase={cfg["phase_samples"]}/{P}, Vt={cfg["threshold_v"]:.6f}, errors={metrics["errors"]}/{bits_count}',flush=True)
        if n==3:
            # A separately labeled deliberate threshold mistuning study.
            stress_threshold=cfg["threshold_v"]+.6*cfg["main_cursor_v"]
            stress_z, stress_ix, stress_decisions=fixed_rate_dfe(y,P,cfg["taps_v"],stress_threshold,cfg["phase_samples"])
            stress_predictions=stress_decisions[locations]
            stressmetrics=conditional_metrics(stress_z[scored],truth,stress_threshold,stress_predictions)
            e=stressmetrics["first_error_payload_bit"]
            stress_example=None
            if e is not None:
                stress_example={"payload_bit":e,"true_bit":int(truth[e]),"decision":int(stress_predictions[e]),
                                "post_channel_v":float(y[scored[e]]),"summer_v":float(stress_z[scored[e]]),
                                "threshold_v":stress_threshold,"source_context_start":max(0,e-5),
                                "source_context":truth[max(0,e-5):e+6].astype(int).tolist()}
            stress_snap_bit=min(max((e if e is not None else 512)-10,0),bits_count-41)
            stress_snapshot=plot_snapshot(plt,y,stress_z,scored,truth,stress_predictions,stress_threshold,
                                          scored[stress_snap_bit],"stress_threshold",out,
                                          "Deliberately mistuned threshold: 1 ns example")
            stress={"label":"Deliberate threshold mistuning, separate from baseline",
                    "threshold_v":stress_threshold,"metrics":stressmetrics,"first_error":stress_example,"snapshot":stress_snapshot}
            # Isolate error propagation by flipping one baseline decision.
            inject_payload=512
            inject_global=int(locations[inject_payload])
            inj_samples, inj_decisions=sample_reference(y[indices],cfg["taps_v"],cfg["threshold_v"],inject_global)
            rows=[]
            for q in range(inject_payload,inject_payload+12):
                loc=locations[q]
                rows.append([q,int(truth[q]),int(predicted[q]),int(inj_decisions[loc]),
                             float(scalar[loc]),float(inj_samples[loc]),float(inj_samples[loc]-scalar[loc])])
            save_csv(out/"injected_error.csv",
                     ["payload_bit","truth","baseline_decision","injected_run_decision",
                      "baseline_summer_v","injected_summer_v","voltage_difference_v"],rows)
            injection={"label":"One forced decision flip; artificial diagnostic",
                       "payload_bit":inject_payload,"forced_from":int(predicted[inject_payload]),
                       "forced_to":int(1-predicted[inject_payload]),
                       "extra_decision_changes_after_injection":int(np.count_nonzero(inj_decisions[locations[inject_payload+1:]]!=predicted[inject_payload+1:])),
                       "next_three_voltage_changes_v":[r[-1] for r in rows[1:4]],
                       "csv":"injected_error.csv"}
    summary={
        "assignment":"HW9 DFE simulation and design",
        "model":{"data_rate_gbps":40,"ui_ps":25,"rise_time_ps":5,"source_vpp_v":2,
                 "sample_step_ps":.5,"samples_per_ui":P,"channel_transfer":"S21/2","reference_impedance_ohm":z0,
                 "frequency_points":len(f),"frequency_step_mhz":float((f[1]-f[0])/1e6),
                 "max_frequency_ghz":float(f[-1]/1e9),"kaiser_alpha":ALPHA,"impulse_truncation_ns":20,
                 "payload_bits":bits_count,"guard_bits_each_side":guard,"seed":seed,
                 "sequence":"Seeded independent Bernoulli(0.5) bits, not an LFSR PRBS",
                 "decision_alphabet":[0,1],"initial_history":0,
                 "noise":"None; no jitter or component nonidealities modeled"},
        "method":{"phase_search":"Exactly 50 consecutive sample latencies centered near the SBR maximum; first maximum wins",
                  "threshold":"Midpoint between correct-history worst-case 0 upper bound and 1 lower bound",
                  "tap_units":"V per logical 1; additive feedback c_j=-a_j; no tap normalization",
                  "truth_alignment":"Source bit k sampled at k*50 + full_latency; phase alone is insufficient",
                  "eye_node":"Summer output, before slicer"},
        "verification":verification,"cases":cases,"stress":stress,"injected_error":injection,
        "limitations":["Worst-case bounds are conditional on correct previous decisions.",
                       "Observed error fraction is a finite noiseless record result, not a physical BER estimate.",
                       "FIR is truncated at 20 ns; spectrum above 80 GHz is zero; Kaiser alpha is 0.5.",
                       "No ADS ChannelSim comparison was performed (optional assignment part)."]
    }
    # Fail the run rather than publish internally inconsistent results.
    assert verification["sectional_vs_whole_fft_max_abs_v"]<1e-10
    for case in cases:
        v=case["verification"]
        assert v["waveform_vs_scalar_max_abs_v"]<1e-10
        assert v["waveform_vs_scalar_decision_mismatches"]==0
        assert v["source_bits_cursor_reconstruction_max_abs_v"]<1e-10
        assert v["column_normalization_max_error"]<1e-12
        assert case["snapshot"]["samples"]==2001
    write_json(out/"summary.json",summary)
    save_csv(out/"metrics.csv",["taps","phase_samples","latency_ps","threshold_v","worst_case_eye_v","observed_eye_v","errors","bits"],
             ([c["design"]["tap_count"],c["design"]["phase_samples"],c["design"]["latency_samples"]*.5,
               c["design"]["threshold_v"],c["design"]["worst_case_eye_v"],c["metrics"]["observed_eye_v"],
               c["metrics"]["errors"],bits_count] for c in cases))
    print(json.dumps({"output":out.name,"cases":len(cases),"verification":"passed"},indent=2),flush=True)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits",type=int,default=65536,help="Payload bits; seeded guards added automatically")
    parser.add_argument("--seed",type=int,default=20261025)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    run(args.bits,args.seed,args.output)


if __name__=="__main__":
    main()
