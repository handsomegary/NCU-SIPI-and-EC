"""HW13: empirical jitter statistics and reproducible two-Gaussian tail fits.

Times are in ps and densities in 1/ps throughout this module. The primary input
is the verified HW12 eye-slice export. Extrapolated BER is a model prediction,
not an observed error rate or a confidence bound.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.special import log_ndtr, logsumexp

BASE = Path(__file__).resolve().parents[1]
CASES = ("case1_no_tx_jitter", "case2_sj", "case3_rj", "case4_sj_rj")
LABELS = ("No Tx jitter", "SJ", "RJ", "SJ + RJ")
UI_PS = 25.0
RHO = 0.5
TI = 12.5
TIP1 = 37.5
TAIL_MIN = 3e-4
TAIL_MAX = 3e-2
TARGETS = (1e-6, 1e-9, 1e-12, 1e-15)
COLORS = ("#2563a6", "#bf6500", "#16816a", "#8456a6")


@dataclass
class Histogram:
    centers: np.ndarray
    pmf: np.ndarray
    observations: float | None = None

    def __post_init__(self):
        self.centers = np.asarray(self.centers, dtype=float)
        self.pmf = np.asarray(self.pmf, dtype=float)
        if (self.centers.ndim != 1 or self.centers.size < 2
                or self.pmf.shape != self.centers.shape
                or not np.all(np.isfinite(self.centers))
                or not np.all(np.isfinite(self.pmf)) or np.any(self.pmf < 0)):
            raise ValueError("Invalid finite, nonnegative histogram")
        spacing = np.diff(self.centers)
        self.dt = float(np.median(spacing))
        if self.dt <= 0 or not np.allclose(spacing, self.dt, rtol=1e-10, atol=1e-12):
            raise ValueError("Histogram needs a uniform increasing time grid")
        total = float(self.pmf.sum())
        if not np.isfinite(total) or total <= 0:
            raise ValueError("Histogram must have positive mass")
        self.pmf = self.pmf / total
        self.edges = np.r_[self.centers - self.dt / 2, self.centers[-1] + self.dt / 2]
        self.cumulative = np.r_[0.0, np.cumsum(self.pmf)]
        self.survival = np.r_[np.cumsum(self.pmf[::-1])[::-1], 0.0]
        self.cumulative[-1] = 1.0
        self.survival[0] = 1.0

    @property
    def pdf(self):
        return self.pmf / self.dt

    def cdf(self, t):
        return np.interp(t, self.edges, self.cumulative, left=0.0, right=1.0)

    def sf(self, t):
        # Reverse summation retains tiny right-tail masses; never use 1-CDF.
        return np.interp(t, self.edges, self.survival, left=1.0, right=0.0)

    def ber(self, t, rho=RHO, ti=TI, tip1=TIP1):
        t = np.asarray(t)
        return rho * (self.sf(t - ti) + self.cdf(t - tip1))


def validate_parameters(p):
    p = np.asarray(p, dtype=float)
    if (p.shape != (5,) or not np.all(np.isfinite(p))
            or p[1] <= 0 or p[3] <= 0 or not 0 <= p[4] <= 1):
        raise ValueError("Expected five finite parameters, positive sigmas and weight in [0,1]")
    return p


def log_weights(w):
    return (float(np.log(w)) if w > 0 else -np.inf,
            float(np.log1p(-w)) if w < 1 else -np.inf)


def log_pdf(t, p):
    """p = (mu1, sigma1, mu2, sigma2, weight1); weight2 = 1-weight1."""
    m1, s1, m2, s2, w = validate_parameters(p)
    lw1, lw2 = log_weights(w)
    t = np.asarray(t)
    return logsumexp(np.stack((lw1 - np.log(s1) - .5 * ((t - m1) / s1) ** 2,
                              lw2 - np.log(s2) - .5 * ((t - m2) / s2) ** 2)),
                     axis=0) - .5 * np.log(2 * np.pi)


def log_tail(t, p, left=True):
    m1, s1, m2, s2, w = validate_parameters(p)
    lw1, lw2 = log_weights(w)
    t = np.asarray(t)
    sign = 1.0 if left else -1.0
    return np.logaddexp(lw1 + log_ndtr(sign * (t - m1) / s1),
                        lw2 + log_ndtr(sign * (t - m2) / s2))


def log_ber(t, p, rho=RHO, ti=TI, tip1=TIP1):
    if not np.all(np.isfinite([rho, ti, tip1])) or not 0 < rho <= 1 or tip1 <= ti:
        raise ValueError("Invalid transition density or bit interval")
    t = np.asarray(t)
    return np.log(rho) + np.logaddexp(log_tail(t - ti, p, False),
                                     log_tail(t - tip1, p, True))


def model_ber(t, p, **kwargs):
    return np.exp(log_ber(t, p, **kwargs))


def tail_masks(hist, low, high):
    if not 0 < low < high < 0.5:
        raise ValueError("Tail limits must satisfy 0 < low < high < 0.5")
    left = (hist.cumulative[1:] >= low) & (hist.cumulative[1:] <= high)
    right = (hist.survival[:-1] >= low) & (hist.survival[:-1] <= high)
    if min(left.sum(), right.sum()) < 3 or left.sum() + right.sum() < 8:
        raise ValueError("Insufficient tail samples for a five-parameter fit")
    return left, right


def fit_tails(hist, low=TAIL_MIN, high=TAIL_MAX):
    left, right = tail_masks(hist, low, high)
    xleft, xright = hist.edges[1:][left], hist.edges[:-1][right]
    yleft, yright = hist.cumulative[1:][left], hist.survival[:-1][right]

    def residual(p):
        return np.r_[log_tail(xleft, p) - np.log(yleft),
                     log_tail(xright, p, False) - np.log(yright)] / np.log(10)

    trials = []
    for weight in (.15, .5, .85):
        for shift in (1.0, 3.0):
            result = least_squares(residual, [-shift, 1.5, shift, 1.5, weight],
                bounds=([-8, .2, 0, .2, .01], [0, 5, 8, 5, .99]),
                max_nfev=3000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
            if result.success and np.all(np.isfinite(result.fun)):
                trials.append(result)
    if not trials:
        raise RuntimeError("No converged tail fit")
    best = min(trials, key=lambda f: float(np.dot(f.fun, f.fun)))
    p = best.x
    pdfmask = (left | right) & (hist.pmf > 0)
    pdf_error = (log_pdf(hist.centers[pdfmask], p) - np.log(hist.pdf[pdfmask])) / np.log(10)
    details = {
        "parameters": dict(zip(("mu1_ps", "sigma1_ps", "mu2_ps", "sigma2_ps", "weight1"), map(float, p))),
        "weight2": float(1 - p[4]), "parameter_vector": p.tolist(),
        "tail_min": low, "tail_max": high, "starts": 6,
        "converged_starts": len(trials), "selected_bins": int(left.sum() + right.sum()),
        "left_center_range_ps": [float(hist.centers[left].min()), float(hist.centers[left].max())],
        "right_center_range_ps": [float(hist.centers[right].min()), float(hist.centers[right].max())],
        "tail_log10_rmse": float(np.sqrt(np.mean(best.fun ** 2))),
        "pdf_tail_log10_rmse": float(np.sqrt(np.mean(pdf_error ** 2))),
        "minimum_tail_observations": None if hist.observations is None else low * hist.observations,
        "bounds": [[-8, .2, 0, .2, .01], [0, 5, 8, 5, .99]],
        "at_bound": bool(np.any(best.active_mask)),
    }
    return details


def metrics(p, targets=TARGETS):
    # Scan first, then refine every local minimum instead of assuming unimodality.
    grid = np.linspace(TI, TIP1, 4097)
    vals = log_ber(grid, p)
    local = np.flatnonzero((vals[1:-1] <= vals[:-2]) & (vals[1:-1] <= vals[2:])) + 1
    candidates = [(float(vals[0]), TI), (float(vals[-1]), TIP1)]
    for i in local:
        opt = minimize_scalar(lambda t: float(log_ber(t, p)), bounds=(grid[i-1], grid[i+1]),
                              method="bounded", options={"xatol": 1e-12})
        candidates.append((float(opt.fun), float(opt.x)))
    value, tmin = min(candidates)
    eyes = {}
    for target in targets:
        if not 0 < target < 1:
            raise ValueError("BER targets must be in (0,1)")
        if value > np.log(target):
            eyes[f"{target:.0e}"] = {"exists": False, "width_ps": 0.0, "width_ui": 0.0,
                                   "left_ps": None, "right_ps": None}
            continue
        # Include the refined optimum so even a sub-grid opening is detected.
        nodes = np.unique(np.r_[grid, tmin])
        f = log_ber(nodes, p) - np.log(target)
        roots = [float(nodes[i]) for i in np.flatnonzero(f == 0)]
        for i in np.flatnonzero(f[:-1] * f[1:] < 0):
            roots.append(float(brentq(lambda t: float(log_ber(t, p) - np.log(target)),
                                      nodes[i], nodes[i+1], xtol=1e-12)))
        cuts = sorted(set([TI, TIP1] + roots))
        intervals = [(a, b) for a, b in zip(cuts[:-1], cuts[1:])
                     if float(log_ber((a + b)/2, p)) <= np.log(target)]
        if not intervals:
            a = b = tmin  # A tangency has zero eye width.
        else:
            a, b = max(intervals, key=lambda ab: ab[1] - ab[0])
        eyes[f"{target:.0e}"] = {"exists": b > a, "width_ps": b-a,
            "width_ui": (b-a)/UI_PS, "left_ps": a, "right_ps": b}
    return {"minimum_ber": float(np.exp(value)), "optimal_time_ps": tmin, "eye_widths": eyes}


def load_primary(path, metadata):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    groups = {}
    for case in CASES:
        for side in ("tx", "rx"):
            selected = data[(data["case"] == case) & (data["side"] == side)]
            if len(selected) != 200:
                raise ValueError("Expected 200 bins per case and side")
            hist = Histogram(selected["time_seconds"] * 1e12, selected["pmf"],
                             metadata["cases"][case][side]["observation_count"])
            if not np.isclose(hist.dt, .125) or not np.allclose(hist.pdf, selected["pdf_per_second"] / 1e12):
                raise ValueError("Inconsistent primary PDF or physical time step")
            groups[f"{case}_{side}"] = hist
    return groups


def load_legacy(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.size != 200:
        raise ValueError("Expected 200 legacy bins")
    # Known display offset from the supplied HW12 code, not a noise-floor guess.
    centers = (np.arange(200) - 100) * (UI_PS / 200)
    groups, removed = {}, {}
    for i, case in enumerate(CASES, 1):
        for side, prefix in (("tx", "jtx"), ("rx", "jty")):
            raw = np.asarray(data[f"{prefix}{i}"])
            if np.any(raw < 0) or not np.all(np.isfinite(raw)):
                raise ValueError("Invalid legacy histogram")
            corrected = np.maximum(raw - 1e-15, 0)
            key = f"{case}_{side}"
            groups[key] = Histogram(centers, corrected)
            removed[key] = int(np.count_nonzero(corrected == 0))
    return groups, {"subtracted_display_offset": 1e-15, "zero_bins_after_correction": removed,
        "original_coordinate_spacing_ps": float(np.median(np.diff(data["tjp"]))),
        "corrected_coordinate_spacing_ps": .125,
        "coordinate_rule": "(index - 100) * 0.125 ps; physical HW12 column spacing",
        "observation_count": None, "provenance": "supplied export; exact generating run is unavailable"}


def positive(a):
    return np.ma.masked_less_equal(a, 0)


def save_figure(fig, output, name):
    fig.savefig(output / name, dpi=180, bbox_inches="tight", metadata={"Software": "HW13 analysis"})
    plt.close(fig)


def make_plots(groups, fits, sensitivity, output):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": .23})
    time = np.linspace(0, 50, 4001)
    for kind in ("cdf", "bathtub", "pdf"):
        fig, ax = plt.subplots(figsize=(8.1, 4.4), layout="constrained")
        for i, case in enumerate(CASES):
            for side in ("tx", "rx"):
                h = groups[f"{case}_{side}"]
                if kind == "cdf":
                    x, y = h.edges, h.cumulative
                elif kind == "bathtub":
                    x, y = time, h.ber(time)
                else:
                    x, y = h.centers, h.pdf
                ax.semilogy(x, positive(y), color=COLORS[i], ls="--" if side == "tx" else "-",
                            lw=1.4, marker="o" if kind == "pdf" and np.count_nonzero(y) == 1 else None,
                            markersize=4, label=f"Case {i+1} {side.upper()}: {LABELS[i]}")
        ax.set(xlabel="Sampling time (ps)" if kind == "bathtub" else "Jitter (ps)",
               ylabel={"cdf": "Cumulative probability", "bathtub": "BER", "pdf": "Probability density (1/ps)"}[kind],
               ylim=(1e-6, 10 if kind == "pdf" else 1),
               xlim=(0, 50) if kind == "bathtub" else (-12.5, 12.5),
               title=f"HW12-derived empirical {kind.upper() if kind != 'bathtub' else 'bathtub curves'}")
        ax.legend(fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(.5, -.16))
        save_figure(fig, output, f"primary_{kind}.png")
    for case in CASES[2:]:
        h = groups[f"{case}_rx"]
        rec = fits[case]
        p = rec["parameter_vector"]
        left, right = tail_masks(h, rec["tail_min"], rec["tail_max"])
        fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.1), layout="constrained")
        t = np.linspace(-12.5, 12.5, 2001)
        a.semilogy(h.centers, positive(h.pdf), color=COLORS[0], lw=1.2, label="Empirical PDF")
        a.semilogy(t, np.exp(log_pdf(t, p)), "--", color=COLORS[1], label="Two-Gaussian model")
        selected = (left | right) & (h.pmf > 0)
        a.scatter(h.centers[selected], h.pdf[selected], s=12, color=COLORS[0], label="Selected tail bins", zorder=3)
        a.set(xlim=(-12.5, 12.5), ylim=(1e-8, 1), xlabel="Jitter (ps)", ylabel="PDF (1/ps)", title="Tail fit; center excluded")
        b.semilogy(time, positive(h.ber(time)), color=COLORS[0], label="Empirical BER")
        b.semilogy(time, model_ber(time, p), "--", color=COLORS[1], label="Extrapolated model")
        b.axhline(1e-12, color=".45", lw=.8, ls=":", label="Target 1e-12")
        m = rec["metrics"]
        b.scatter([m["optimal_time_ps"]], [m["minimum_ber"]], s=26, marker="x", color=COLORS[1])
        b.set(xlim=(10, 40), ylim=(1e-16, 1), xlabel="Sampling time (ps)", ylabel="BER",
              title=f"Model minimum: {m['minimum_ber']:.3g}")
        for ax in (a, b):
            ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.17))
        save_figure(fig, output, f"primary_{case.split('_')[0]}_fit.png")
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4), layout="constrained")
    for ax, case in zip(axes, CASES[2:]):
        for j, high in enumerate((.01, .03, .1)):
            rows = [r for r in sensitivity if r["case"] == case and r["tail_max"] == high]
            ax.loglog([r["tail_min"] for r in rows], [r["metrics"]["minimum_ber"] for r in rows],
                      marker="o", color=COLORS[j], label=f"Upper tail probability {high:g}")
        ax.axhline(1e-12, color=".4", ls=":", lw=1)
        ax.set(xlabel="Lower tail probability included in fit", ylabel="Predicted minimum BER",
               title=f"{case.split('_')[0].title()} Rx: window sensitivity")
        ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.18))
    save_figure(fig, output, "sensitivity.png")


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    meta = json.loads((BASE / "inputs/hw12_source_summary.json").read_text(encoding="utf-8"))
    primary_path = BASE / "inputs/hw12_verified_jitter.csv"
    groups = load_primary(primary_path, meta)
    legacy, preprocessing = load_legacy(BASE / "inputs/jitter_variables_for_hw13.csv")
    fits, legacy_fits, sensitivity = {}, {}, []
    for case in CASES[2:]:
        fit = fit_tails(groups[f"{case}_rx"])
        fit["metrics"] = metrics(fit["parameter_vector"])
        fits[case] = fit
        old = fit_tails(legacy[f"{case}_rx"])
        old["metrics"] = metrics(old["parameter_vector"])
        legacy_fits[case] = old
        for low in (1e-4, 3e-4, 1e-3):
            for high in (.01, .03, .1):
                rec = fit if (low, high) == (TAIL_MIN, TAIL_MAX) else fit_tails(groups[f"{case}_rx"], low, high)
                sensitivity.append({"case": case, "tail_min": low, "tail_max": high,
                    "parameter_vector": rec["parameter_vector"], "tail_log10_rmse": rec["tail_log10_rmse"],
                    "metrics": metrics(rec["parameter_vector"])})
    sample_stats = {}
    for key, h in groups.items():
        sample_stats[key] = {"observations": h.observations, "pmf_sum": float(h.pmf.sum()),
            "pdf_integral": float(h.pdf.sum() * h.dt), "dt_ps": h.dt,
            "zero_bins": int(np.count_nonzero(h.pmf == 0)),
            "smallest_nonzero_pmf": float(h.pmf[h.pmf > 0].min())}
    historical = {"archive_case3": [-1.2, 1.8, 1.2, 1.5, .5], "archive_case4": [-2.7, 2, 2.1, 1.8, .5],
                  "standalone_case3": [-1.3, 1.7, 1.4, 1.3, .6], "standalone_case4": [-2.8, 2.1, 2, 1.7, .4]}
    summary = {"schema_version": 1, "parameters": {"ui_ps": UI_PS, "dt_ps": .125,
        "rho_transition": RHO, "ideal_boundaries_ps": [TI, TIP1], "tail_min": TAIL_MIN, "tail_max": TAIL_MAX,
        "fit_objective": "unweighted log10 left-CDF and right-SF residuals on selected tails",
        "fit_uncertainty": "window sensitivity is not a confidence interval",
        "histogram_interpolation": "uniform density within each centered 0.125 ps cell"},
        "source": {"data": "inputs/hw12_verified_jitter.csv", "sha256": hashlib.sha256(primary_path.read_bytes()).hexdigest(),
                   "bits_per_case": meta["parameters"]["bits"], "description": "reused verified HW12 baseline; channel simulation not rerun by HW13"},
        "sample_statistics": sample_stats, "cases": fits, "sensitivity": sensitivity,
        "legacy": {"preprocessing": preprocessing, "refitted_cases": legacy_fits,
                   "supplied_parameter_predictions": {key: {"parameter_vector": p, "metrics": metrics(p)} for key, p in historical.items()}},
        "limitations": ["Low-BER predictions rely on Gaussian tails and a two-component model.",
            "CDF tail points are correlated; fit RMSE is descriptive, not a statistical uncertainty estimate.",
            "Eye-slice occupancy weights slopes differently from equal-weight crossing-event histograms.",
            "Zero empirical errors do not establish zero physical BER.",
            "The one-UI histogram folds extreme jitter; no low-BER measurement or ADS validation was performed."]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_csv(output / "empirical_distributions.csv", ["dataset", "group", "center_ps", "pmf", "pdf_per_ps", "cdf_at_upper_edge", "sf_at_lower_edge"],
        ((dataset, key, x, mass, density, cdf, sf) for dataset, allgroups in (("primary", groups), ("legacy_corrected", legacy))
         for key, h in allgroups.items() for x, mass, density, cdf, sf in zip(h.centers, h.pmf, h.pdf, h.cumulative[1:], h.survival[:-1])))
    time = np.linspace(0, 50, 4001)
    keys = list(groups)
    curve_columns = [groups[k].ber(time) for k in keys] + [model_ber(time, fits[c]["parameter_vector"]) for c in CASES[2:]]
    write_csv(output / "bathtub_curves.csv", ["time_ps"] + keys + [c + "_rx_model" for c in CASES[2:]], zip(time, *curve_columns))
    make_plots(groups, fits, sensitivity, output)
    for case, rec in fits.items():
        m = rec["metrics"]
        print(f"{case}: min BER={m['minimum_ber']:.8g}; EW@1e-12={m['eye_widths']['1e-12']['width_ps']:.6f} ps")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE / "results")
    arguments = parser.parse_args()
    run(arguments.output.resolve())
