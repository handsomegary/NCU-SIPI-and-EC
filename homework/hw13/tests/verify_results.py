"""Verify saved outputs independently. Does not import production analysis code.

Gaussian probabilities use Python math.erfc; histograms are integrated directly
as cell masses; minima and roots use independent scalar search routines.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

BASE = Path(__file__).resolve().parents[1]
RESULTS = BASE / "results"
SQRT2 = math.sqrt(2)
CASES = ("case1_no_tx_jitter", "case2_sj", "case3_rj", "case4_sj_rj")
CHECKS = []


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def close(actual, expected, *, relative=1e-9, absolute=1e-13):
    return math.isfinite(actual) and math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute)


def rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def tail(x, p, left=True):
    m1, s1, m2, s2, w = p
    sign = -1 if left else 1
    return w * .5 * math.erfc(sign * (x-m1)/(s1*SQRT2)) + (1-w) * .5 * math.erfc(sign * (x-m2)/(s2*SQRT2))


def density(x, p):
    m1, s1, m2, s2, w = p
    return (w/s1*math.exp(-.5*((x-m1)/s1)**2) + (1-w)/s2*math.exp(-.5*((x-m2)/s2)**2)) / math.sqrt(2*math.pi)


def ber(t, p):
    return .5 * (tail(t-12.5, p, False) + tail(t-37.5, p))


def golden_min(fn, a, b):
    ratio = (math.sqrt(5)-1)/2
    c, d = b-ratio*(b-a), a+ratio*(b-a)
    fc, fd = fn(c), fn(d)
    for _ in range(100):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b-ratio*(b-a); fc = fn(c)
        else:
            a, c, fc = c, d, fd
            d = a+ratio*(b-a); fd = fn(d)
    t = (a+b)/2
    return t, fn(t)


def independent_minimum(p):
    grid = [12.5 + 25*i/1000 for i in range(1001)]
    vals = [ber(t, p) for t in grid]
    candidates = [(grid[0], vals[0]), (grid[-1], vals[-1])]
    for i in range(1, 1000):
        if vals[i] <= vals[i-1] and vals[i] <= vals[i+1]:
            candidates.append(golden_min(lambda t: ber(t, p), grid[i-1], grid[i+1]))
    return min(candidates, key=lambda pair: pair[1])


def bisect(fn, a, b):
    fa = fn(a)
    for _ in range(100):
        c = (a+b)/2
        if fn(c)*fa > 0:
            a, fa = c, fn(c)
        else:
            b = c
    return (a+b)/2


def empirical_tail(t, hist, left=True):
    # Integrate each uniform-density histogram cell directly, without CDF interpolation.
    if left:
        return math.fsum(m*min(1., max(0., (t-(x-.0625))/.125)) for x, m in hist)
    return math.fsum(m*min(1., max(0., ((x+.0625)-t)/.125)) for x, m in hist)


def empirical_ber(t, hist):
    return .5*(empirical_tail(t-12.5, hist, False)+empirical_tail(t-37.5, hist))


def verify_model(label, record):
    p, saved = record["parameter_vector"], record["metrics"]
    require(len(p) == 5 and all(math.isfinite(v) for v in p) and p[1] > 0 and p[3] > 0 and 0 <= p[4] <= 1, label+": valid parameters")
    tmin, minimum = independent_minimum(p)
    require(close(saved["minimum_ber"], minimum, relative=3e-9, absolute=0), label+": saved minimum")
    require(abs(saved["optimal_time_ps"]-tmin) < 3e-5, label+": optimal time")
    require(close(ber(saved["optimal_time_ps"], p), saved["minimum_ber"], relative=3e-10, absolute=0), label+": BER at optimum")
    # dBER/dt = rho*(PDF(t-T[i+1]) - PDF(t-T[i])).
    a, b = density(saved["optimal_time_ps"]-37.5, p), density(saved["optimal_time_ps"]-12.5, p)
    stationarity = abs(a-b)/max(a, b, 1e-300)
    require(stationarity < 2e-5, label+": minimum stationarity")
    for target_text, eye in saved["eye_widths"].items():
        target = float(target_text)
        expected_open = minimum < target
        require(eye["exists"] == expected_open, label+": opening existence at "+target_text)
        if not expected_open:
            require(eye["width_ps"] == 0 and eye["width_ui"] == 0 and eye["left_ps"] is None and eye["right_ps"] is None, label+": explicit zero-width closed eye")
            continue
        left = bisect(lambda t: ber(t,p)-target, 12.5, tmin)
        right = bisect(lambda t: ber(t,p)-target, tmin, 37.5)
        require(close(eye["left_ps"], left, absolute=2e-8), label+": left root")
        require(close(eye["right_ps"], right, absolute=2e-8), label+": right root")
        require(close(eye["width_ps"], right-left, absolute=2e-8), label+": eye width")
        require(close(eye["width_ui"], (right-left)/25, absolute=1e-9), label+": UI width")
        for root in (eye["left_ps"], eye["right_ps"]):
            require(close(ber(root,p), target, relative=2e-8, absolute=0), label+": root residual")
    return {"model": label, "minimum_ber_erfc": minimum, "optimal_time_ps_erfc": tmin,
            "relative_stationarity_residual": stationarity, "targets_checked": len(saved["eye_widths"])}


def verify_fit(label, rec, hist):
    p = rec["parameter_vector"]
    low, high = rec.get("tail_min", 3e-4), rec.get("tail_max", .03)
    errors, pdf_errors, selected_left, selected_right = [], [], [], []
    for i, (x, mass) in enumerate(hist):
        c, s = math.fsum(m for _,m in hist[:i+1]), math.fsum(m for _,m in hist[i:])
        left, right = low <= c <= high, low <= s <= high
        if left:
            errors.append(math.log10(tail(x+.0625,p)/c)); selected_left.append(x)
        if right:
            errors.append(math.log10(tail(x-.0625,p,False)/s)); selected_right.append(x)
        if (left or right) and mass > 0:
            pdf_errors.append(math.log10(density(x,p)/(mass/.125)))
    rmse = math.sqrt(math.fsum(e*e for e in errors)/len(errors))
    require(close(rec["tail_log10_rmse"], rmse, relative=1e-7, absolute=2e-10), label+": fit objective")
    if "selected_bins" in rec:
        require(rec["selected_bins"] == len(errors), label+": selected bins")
        require(close(rec["pdf_tail_log10_rmse"], math.sqrt(math.fsum(e*e for e in pdf_errors)/len(pdf_errors)), absolute=2e-10), label+": PDF tail residual")
        require(rec["left_center_range_ps"] == [min(selected_left),max(selected_left)] and rec["right_center_range_ps"] == [min(selected_right),max(selected_right)], label+": tail selection ranges")


def verify_all():
    summary = json.loads((RESULTS/"summary.json").read_text(encoding="utf-8"))
    input_path = BASE/"inputs/hw12_verified_jitter.csv"
    source_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
    require(summary["source"]["sha256"] == source_hash, "primary source hash")
    upstream = BASE.parent/"hw12/results/jitter_distributions.csv"
    require(input_path.read_bytes() == upstream.read_bytes(), "primary bytes equal upstream HW12 export")
    metadata = json.loads((BASE/"inputs/hw12_source_summary.json").read_text(encoding="utf-8"))
    require(summary["source"]["bits_per_case"] == metadata["parameters"]["bits"] == 2**20, "full HW12 bit count")
    groups, raw_pdfs = {}, {}
    for row in rows(input_path):
        key = row["case"]+"_"+row["side"]
        groups.setdefault(key,[]).append((float(row["time_seconds"])*1e12,float(row["pmf"])))
        raw_pdfs.setdefault(key,[]).append(float(row["pdf_per_second"])/1e12)
    require(len(groups) == 8, "all eight primary groups")
    for key, hist in groups.items():
        require(len(hist) == 200 and all(close(hist[i+1][0]-hist[i][0],.125,absolute=1e-12) for i in range(199)), key+": physical grid")
        require(all(m >= 0 and math.isfinite(m) for _,m in hist) and close(math.fsum(m for _,m in hist),1,absolute=2e-14), key+": normalized PMF")
        require(all(close(m/.125,pdf,absolute=2e-13) for (_,m),pdf in zip(hist,raw_pdfs[key])), key+": PDF time units")
        stats = summary["sample_statistics"][key]
        case, side = key.rsplit("_",1)
        require(stats["observations"] == metadata["cases"][case][side]["observation_count"], key+": observation count")
        require(stats["zero_bins"] == sum(m==0 for _,m in hist), key+": zeros preserved")
        require(close(stats["pmf_sum"],1) and close(stats["pdf_integral"],1) and close(stats["dt_ps"],.125), key+": saved integral statistics")
    CHECKS.append({"check":"primary provenance and PMF/PDF units", "groups":8, "sha256":source_hash})

    legacy = {}
    raw_legacy = rows(BASE/"inputs/jitter_variables_for_hw13.csv")
    for index, case in enumerate(CASES,1):
        for side,prefix in (("tx","jtx"),("rx","jty")):
            masses = [max(0.,float(r[f"{prefix}{index}"])-1e-15) for r in raw_legacy]
            total = math.fsum(masses)
            key = f"{case}_{side}"
            legacy[key] = [((i-100)*.125,m/total) for i,m in enumerate(masses)]
            require(summary["legacy"]["preprocessing"]["zero_bins_after_correction"][key] == sum(m==0 for m in masses), key+": legacy offset zero count")
    exported = {}
    for row in rows(RESULTS/"empirical_distributions.csv"):
        exported.setdefault((row["dataset"],row["group"]),[]).append(row)
    require(len(exported) == 16, "all empirical export groups")
    for dataset,hists in (("primary",groups),("legacy_corrected",legacy)):
        for key,hist in hists.items():
            out = exported[(dataset,key)]
            require(len(out) == 200,key+": exported length")
            for i,((x,mass),row) in enumerate(zip(hist,out)):
                require(close(float(row["center_ps"]),x,absolute=1e-12),key+": exported time")
                require(close(float(row["pmf"]),mass,relative=2e-12,absolute=1e-18),key+": exported PMF")
                require(close(float(row["pdf_per_ps"]),mass/.125,relative=2e-12,absolute=1e-17),key+": exported PDF")
                require(close(float(row["cdf_at_upper_edge"]),math.fsum(m for _,m in hist[:i+1]),absolute=3e-14),key+": exported CDF")
                require(close(float(row["sf_at_lower_edge"]),math.fsum(m for _,m in hist[i:]),absolute=3e-14),key+": exported SF")
    CHECKS.append({"check":"all primary and corrected legacy empirical exports", "groups":16,"bins":3200})

    model_checks = []
    for case,rec in summary["cases"].items():
        model_checks.append(verify_model("primary/"+case,rec)); verify_fit(case,rec,groups[case+"_rx"])
    for i,rec in enumerate(summary["sensitivity"]):
        label=f"sensitivity/{i}/{rec['case']}"
        model_checks.append(verify_model(label,rec)); verify_fit(label,rec,groups[rec["case"]+"_rx"])
    for case,rec in summary["legacy"]["refitted_cases"].items():
        model_checks.append(verify_model("legacy_refit/"+case,rec)); verify_fit(case,rec,legacy[case+"_rx"])
    for label,rec in summary["legacy"]["supplied_parameter_predictions"].items():
        model_checks.append(verify_model("historical/"+label,rec))
    CHECKS.append({"check":"all saved model minima, stationarity, roots and widths with math.erfc", "models":len(model_checks),"details":model_checks})

    # Broad regression bounds come from the independent audit, not parameter freezing.
    c3,c4 = summary["cases"]["case3_rj"],summary["cases"]["case4_sj_rj"]
    require(close(c3["metrics"]["minimum_ber"],2.94553e-13,relative=.02,absolute=0), "independent Case 3 baseline")
    require(abs(c3["metrics"]["eye_widths"]["1e-12"]["width_ps"]-.83002) < .04, "independent Case 3 eye width baseline")
    require(close(c4["metrics"]["minimum_ber"],2.30408e-8,relative=.02,absolute=0), "independent Case 4 baseline")
    require(not c4["metrics"]["eye_widths"]["1e-12"]["exists"], "Case 4 no 1e-12 opening")
    CHECKS.append({"check":"independent audit baselines and explicit closed-eye semantics"})

    curves = rows(RESULTS/"bathtub_curves.csv")
    require(len(curves) == 4001, "saved curve sample count")
    max_empirical_error,max_model_relative_error = 0.,0.
    for row in curves:
        t = float(row["time_ps"])
        for key,hist in groups.items():
            expected,actual = empirical_ber(t,hist),float(row[key])
            max_empirical_error=max(max_empirical_error,abs(expected-actual))
            require(close(actual,expected,relative=2e-12,absolute=3e-14),key+": saved sampled bathtub")
        for case,rec in summary["cases"].items():
            expected,actual = ber(t,rec["parameter_vector"]),float(row[case+"_rx_model"])
            require(close(actual,expected,relative=3e-12,absolute=1e-300),case+": saved model curve")
            max_model_relative_error=max(max_model_relative_error,abs(actual-expected)/expected)
    CHECKS.append({"check":"every saved bathtub ordinate", "sampled_curves":8,"model_curves":2,"time_points":4001,
                   "maximum_sampled_absolute_error":max_empirical_error,"maximum_model_relative_error":max_model_relative_error})
    return {"passed":True,"python":platform.python_version(),"independent_oracles":["math.erfc Gaussian tails","analytic uniform-cell mass integration","golden-section minima and bisection roots"],
            "production_module_imported":False,"checks":CHECKS}


if __name__ == "__main__":
    try:
        result = verify_all()
    except Exception as error:
        result={"passed":False,"python":platform.python_version(),"error":str(error),"checks_completed":CHECKS}
        (RESULTS/"verification.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
        raise
    (RESULTS/"verification.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"passed":True,"check_groups":len(CHECKS),"saved_models_checked":26,"curve_ordinates_checked":40010},indent=2))
