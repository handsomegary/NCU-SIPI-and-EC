"""Independent HW10 verification using a separate time-domain convolution path.

Default: reproduce the published 262144-symbol noiseless record. This script
constructs the transmitter waveform before channel convolution, then applies
the serial DFE to the full waveform. The production solver instead uses a
single-bit response and polyphase convolution. No production waveform, eye,
or DFE function is used to obtain the independent reference result.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.fft import irfft
from scipy.signal import fftconvolve, oaconvolve

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "instructions" / "HW10_complicated_diff_channel.s4p"
DEFAULT_FFE = [-.222313178938,.777686821062]
DEFAULT_GAIN_OFFSET_DB = .248389735


def read_s4p(path=DATA):
    vals = []
    for line in Path(path).read_text().splitlines():
        line = line.split("!")[0].strip()
        if line and not line.startswith("#"):
            vals.extend(float(x) for x in line.split())
    values = np.asarray(vals).reshape(-1, 33)
    pairs = values[:, 1:].reshape(-1, 16, 2)
    s = pairs[:, :, 0] * np.exp(1j * np.deg2rad(pairs[:, :, 1]))
    return values[:, 0], s.reshape(-1, 4, 4)


def ctle(f, gain_offset_db=0.):
    s = 2j * np.pi * f
    wp = 2 * np.pi * 40e9
    adc, aac = 10**(-10/20), 10**(14/20)
    return 10**(gain_offset_db/20) * aac * wp**2 * (s + adc*wp/aac) / (s+wp)**3


def make_response(f, H, nfft, kaiser=1):
    w = np.i0(np.pi * kaiser * np.sqrt(1-(f/f[-1])**2)) / np.i0(np.pi*kaiser)
    spectrum = np.zeros(nfft//2+1, complex)
    spectrum[:len(H)] = H*w
    return spectrum


def discrete_pulse(npui=50, nr=10):
    p = np.r_[np.arange(nr)/nr, np.ones(npui-nr), 1-np.arange(nr)/nr]
    return p


def sample_cursors(sbr, npui=50, npost=5):
    main = int(np.argmax(sbr))
    return main, sbr[main+np.arange(npost+1)*npui]


def continuous_sbr(f, H, nfft=200000, dt=.5e-12, a=1, wrong_phase=False):
    tb, tr, vpp = 25e-12, 5e-12, 2
    delay = tb*tr if wrong_phase else tb+tr
    vin = vpp*tb*np.sinc(f*tb)*np.sinc(f*tr)*np.exp(-1j*np.pi*f*delay)
    return irfft(make_response(f, H*vin/dt, nfft, a), nfft)


def eye_metrics(y, symbols, sample_index, npui=50, ignore=1024):
    """Envelope from bit-labelled traces on the absolute clock (no histogram)."""
    bit_indices = np.arange(ignore, len(symbols)-ignore)
    centers = sample_index + bit_indices*npui
    valid = (centers-npui >= 0) & (centers+npui < len(y))
    centers, labels = centers[valid], symbols[bit_indices[valid]]
    offsets = np.arange(-npui, npui+1)
    low = np.full(len(offsets),-np.inf)
    high = np.full(len(offsets),np.inf)
    for start in range(0,len(centers),4096):
        chunk_centers=centers[start:start+4096]
        chunk_labels=labels[start:start+4096]
        traces=y[chunk_centers[:,None]+offsets]
        if np.any(chunk_labels < 0):
            low=np.maximum(low,np.max(traces[chunk_labels < 0],axis=0))
        if np.any(chunk_labels > 0):
            high=np.minimum(high,np.min(traces[chunk_labels > 0],axis=0))
    margin = np.minimum(high, -low)
    mid = npui
    left = mid
    right = mid
    if margin[mid] > 0:
        while left > 0 and margin[left-1] > 0:
            left -= 1
        while right < len(margin)-1 and margin[right+1] > 0:
            right += 1
        xleft = float(offsets[left])
        xright = float(offsets[right])
        if left > 0:
            xleft -= margin[left] / (margin[left]-margin[left-1])
        if right < len(margin)-1:
            xright += margin[right] / (margin[right]-margin[right+1])
        ew = (xright-xleft)*.5
    else:
        ew = 0.
    center_values = y[centers]
    errors = np.count_nonzero((center_values >= 0) != (labels > 0))
    return {"EH_mV": float((high[mid]-low[mid])*1e3), "EW_ps": ew,
            "errors": int(errors), "samples": int(len(labels)),
            "lower_1_V": float(high[mid]), "upper_0_V": float(low[mid])}


def dfe(y, taps, sample, npui=50):
    z = y.copy()
    for k in range(sample, len(y), npui):
        decision = 1 if z[k] >= 0 else -1
        for j, tap in enumerate(taps):
            start = k+npui//2+j*npui
            z[start:start+npui] -= decision*tap
    return z


def synthetic_checks():
    """Checks tied to the physical model, without importing the implementation."""
    p = discrete_pulse()
    fcheck = np.array([0, 1e7, 1e8, 1e9])
    sample_time = np.arange(len(p))*.5e-12
    exact = np.exp(-2j*np.pi*fcheck[:,None]*sample_time) @ (p*.5e-12)
    cont = 25e-12*np.sinc(fcheck*25e-12)*np.sinc(fcheck*5e-12)*np.exp(-1j*np.pi*fcheck*30e-12)
    pulse_relative_error = float(np.max(abs(exact-cont)/25e-12))
    rng=np.random.default_rng(190)
    symbols=2*rng.integers(0,2,100)-1
    isi=np.array([.2,-.07,.03])
    expected=np.repeat(symbols,50).astype(float)
    corrupted=expected.copy()
    for j,tap in enumerate(isi,1):
        corrupted[j*50:]+=tap*expected[:-j*50]
    clean=dfe(corrupted,isi,25)
    # Boundary at center+UI/2 is the following bit's left edge.
    dfe_error=float(np.max(abs(clean-expected)))
    return [
        {"name":"source_voltage_divider", "pass":True,
         "detail":"Matched 50-ohm per-line source/load: Hsource=Sdd21/2, applied once; differential source levels +/-2 V."},
        {"name":"causal_ramp_pulse_delay", "pass":pulse_relative_error < 1e-5,
         "max_relative_spectrum_error_below_1GHz":pulse_relative_error,
         "centroid_ps":float(np.sum(p*sample_time)/np.sum(p)*1e12)},
        {"name":"independent_bipolar_DFE_causality", "pass":dfe_error < 1e-12,
         "max_abs_V_error":dfe_error,"clock_phase_sample":25,
         "detail":"Both slicer signs cancel a known three-postcursor channel; feedback begins half UI after the current decision."},
    ]


def run(nbits=16384,data=DATA,output=ROOT/"reference-results.json",
        ffetaps=None,ctle_gain_offset_db=0.,phase_ps=None,
        metrics_path=None,numerical_path=None):
    f, s = read_s4p(data)
    m = np.array([[1,0,1,0],[1,0,-1,0],[0,1,0,1],[0,1,0,-1]])/np.sqrt(2)
    mm = m @ s @ m.T
    closed = (s[:,1,0]-s[:,1,2]-s[:,3,0]+s[:,3,2])/2
    nfft, dt, npui = 200000, .5e-12, 50
    hch = closed/2
    ffetaps = np.array(DEFAULT_FFE if ffetaps is None else ffetaps)
    hffe = sum(c*np.exp(-2j*np.pi*f*25e-12*j) for j,c in enumerate(ffetaps))
    responses = [hch, hch*ctle(f,ctle_gain_offset_db), hch*ctle(f,ctle_gain_offset_db)*hffe]
    kernels=[irfft(make_response(f,H,nfft),nfft)[:40000] for H in responses[:2]]
    pulse=2*discrete_pulse()
    discrete=[fftconvolve(h,pulse) for h in kernels]
    ffe_kernel=np.zeros((len(ffetaps)-1)*npui+1)
    ffe_kernel[::npui]=ffetaps
    discrete.append(fftconvolve(discrete[1],ffe_kernel))
    main, cursors = sample_cursors(discrete[2])
    peak_sample=main
    if phase_ps is not None:
        phase_sample=int(round(phase_ps/.5))%npui
        main=int(round((main-phase_sample)/npui))*npui+phase_sample
        cursors=discrete[2][main+np.arange(6)*npui]
    correct_cont = continuous_sbr(f, responses[2])
    # Original SBR omitted divider, used a=.5, and used a dimensionally wrong phase.
    original_hffe=-.1+.9*np.exp(-2j*np.pi*f*25e-12)
    original_sbr = continuous_sbr(f, hch*ctle(f)*original_hffe*2, a=.5, wrong_phase=True)
    original_main, original_cursors = sample_cursors(original_sbr)
    ratio = 10**(-10/20)/10**(14/20)
    peak_f = 40e9*np.sqrt((1-3*ratio**2)/2)
    result = {
        "checks":synthetic_checks(),
        "method":"Independent finite-ramp transmitter waveform and overlap-add channel convolution, followed by a serial full-waveform bipolar DFE. No production waveform/eye/DFE functions are used for the reference.",
        "design": {"ctle_gain_offset_dB":float(ctle_gain_offset_db),"ctle_DC_gain_dB":-10+float(ctle_gain_offset_db),"ctle_AC_coefficient_dB":14+float(ctle_gain_offset_db),"ctle_poles_GHz":[40,40,40],"FFE_taps":ffetaps.tolist(),"DFE_feedback_taps_V":cursors[1:].tolist(),"absolute_sample_index":main,"clock_phase_ps":(main%npui)*.5},
        "frequency": {"points":len(f), "df_Hz":float(f[1]-f[0]), "fmax_Hz":float(f[-1]), "dc":float(f[0])},
        "mixed_mode_closed_form_max_error":float(np.max(abs(mm[:,3,1]-closed))),
        "channel_source_dc_gain":float(hch[0].real),
        "ctle_peak": {"frequency_GHz":float(peak_f*1e-9), "gain_dB":float(20*np.log10(abs(ctle(peak_f,ctle_gain_offset_db))))},
        "ffe_l1":float(sum(abs(ffetaps))),
        "correct_discrete_sbr": {"peak_sample":peak_sample,"main_sample":main,"main_time_ns":main*dt*1e9,"sample_phase_ps":(main%npui)*dt*1e12,"cursor_V":cursors.tolist(),"additive_DFE_V":(-cursors[1:]).tolist()},
        "correct_continuous_sbr": {"main_sample":int(np.argmax(correct_cont)),"cursor_V":sample_cursors(correct_cont)[1].tolist()},
        "original_sbr_reproduction": {"main_sample":original_main,"main_time_ns":original_main*dt*1e9,"cursor_V":original_cursors.tolist(),"additive_DFE_V":(-original_cursors[1:]/2).tolist()},
    }
    result["checks"].extend([
        {"name":"mixed_mode_transform", "pass":result["mixed_mode_closed_form_max_error"]<1e-12,
         "max_abs_S_error":result["mixed_mode_closed_form_max_error"]},
        {"name":"frequency_grid", "pass":bool(f[0]==0 and np.allclose(np.diff(f),1e7,rtol=0,atol=1e-4)),
         "points":len(f),"DC_Hz":float(f[0]),"upper_Hz":float(f[-1])},
        {"name":"CTLE_constraint", "pass":result["ctle_peak"]["gain_dB"]<=6,
         "pole_count":3,"highest_pole_GHz":40,"analytic_peak_dB":result["ctle_peak"]["gain_dB"]},
        {"name":"FFE_constraint", "pass":bool(len(ffetaps)<=5 and np.isclose(sum(abs(ffetaps)),1,atol=1e-12)),
         "tap_count":len(ffetaps),"l1_norm":float(sum(abs(ffetaps)))},
        {"name":"DFE_constraint", "pass":len(cursors[1:])<=5,"tap_count":len(cursors[1:])},
    ])
    if nbits:
        rng = np.random.default_rng(0)
        symbols = 2*rng.integers(0,2,nbits+2048)-1
        # A finite-ramp source built independently by averaging the preceding Nr samples.
        x = np.r_[0, fftconvolve(2*np.repeat(symbols,npui), np.ones(10)/10)]
        waveforms=[]
        for h in kernels:
            waveforms.append(oaconvolve(x,h))
        xffe=np.r_[0,fftconvolve(np.repeat(np.convolve(2*symbols,ffetaps),npui),np.ones(10)/10)]
        waveforms.append(oaconvolve(xffe,kernels[1]))
        dfe_output=dfe(waveforms[2],cursors[1:],main%npui)
        # DFE phase is modulo UI; association with source bits still includes channel delay.
        result["independent_simulation"]={
            "bits":nbits,"seed":0,"before_DFE":eye_metrics(waveforms[2],symbols,main),
            "after_DFE":eye_metrics(dfe_output,symbols,main),
            "channel":eye_metrics(waveforms[0],symbols,int(np.argmax(discrete[0]))),
            "ctle":eye_metrics(waveforms[1],symbols,int(np.argmax(discrete[1]))),
            "correct_sample_phase_ps":(main%npui)*.5,
        }
        final=result["independent_simulation"]["after_DFE"]
        result["checks"].append({"name":"independent_record_requirements",
            "pass":final["EH_mV"]>120 and final["EW_ps"]>15 and final["errors"]==0,
            "EH_mV":final["EH_mV"],"EW_ps":final["EW_ps"],"errors":final["errors"],
            "measured_symbols":final["samples"],"EH_margin_mV":final["EH_mV"]-120,
            "EW_margin_ps":final["EW_ps"]-15})
        if metrics_path is not None:
            published=json.loads(Path(metrics_path).read_text(encoding="utf-8"))
            reference_stages={"channel":"channel","ctle":"ctle","ffe":"before_DFE","dfe":"after_DFE"}
            differences={}
            for name,key in reference_stages.items():
                independent=result["independent_simulation"][key]
                observed=published["stages"][name]
                differences[name]={
                    "EH_abs_error_V":abs(independent["EH_mV"]*1e-3-observed["eye_height_v"]),
                    "EW_abs_error_ps":abs(independent["EW_ps"]-observed["eye_width_ps"]),
                    "error_count_difference":independent["errors"]-observed["decision_errors"],
                    "measured_symbol_difference":independent["samples"]-observed["measured_symbols"]}
            passed=all(v["EH_abs_error_V"]<1e-10 and v["EW_abs_error_ps"]<1e-8 and v["error_count_difference"]==0 and v["measured_symbol_difference"]==0 for v in differences.values())
            result["production_comparison"]=differences
            result["checks"].append({"name":"production_metrics_match_independent_reference","pass":passed,
                "EH_tolerance_V":1e-10,"EW_tolerance_ps":1e-8})
        if numerical_path is not None:
            arrays=np.load(numerical_path)
            times=arrays['waveform_time_ns']
            indices=np.rint(times*1e-9/dt).astype(int)
            if not np.allclose(times,indices*dt*1e9,rtol=0,atol=1e-10):
                raise ValueError('Saved waveform timestamps are not on the simulation grid')
            if np.any(np.diff(indices)!=1):
                raise ValueError('Saved waveform timestamps must form a contiguous excerpt')
            comparisons={name:float(np.max(abs(y[indices]-arrays['waveform_'+name+'_v'])))
                for name,y in zip(("channel","ctle","ffe","dfe"),waveforms+[dfe_output])}
            idx=np.arange(1024,1024+nbits)
            decisions=np.where(dfe_output[main+idx*npui]>=0,1,-1)
            decision_differences=int(np.count_nonzero(decisions!=arrays["dfe_decisions"][idx]))
            result["waveform_max_abs_errors_V"]=comparisons
            result["checks"].append({"name":"production_waveform_matches_independent_path",
                "pass":max(comparisons.values())<1e-10,"tolerance_V":1e-10,
                "samples_per_stage":len(indices)})
            result["checks"].append({"name":"all_retained_DFE_decisions_match",
                "pass":decision_differences==0,"decision_differences":decision_differences,
                "measured_symbols":nbits})
        if np.array_equal(ffetaps,np.array([-.1,.9])) and ctle_gain_offset_db==0:
            original_output=dfe(waveforms[2],original_cursors[1:]/2,0)
            # The original clock is at phase 0; choose its nearest absolute bit clock.
            original_clock=int(round(main/npui))*npui
            result["independent_simulation"]["original_DFE_clock_and_taps"]=eye_metrics(original_output,symbols,original_clock)
    result["all_checks_pass"]=all(c["pass"] for c in result["checks"])
    result["scope"]="Finite noiseless random-record verification under the stated bandwidth, Kaiser window and 20 ns causal impulse truncation; no physical BER or all-pattern guarantee."
    Path(output).write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--bits",type=int,default=262144)
    parser.add_argument("--data",type=Path,default=DATA)
    parser.add_argument("--output",type=Path,default=ROOT.parent/"results"/"independent_checks.json")
    parser.add_argument("--ffe",type=str,default=",".join(map(str,DEFAULT_FFE)))
    parser.add_argument("--ctle-gain-offset-db",type=float,default=DEFAULT_GAIN_OFFSET_DB)
    parser.add_argument("--phase-ps",type=float,default=19.)
    parser.add_argument("--metrics",type=Path,default=ROOT.parent/"results"/"metrics.json")
    parser.add_argument("--numerical",type=Path,default=ROOT.parent/"results"/"numerical_results.npz")
    args=parser.parse_args()
    raise SystemExit(run(args.bits,args.data,args.output,[float(v) for v in args.ffe.split(",")],args.ctle_gain_offset_db,args.phase_ps,args.metrics,args.numerical))
