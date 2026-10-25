"""Reproducible matched differential-channel simulation for HW10.

Voltages are measured after the matched source divider: H = Sdd21 / 2.
The source differential levels are +/-2 V. A sampled causal linear ramp,
the channel and CTLE are applied once, followed by transmitter FFE and a
serial decision-feedback equalizer. No transmitted bits enter the DFE.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.fft import irfft, rfft, next_fast_len
from scipy.signal import fftconvolve
from scipy.optimize import minimize_scalar

BIT_RATE = 40e9
UI = 1 / BIT_RATE
RISE = 5e-12
VPP = 2.0
ALPHA = 1.0
TRUNCATION = 20e-9
DEFAULT_FFE = np.array([-0.222313178938, 0.777686821062])
DEFAULT_DC_DB = -9.751610265
DEFAULT_AC_DB = 14.248389735
DEFAULT_CLOCK_PS = 4069.0


def read_touchstone(path):
    """Read Touchstone 1.x four-port rows (MA, RI or DB), with comments."""
    options = ['ghz', 's', 'ma', 'r', '50']
    numbers = []
    for raw in Path(path).read_text(encoding='utf-8-sig').splitlines():
        line = raw.split('!', 1)[0].strip()
        if not line:
            continue
        if line.startswith('#'):
            options = line[1:].lower().split()
        elif line.startswith('['):
            raise ValueError('Only four-port Touchstone 1.x is supported')
        else:
            numbers.extend(float(x.replace('D', 'E').replace('d', 'e')) for x in line.split())
    if not numbers or len(numbers) % 33:
        raise ValueError('A four-port record must contain 33 numbers')
    if 's' not in options:
        raise ValueError('Only scattering parameters are supported')
    units = {'hz': 1, 'khz': 1e3, 'mhz': 1e6, 'ghz': 1e9}
    unit = next((x for x in options if x in units), 'ghz')
    fmt = next((x for x in options if x in ('ma', 'ri', 'db')), 'ma')
    z0 = float(options[options.index('r') + 1]) if 'r' in options else 50.0
    a = np.asarray(numbers).reshape(-1, 33)
    f = a[:, 0] * units[unit]
    q = a[:, 1:].reshape(-1, 4, 4, 2)
    if fmt == 'ri':
        s = q[..., 0] + 1j * q[..., 1]
    else:
        mag = 10 ** (q[..., 0] / 20) if fmt == 'db' else q[..., 0]
        s = mag * np.exp(1j * np.deg2rad(q[..., 1]))
    if len(f) < 2 or not np.isfinite(a).all() or f[0] != 0 or not np.all(np.diff(f) > 0):
        raise ValueError('Finite increasing frequency samples including DC are required')
    if not np.allclose(np.diff(f), f[1] - f[0], rtol=1e-9, atol=1e-4):
        raise ValueError('The frequency grid must be uniform')
    return f, s, z0


def mixed_mode(s):
    """Pairs (1,3),(2,4), mode order c1,d1,c2,d2, 100/25 ohm refs."""
    m = np.array([[1, 0, 1, 0], [1, 0, -1, 0],
                  [0, 1, 0, 1], [0, 1, 0, -1]]) / np.sqrt(2)
    return m @ np.asarray(s) @ m.T


def ctle_response(f, dc_db=DEFAULT_DC_DB, ac_db=DEFAULT_AC_DB, poles_ghz=(40., 40., 40.)):
    """One zero / three poles. ac_db is a numerator parameter, not peak gain."""
    poles = 2 * np.pi * np.asarray(poles_ghz) * 1e9
    dc, ac = 10 ** (dc_db / 20), 10 ** (ac_db / 20)
    z = dc * poles[0] / ac
    s = 2j * np.pi * np.asarray(f)
    return ac * poles[1] * poles[2] * (s + z) / np.prod(s[..., None] + poles, axis=-1)


def ctle_peak(dc_db=DEFAULT_DC_DB, ac_db=DEFAULT_AC_DB, poles_ghz=(40., 40., 40.)):
    result = minimize_scalar(lambda x: -abs(ctle_response(10. ** x, dc_db, ac_db, poles_ghz)),
                             bounds=(3, 13), method='bounded', options={'xatol': 1e-12})
    return float(20 * np.log10(-result.fun)), float(10. ** result.x)


def impulse_response(f, transfer, samples_per_ui=50, alpha=ALPHA, duration=TRUNCATION):
    """Discrete convolution kernel; real Hermitian IDFT, then causal truncation."""
    dt = UI / samples_per_ui
    df = f[1] - f[0]
    n = round(1 / (dt * df))
    if not np.isclose(n * dt * df, 1, rtol=1e-11):
        raise ValueError('Time and frequency grids do not form an integral FFT length')
    if f[-1] > .5 / dt or len(f) > n // 2 + 1:
        raise ValueError('Measured band exceeds the time-domain Nyquist frequency')
    win = np.i0(np.pi * alpha * np.sqrt(np.maximum(0, 1 - (f / f[-1]) ** 2))) / np.i0(np.pi * alpha)
    spectrum = np.zeros(n // 2 + 1, complex)
    spectrum[:len(f)] = transfer * win
    spectrum[0] = spectrum[0].real
    full = irfft(spectrum, n=n)
    stop = round(duration / dt)
    if not 0 < stop < n:
        raise ValueError('Impulse duration must be shorter than 1/df')
    return full[:stop], full, win


def trapezoid(samples_per_ui=50):
    """A +/-1 symbol produces this causal +2 V single-bit source pulse.

    Samples at 0 and UI+rise are zero. A one-sample shift gives the requested
    ramp with x(0)=0, x(rise)=2, x(UI)=2 and x(UI+rise)=0.
    """
    nr = round(RISE / UI * samples_per_ui)
    if nr < 1 or not np.isclose(nr * UI / samples_per_ui, RISE):
        raise ValueError('Rise time must be an integer number of samples')
    return np.r_[0., VPP * np.convolve(np.ones(samples_per_ui), np.ones(nr) / nr), 0.]


def apply_ffe(sbr, taps, samples_per_ui):
    taps = np.asarray(taps, float)
    if not 1 <= len(taps) <= 5 or not np.isclose(np.abs(taps).sum(), 1, atol=1e-12):
        raise ValueError('FFE needs 1..5 taps with sum(abs(taps))=1')
    out = np.zeros(len(sbr) + (len(taps) - 1) * samples_per_ui)
    for j, c in enumerate(taps):
        out[j * samples_per_ui:j * samples_per_ui + len(sbr)] += c * sbr
    return out


def symbol_waveform(symbols, sbr, samples_per_ui):
    """Linear convolution using a polyphase FFT, with no circular wraparound."""
    rows = len(symbols) + (len(sbr) - 1) // samples_per_ui
    nfft = next_fast_len(rows)
    spec = rfft(symbols, n=nfft)
    output = np.empty((rows, samples_per_ui))
    for p in range(samples_per_ui):
        sub = sbr[p::samples_per_ui]
        output[:, p] = irfft(spec * rfft(sub, n=nfft), n=nfft)[:rows]
    return output.ravel()[:(len(symbols) - 1) * samples_per_ui + len(sbr)]


def dfe_decisions(samples, taps):
    """Actual recursive +/-1 slicer. Transmitted symbols are never consulted."""
    taps = np.asarray(taps, float)
    if len(taps) > 5:
        raise ValueError('DFE may use at most five taps')
    samples = np.asarray(samples, float)
    decisions = np.empty(len(samples), np.int8)
    feedback = np.zeros(len(samples))
    for n, value in enumerate(samples):
        q = min(n, len(taps))
        if q:
            feedback[n] = np.dot(taps[:q], decisions[n-q:n][::-1])
        decisions[n] = 1 if value - feedback[n] >= 0 else -1
    return decisions, feedback


def _opening_width(phase_ps, one_lower, zero_upper, center):
    """Contiguous zero-threshold opening; interpolate its limiting margins."""
    margin = np.minimum(one_lower, -zero_upper)
    if margin[center] <= 0:
        return 0., float(phase_ps[center]), float(phase_ps[center])
    l = r = center
    while l > 0 and margin[l - 1] > 0:
        l -= 1
    while r < len(margin) - 1 and margin[r + 1] > 0:
        r += 1
    def root(a, b):
        return float(phase_ps[a] - margin[a] * (phase_ps[b] - phase_ps[a]) / (margin[b] - margin[a]))
    left = root(l - 1, l) if l > 0 else float(phase_ps[l])
    right = root(r, r + 1) if r < len(margin) - 1 else float(phase_ps[r])
    return right - left, left, right


def measure_eye(y, symbols, clock_sample, samples_per_ui, start_symbol, count,
                feedback=None, voltage_bins=600, histogram=True):
    """Use every retained symbol and every sample; no voltage percentiles.

    Feedback is held on [clock-UI/2, clock+UI/2) for the current symbol.
    The right endpoint is the next symbol's hold value, so discontinuities
    cannot artificially extend the measured DFE eye.
    """
    p = samples_per_ui
    offsets = np.arange(-p // 2, p // 2 + 1)
    idx = np.arange(start_symbol, start_symbol + count)
    truth = symbols[idx]
    positive = truth == 1
    lower, upper = [], []
    vmin = float(np.min(y)); vmax = float(np.max(y))
    if feedback is not None:
        limit = float(np.max(np.abs(feedback)))
        vmin -= limit; vmax += limit
    padding = max(.01, (vmax-vmin)*.02)
    edges = np.linspace(vmin-padding, vmax+padding, voltage_bins+1)
    hist = np.zeros((voltage_bins, len(offsets)), np.int64) if histogram else None
    centers = clock_sample + idx * p
    error_count = 0
    for j, offset in enumerate(offsets):
        voltage = y[centers + offset].copy()
        if feedback is not None:
            correction_index = idx + (1 if offset == p//2 else 0)
            voltage -= feedback[correction_index]
        lower.append(np.min(voltage[positive])); upper.append(np.max(voltage[~positive]))
        if offset == 0:
            error_count = int(np.count_nonzero(np.where(voltage >= 0, 1, -1) != truth))
        if histogram:
            hist[:, j] = np.histogram(voltage, bins=edges)[0]
    lower, upper = np.array(lower), np.array(upper)
    phase_ps = offsets * UI / p * 1e12
    center = p//2
    width, left, right = _opening_width(phase_ps, lower, upper, center)
    metric = {
        'eye_height_v': float(lower[center] - upper[center]),
        'eye_width_ps': float(width), 'left_crossing_ps':left, 'right_crossing_ps':right,
        'lower_v':float(lower[center]), 'upper_v':float(upper[center]),
        'decision_errors':error_count, 'clock_sample':int(clock_sample),
        'clock_time_ps':float(clock_sample * UI / p * 1e12),
        'clock_phase_ps':float((clock_sample % p) * UI / p * 1e12),
        'measured_symbols':int(count),
    }
    arrays = {'phase_ps':phase_ps, 'one_lower_v':lower, 'zero_upper_v':upper}
    if histogram:
        arrays.update(histogram=hist, voltage_edges_v=edges)
    return metric, arrays


def build_models(channel, samples_per_ui=50, taps=DEFAULT_FFE, dc_db=DEFAULT_DC_DB, ac_db=DEFAULT_AC_DB):
    f, s, z0 = read_touchstone(channel)
    if not np.isclose(z0, 50):
        raise ValueError('This assignment requires 50 ohm single-ended ports')
    sm = mixed_mode(s)
    hch = sm[:,3,1] / 2
    ctle = ctle_response(f, dc_db, ac_db)
    pulse = trapezoid(samples_per_ui)
    h0, full0, window = impulse_response(f, hch, samples_per_ui)
    h1, full1, _ = impulse_response(f, hch * ctle, samples_per_ui)
    sbr0 = fftconvolve(h0, pulse)
    sbr1 = fftconvolve(h1, pulse)
    sbr2 = apply_ffe(sbr1, taps, samples_per_ui)
    return {
        'f':f, 's21':s[:,1,0], 'sdd21':sm[:,3,1], 'ctle':ctle,
        'h_channel':h0,'h_ctle':h1,'pulse':pulse,
        'sbr_channel':sbr0,'sbr_ctle':sbr1,'sbr_ffe':sbr2,
        'window':window, 'full_h_channel':full0,'full_h_ctle':full1,
    }


def run_record(models, samples_per_ui=50, n_bits=2**18, seed=0, clock_offset=0,
               histogram=True, all_stages=True):
    guard = 1024
    total = n_bits + 2*guard
    symbols = (2*np.random.default_rng(seed).integers(0, 2, total)-1).astype(np.int8)
    stages, arrays = {}, {}
    names = ('channel','ctle','ffe') if all_stages else ('ffe',)
    snapshot_start = (guard + 256)*samples_per_ui
    snapshot_count = 40*samples_per_ui
    arrays['waveform_time_ns'] = np.arange(snapshot_start, snapshot_start+snapshot_count)*UI/samples_per_ui*1e9
    for name in names:
        sbr = models['sbr_'+name]
        clock = int(np.argmax(sbr)) + (clock_offset if name == 'ffe' else 0)
        y = symbol_waveform(symbols, sbr, samples_per_ui)
        metric, eye = measure_eye(y, symbols, clock, samples_per_ui, guard, n_bits, histogram=histogram)
        stages[name] = metric
        arrays.update({name+'_'+key:value for key,value in eye.items()})
        arrays['waveform_'+name+'_v'] = y[snapshot_start:snapshot_start+snapshot_count]
        if name == 'ffe':
            taps = np.array([sbr[clock+i*samples_per_ui] for i in range(1,6)])
            decision_samples = y[clock+np.arange(total)*samples_per_ui]
            decisions, feedback = dfe_decisions(decision_samples, taps)
            metric, eye = measure_eye(y, symbols, clock, samples_per_ui, guard, n_bits, feedback, histogram=histogram)
            stages['dfe'] = metric
            arrays.update({'dfe_'+key:value for key,value in eye.items()})
            positions = np.arange(snapshot_start,snapshot_start+snapshot_count)
            current = np.floor_divide(positions-clock+samples_per_ui//2,samples_per_ui)
            arrays['waveform_dfe_v'] = y[snapshot_start:snapshot_start+snapshot_count]-feedback[current]
            arrays['dfe_taps_v'] = taps
            arrays['symbols'] = symbols
            arrays['dfe_decisions'] = decisions
    return stages, arrays


def plot_results(output, models, stages, arrays, samples_per_ui):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,
                         'grid.alpha':.25,'figure.dpi':140,'savefig.dpi':180})
    labels = {'channel':'Channel only','ctle':'Channel + CTLE','ffe':'FFE + channel + CTLE','dfe':'FFE + channel + CTLE + DFE'}
    colors = {'channel':'#4c78a8','ctle':'#f58518','ffe':'#54a24b','dfe':'#b279a2'}
    def save(fig,name):
        fig.tight_layout(); fig.savefig(output/name); plt.close(fig)
    def db(v):return 20*np.log10(np.maximum(abs(v),1e-12))
    fig,ax=plt.subplots(figsize=(7.2,4.3))
    ax.plot(models['f']/1e9,db(models['s21']),label='Single-ended S21')
    ax.plot(models['f']/1e9,db(models['sdd21']),label='Differential Sdd21')
    ax.set(xlabel='Frequency (GHz)',ylabel='Magnitude (dB)',title='Single-ended and differential transmission');ax.legend()
    save(fig,'sparameters.png')
    fig,ax=plt.subplots(figsize=(7.2,4.3))
    ax.plot(models['f']/1e9,db(models['ctle']),label='CTLE')
    ax.plot(models['f']/1e9,db(models['sdd21']/2),label='Matched channel: Sdd21 / 2')
    ax.plot(models['f']/1e9,db(models['ctle']*models['sdd21']/2),label='Matched channel + CTLE')
    ax.axhline(6,ls=':',color='black',label='CTLE gain limit')
    ax.set(xlabel='Frequency (GHz)',ylabel='Magnitude (dB)',title='CTLE and channel frequency responses',ylim=(-70,10));ax.legend()
    save(fig,'frequency_response.png')
    fig,ax=plt.subplots(figsize=(7.2,4.3))
    for name in ('channel','ctle','ffe'):
        sbr=models['sbr_'+name]
        ax.plot(np.arange(len(sbr))*UI/samples_per_ui*1e9,sbr,label=labels[name],color=colors[name])
    clock=stages['ffe']['clock_sample'];tapindices=clock+np.arange(6)*samples_per_ui
    ax.plot(tapindices*UI/samples_per_ui*1e9,models['sbr_ffe'][tapindices],'o',color=colors['dfe'],label='Main / DFE cursors')
    ax.set(xlabel='Time (ns)',ylabel='Differential voltage (V)',title='Consistent finite-rise single-bit responses',xlim=(3.8,4.4));ax.legend()
    save(fig,'sbr.png')
    for name in labels:
        hist=arrays[name+'_histogram'];edges=arrays[name+'_voltage_edges_v'];phase=arrays[name+'_phase_ps']
        probability=hist/hist.sum(axis=0,keepdims=True)
        fig,ax=plt.subplots(figsize=(7.2,4.3))
        im=ax.imshow(probability,origin='lower',extent=[phase[0],phase[-1],edges[0],edges[-1]],aspect='auto',cmap='viridis',norm=LogNorm(vmin=1/stages[name]['measured_symbols'],vmax=max(.1,float(probability.max()))),interpolation='nearest')
        ax.plot(phase,arrays[name+'_one_lower_v'],color='white',lw=.8,ls='--',label='Conditional eye boundaries')
        ax.plot(phase,arrays[name+'_zero_upper_v'],color='white',lw=.8,ls='--')
        ax.axvline(0,color='white',lw=.8,ls=':')
        ax.set(xlabel='Time relative to decision clock (ps)',ylabel='Differential voltage (V)',title=labels[name]+'\nEH %.2f mV | EW %.2f ps'%(1000*stages[name]['eye_height_v'],stages[name]['eye_width_ps']))
        fig.colorbar(im,ax=ax,label='Probability per voltage bin');save(fig,'eye_'+name+'.png')
    fig,ax=plt.subplots(figsize=(9,4.7))
    for name in labels:
        ax.plot(arrays['waveform_time_ns'],arrays['waveform_'+name+'_v'],label=labels[name],color=colors[name],lw=.95)
    ax.set(xlabel='Absolute time (ns)',ylabel='Differential voltage (V)',title='Four output waveforms for the same transmitted bitstream')
    ax.legend(ncol=2,loc='upper center',bbox_to_anchor=(.5,1.0));save(fig,'waveforms.png')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--channel',type=Path,default=Path(__file__).resolve().parents[1]/'instructions'/'HW10_complicated_diff_channel.s4p')
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'results')
    parser.add_argument('--bits',type=int,default=2**18)
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--samples-per-ui',type=int,default=50)
    parser.add_argument('--ffe',type=float,nargs='+',default=DEFAULT_FFE.tolist())
    parser.add_argument('--dc-db',type=float,default=DEFAULT_DC_DB)
    parser.add_argument('--ac-db',type=float,default=DEFAULT_AC_DB)
    parser.add_argument('--clock-time-ps',type=float,default=DEFAULT_CLOCK_PS,help='Absolute decision time for source symbol 0; remains fixed under time-grid refinement')
    parser.add_argument('--clock-offset',type=int,default=None,help='Override clock time with a sample offset from the FFE SBR maximum')
    parser.add_argument('--no-plots',action='store_true')
    parser.add_argument('--validation',action='store_true',help='Run additional seeds and doubled time resolution')
    args=parser.parse_args(argv)
    if args.samples_per_ui%2 or args.bits<1024:parser.error('Use even samples/UI and at least 1024 measured bits')
    args.output.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    models=build_models(args.channel,args.samples_per_ui,np.array(args.ffe),args.dc_db,args.ac_db)
    peak_sample=int(np.argmax(models['sbr_ffe']))
    clock_offset=(round(args.clock_time_ps*1e-12/(UI/args.samples_per_ui))-peak_sample
                  if args.clock_offset is None else args.clock_offset)
    physical_clock=(peak_sample+clock_offset)*UI/args.samples_per_ui
    if peak_sample+clock_offset < args.samples_per_ui//2:
        parser.error('Decision clock must allow a complete half-UI before symbol 0')
    stages,arrays=run_record(models,args.samples_per_ui,args.bits,args.seed,clock_offset)
    peak,peakf=ctle_peak(args.dc_db,args.ac_db)
    dfe=stages['dfe']
    metrics={
        'parameters':{'bit_rate_bps':BIT_RATE,'ui_ps':UI*1e12,'rise_ps':RISE*1e12,'driver_vpp_v':VPP,'source_differential_levels_v':[-VPP,VPP],'single_ended_reference_ohm':50,'differential_reference_ohm':100,'common_mode_reference_ohm':25,'samples_per_ui':args.samples_per_ui,'dt_ps':UI/args.samples_per_ui*1e12,'kaiser_alpha':ALPHA,'impulse_duration_ns':TRUNCATION*1e9,'n_bits':args.bits,'seed':args.seed,'warmup_symbols':1024,'trailing_guard_symbols':1024,'total_generated_symbols':args.bits+2048,'frequency_step_hz':float(models['f'][1]),'frequency_max_hz':float(models['f'][-1]),'source_sha256':hashlib.sha256(args.channel.read_bytes()).hexdigest()},
        'equalizers':{'ctle':{'poles_ghz':[40.,40.,40.],'zeros_ghz':[40*10**((args.dc_db-args.ac_db)/20)],'dc_gain_db':args.dc_db,'numerator_gain_db':args.ac_db,'peak_gain_db':peak,'peak_frequency_ghz':peakf/1e9},'ffe':{'taps':args.ffe,'l1_norm':float(np.abs(args.ffe).sum())},'dfe':{'taps_v':arrays['dfe_taps_v'].tolist(),'feedback_operation':'subtract','threshold_v':0.,'slicer_levels':[-1,1],'decision_phase_sample':int(dfe['clock_sample']%args.samples_per_ui),'clock_sample':dfe['clock_sample'],'feedback_hold':'from half a UI before through half a UI after each decision; right endpoint belongs to the next hold'}},
        'stages':stages,
        'design':{'training_seed':70,'training_bits':16384,'method':'Finite-training linear programming of FFE coefficients and clock search, followed by separate full-record verification','objective':'Maximize the minimum signed central voltage, with signed margins >=4mV at +/-8ps and >=15mV at +/-4ps on the training record','selected_clock_time_ps':physical_clock*1e12,'held_out_seeds':[0,1,2],'ctle_scaling':'Original -10/14dB numerator settings shifted together to a continuous peak of approximately 5.99dB; zero frequency unchanged','ffe':'Two nonzero taps after removing zero leading/trailing taps from the five-tap search'},
        'measurement':{'eye_height':'minimum logic-one voltage minus maximum logic-zero voltage at the selected decision clock; every retained symbol is included, with no percentile clipping','eye_width':'contiguous interval containing the clock for which min(logic-one)>0 and max(logic-zero)<0; boundaries linearly interpolated from the more restrictive conditional margin','eye_histogram':'600 voltage bins for each of samples_per_ui+1 offsets spanning one UI','waveform':'causal linear-ramp samples, real Hermitian channel IDFT with Kaiser1, 20ns causal impulse truncation, polyphase linear convolution','dfe':'recursive decisions from corrected samples only; transmitted bits used only for error and eye measurement','scope':'Finite noiseless random-pattern evidence; no BER extrapolation, jitter model, ADS validation, or all-pattern worst-case guarantee'},
        'validation':{},
        'compliance':{'ctle_poles_pass':True,'ctle_peak_gain_pass':peak<=6+1e-9,'ctle_max_pole_pass':True,'ffe_taps_pass':len(args.ffe)<=5,'ffe_l1_pass':bool(np.isclose(np.abs(args.ffe).sum(),1)),'dfe_taps_pass':len(arrays['dfe_taps_v'])<=5,'eye_height_pass':dfe['eye_height_v']>.120,'eye_width_pass':dfe['eye_width_ps']>15,'decision_error_pass':dfe['decision_errors']==0}
    }
    if args.validation:
        checks=[]
        for seed in (1,2):
            st,_=run_record(models,args.samples_per_ui,args.bits,seed,clock_offset,histogram=False,all_stages=False)
            passed=st['dfe']['eye_height_v']>.120 and st['dfe']['eye_width_ps']>15 and st['dfe']['decision_errors']==0
            checks.append({'seed':seed,'n_bits':args.bits,'pass':passed,**st['dfe']})
        refined=build_models(args.channel,args.samples_per_ui*2,np.array(args.ffe),args.dc_db,args.ac_db)
        refined_offset=round(physical_clock/(UI/(args.samples_per_ui*2)))-int(np.argmax(refined['sbr_ffe']))
        st,_=run_record(refined,args.samples_per_ui*2,args.bits,args.seed,refined_offset,histogram=False,all_stages=False)
        passed=st['dfe']['eye_height_v']>.120 and st['dfe']['eye_width_ps']>15 and st['dfe']['decision_errors']==0
        metrics['validation']={'seed_checks':checks,'refinement':{'samples_per_ui':args.samples_per_ui*2,'dt_ps':UI/(args.samples_per_ui*2)*1e12,'n_bits':args.bits,'pass':passed,**st['dfe']}}
        metrics['compliance']['validation_records_pass']=all(c['pass'] for c in checks) and passed
    compact={key:value for key,value in models.items() if not key.startswith('full_h')}
    compact.update(arrays)
    np.savez_compressed(args.output/'numerical_results.npz',**compact)
    if not args.no_plots:plot_results(args.output,models,stages,arrays,args.samples_per_ui)
    metrics['runtime_seconds']=time.perf_counter()-started
    (args.output/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'stages':stages,'compliance':metrics['compliance'],'runtime_seconds':metrics['runtime_seconds']},indent=2),flush=True)
    return 0 if not args.validation or all(metrics['compliance'].values()) else 2


if __name__=='__main__':
    raise SystemExit(main())
