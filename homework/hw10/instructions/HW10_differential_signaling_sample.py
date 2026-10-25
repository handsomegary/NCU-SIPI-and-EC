#%% HW 10 sample code

#%% define necessary functions here, including
# imp_resp
# get_sbr
# get_cursor
# fixed_rate_unifdt_stateye
# ZFtap
# FFEfreqresp

#%% define function: single-ended to mixed mode

def SE2MM(spSE):
    '''
    Convert 4-port single-ended S parameters to mixed-mode.
    Ports 1&3 and 2&4 are treated as differential pairs.
    
    Parameters
    ----------
    spSE : np array of shape (Nfreq, 16)
        single-ended S parameters

    Returns
    -------
    spMM : np array of shape (Nfreq, 16)
        mixed-mode S parameters
    Each row of spMM is the flattened array of the matrix
    [[Scc11 Scd11 Scc12 Scd12]
     [Sdc11 Sdd11 Sdc12 Sdd12]
     [Scc21 Scd21 Scc22 Scd22]
     [Sdc21 Sdd21 Sdc22 Sdd22]]
    
    '''
    
    import numpy as np
    
    # define conversion matrix (M and M_inverse)
    Minv = np.array(
        
        
    M = 
    
    spMM = np.zeros((spSE.shape),
    
    for k in range(
        spMM[k,:] = (
    
    return spMM




#%% define functions: sectional convolution with Tx FFE, for both SE and DP


def sctn_conv_random_NRZ_TxFFE_SEDP(H, df, T2, Tb, Tr, NpUI, cFFE=[1], Vpp=1,
                            Nbp:int=2**16, sig:str='SE', a=0.5, seed=0):
    '''
    Perform sectional convolution to obtain output response of a given channel
    with transfer function H, assuming the input is purely random NRZ with the
    given parameters. 
    
    The Tx FFE is specified by the variable 'cFFE'. To consider the effects of 
    FFE, there are two ways: (1) construct the corresponding time-domain input
    waveform of FFE and than do sectional convolution. (2) compute the 
    frequency response of the FFE and multiply with the channel's H to obtain 
    the overall transfer function, and than assume the input is ideal NRZ. In
    this homework, we use method (1), because it is more flexible and can be
    extended to include other nonideal effects such as Tx jitter. The tap 
    spacing is fixed to be Tb, which is the most common case.
    
    This function can be used for both single-ended (SE) case and differential
    pair (DP), specified by the input parameter 'sig'.
    
    Parameters
    ----------
    H : Overall channel transfer function. Must have uniform df and include DC
    df: frequency step of H in Hz
    cFFE: FFE taps, list or numpy 1-dim array,  [c0, c1, c2, ...]
       Note: c0 is not necessarily the main cursor tap. Instead,
       c0 is the "left-most" cursor tap. It may be the main cursor
       or it may be a precursor tap. If it is the 1st precursor, 
       than c1 is the main cursor, c2 is the 1st post cursor, and
       so on. If cFFE=[1], then there is no FFE. The absolute values of
       the elements in cFFE are not important; they will be normalized
       to meet the Vpp constraint. Only the relative ratios in cFFE are
       relevant.
    T2: truncation time of channel impulse response in seconds
    Tb: unit interval in seconds
    Tr: risefall time in seconds
    NpUI: number of points per unit interval
    Vpp: peak-to-peak voltage of NRZ signal. The post-FFE waveform
         at Tx still has peak-to-peak voltage of Vpp. It may not be
         in the range [0,Vpp]. Instead, it may be in the range of
         [-q,Vpp-q] for some positive number q, depending on cFFE.
         
         Important: for differential signaling (sig='DP'), this Vpp
         is the peak-to-peak voltage of each single-ended driver. 
         Therefore, the actual p2p voltage of the differential signal
         at Tx will be 2*Vpp.
    
    Nbp: minimum number of random bits to simulate
    sig: signaling scheme, a string, must be either 'SE' or 'DP'.
    a: Kaiser window parameter
    seed: seed value for random number generator

    Returns
    -------
    ty: time vector of y
    y: total output waveform
    txk: time vector of xk
    xk: the input waveform of the kth section (for checking)
    bp: the overall input bit pattern for validation if needed
        It is always in {0,1}, regardless of 'sig' (not {1,-1}).

    '''
    import numpy as np
    from numpy import pi, sqrt
    from scipy.fft import fft, ifft
    from scipy.special import i0
    from time import time
    
    rng = np.random.default_rng(seed=seed)
    
    # There is no hint below. Try your best to figure it out.
    # In fact, this function is pretty much the same as the single-ended 
    # version which we have written in previous homework. You just have 
    # to do a little modification so that the case of differential signaling 
    # can also be simulated. When the input argument sig == 'SE', this 
    # function should be completely the same as before, in which the input
    # NRZ signal swings between 0 and Vpp. When sig == 'DP', the differential
    # input is the difference between two single-ended sources, which thus
    # swings between -Vpp and +Vpp. That's the only difference that you
    # need to implement.
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    y = y[Nbps*NpUI:Nbps*Ns*NpUI] # truncate y (掐頭去尾)
    ty = ty[Nbps*NpUI:Nbps*Ns*NpUI] # truncate ty
    # This is why the 'Ns' parameter was set to Nbp//Nbps+2. 
    # We add two more sections, then truncated.
    
    txk = np.arange(xk.size)*dt
    print('\nTotal convolution time: %d sec'%(time()-tic))

    return ty, y, txk, xk, bp
    # return xk and bp for checking


#%% define function: simulate fixed-rate DFE output for both SE and DP

def fixed_rate_DFE_SEDP(y,NpUI,c,Vt,sig:str='SE',sample_point=0):
    '''
    Simulate a decision feedback equalizer assuming input has fixed data rate.
    
    The slicer output of the DFE is [0,1] if sig=='SE' (single-ended), or
    [-1,1] if sig=='DP'. 
    
    For differential signaling, Vt is usually 0.
    
    The user has to supply a sample_point for this function, which is the point
    where the DFE makes decision (logic 0 or 1). The DFE will make decisions at
    points of the form (sample_point + i*NpUI) for all possible i (>0 and <0).
    
    This code assumes the sample_point is located at the center of the eye. In
    other words, the clock trigger points are assumed to be at (sample_point 
    + NpUI//2 + i*NpUI) for all possible i (>0 and <0).
    
    Parameters
    ----------
    y : 1-dim np array
        The input waveform at the front of DFE. Should be long enough.
        Should have uniform time step.
    NpUI : int
        Number of points per UI
    c : DFE tap values
        Python list or 1-dim np array
    Vt : float
        threshold voltage for decision making
    sig: signaling scheme, a string, must be either 'SE' or 'DP'.
    sample_point : int
        The sampling point at which the DFE makes decisions.
        The DFE will make decisions at (sample_point + i*NpUI).

    Returns
    -------
    z : 1-dim np array, same length as y
        Post-DFE waveform, i.e., waveform at the summer output of the DFE.

    '''
    import numpy as np
    
    # No hint here. Again, this function is pretty much the same as the
    # single-ended version. There is one important thing that you should be
    # careful: when the decision is False, for single-ended case you don't
    # need to do anything, but for differential-pair, the slicer output is
    # -1 and so you still need to do some addition to the waveform.
    
    return z




#%%
import 
import 
from 
def db(

f, spSE, _ = snpreader('HW10

# convert single-ended to mixed mode
spMM = SE2MM(

#%% compare single-ended S21 and Sdd21
plt.figure(
plt.plot(
plt.plot(
plt.xlabel(
plt.ylabel(
plt.legend(
plt.grid(

#%% get transfer function and evaluate impulse response
H_ch = 
# transfer function of channel only, including the effect of impedance division


# check impulse response
t, h = imp_resp(f,H_ch)


#%% set parameters
df = 
T2 =  # truncation length (in seconds) of impulse response
Tb = 
Tr = 
NpUI = 
Vpp =  # p2p voltage of each single-ended driver, before impedance division
Nbp = 
sig = 
aK =  # Kaiser parameter

#%% no any equalizer

# sectional convolution and build stateye
ty, y1 = sctn_conv_random_NRZ_TxFFE_SEDP(

# build statistical eye
stateye1, extent1 = 

fig, axs = plt.subplots(
im = axs.imshow(
fig.colorbar(
axs.set_xlabel(
axs.set_ylabel(
axs.set_title(
axs.grid(
axs.set_ylim(


#%% Rx CTLE design

from numpy import pi

w = 2 * pi * f

# === CTLE === #

# You can use any kind of CTLE function.
# The design constraints are detailed in the homework description.



S21CTLE = 

H_ch_CTLE =  # transfer function of channel+CTLE


# sectional convolution and build stateye
y2 = sctn_conv_random_NRZ_TxFFE_SEDP(


# build statistical eye
stateye2, extent2 = 

#%% plot stateye
%matplotlib inline
fig, ax = plt.subplots(
ax.plot(f*1e-9,db(S21CTLE),label=r'$S_\mathrm{dd21}^\mathrm{CTLE}$')
ax.plot(f*1e-9,db(Vpp*H_ch),label=r'$S_\mathrm{dd21}^\mathrm{channel}$')
ax.plot(f*1e-9,db(Vpp*H_ch_CTLE),label=r'$S_\mathrm{dd21}^\mathrm{ch+CTLE}$')
ax.grid(
ax.set_xlabel(
ax.set_ylabel(
ax.legend(
ax.set_ylim(

fig, ax = plt.subplots(
im = ax.imshow(
fig.colorbar(
ax.set_xlabel(
ax.set_ylabel(
ax.set_title(
ax.grid(
ax.set_ylim(
ax.set_title(




#%% get SBR for FFE design
# SBR of channel-only (this is just for comparison, not for FFE design)
tsbr, sbr_ch = get_sbr(f,H_ch,

# SBR of channel+CTLE (this is for FFE design)
tsbr, sbr_ch_CTLE = get_sbr(f,


# get cursor values, design FFE taps by zero-forcing, 
# and compute equivalent FFE frequency response

Npre = 
Npost = 
a = get_cursor(
cFFE = ZFtap(

# Instead, you may consider manually set FFE taps yourself
cFFE =

H_FFE = FFEfreqresp(
H_FFE_ch_CTLE =  # transfer function of FFE+Channel+CTLE

#%% recheck SBR (this SBR will also be used for DFE design)

# SBR of FFE+channel+CTLE
tsbr, sbr_FFE_ch_CTLE = get_sbr(

%matplotlib qt
plt.plot(tsbr*1e12,sbr_ch)
plt.plot(tsbr*1e12,sbr_ch_CTLE)
plt.plot(tsbr*1e12,sbr_FFE_ch_CTLE)
plt.xlabel('Time (ps)')
plt.ylabel('SBR (V)')
plt.grid(color='lightgray')
plt.xlim(0,10000)


#%% FFE+channel+CTLE simulation

# sectional convolution and build stateye
y3, txk, xk = sctn_conv_random_NRZ_TxFFE_SEDP(

# build statistical eye
stateye3, extent3 = 



#%% plot stateye

%matplotlib inline
fig, ax = plt.subplots(figsize=(6,4),dpi=300)
# fig = plt.figure(dpi=300,figsize=(12,4))
ax.plot(f*1e-9,db(H_FFE),label='FFE')
ax.plot(f*1e-9,db(Vpp*H_ch),label='channel')
ax.plot(f*1e-9,db(Vpp*H_ch_CTLE),label='channel+CTLE')
ax.plot(f*1e-9,db(Vpp*H_FFE_ch_CTLE),label='FFE+channel+CTLE')
ax.grid(color='lightgray')
ax.set_xlabel(
ax.set_ylabel(
ax.legend(
ax.set_ylim(

fig, ax = plt.subplots(
im = ax.imshow(
fig.colorbar(
ax.set_xlabel(
ax.set_ylabel(
ax.set_title(
ax.grid(
ax.set_ylim(


#%% check input waveform
%matplotlib qt
plt.figure(dpi=300,figsize=(7,2))
plt.plot(txk*1e9,xk)
plt.xlim(2,2.2)
plt.xlabel('Time (ns)')
plt.ylabel('$v_\mathrm{in}$ with FFE (V)')
plt.grid()

#%% plot time-domain waveform of FFE+channel+CTLE
tstart = 
tend = 
tidx = (ty>=tstart) & (ty<=tend)
plt.figure(dpi=300,figsize=(9,3))
plt.plot(ty[tidx], y1[tidx],color=   ,linewidth=2,label='channel')
plt.plot(ty[tidx], y2[
plt.plot(ty[tidx], y3[
plt.xlabel(
plt.ylabel(
plt.legend(
plt.grid(




#%% get cursor values of sbr_FFE_ch_CTLE and design DFE

Npre = 
Npost = 
a = get_cursor(
print(a)
cDFE = 
print(cDFE)

Vt =   # threshold voltage
sample_pt = # sampling point



#%% DFE simulation
y4 = fixed_rate_DFE_SEDP(

# build statistical eye
stateye4, extent4 = 


#%% plot stateye

%matplotlib inline
fig, axs = plt.subplots(
im = axs.imshow(
    
    
    
    
    
    
    

#%% plot time-domain waveform
tstart = 
tend = 
tidx = 
plt.
plt.plot(ty[tidx],
plt.plot(
plt.xlabel(
plt.ylabel(
plt.legend(
plt.grid(

