

#%% HW8: FFE design by zero forcing

#%% define function: get single-bit response

def get_sbr(f,H,Tb,Tr,dt,Vpp=1,a=0.5):
    '''
    Compute the single-bit response of a given channel response

    Parameters
    ----------
    f : frequency vector of H. Must be uniform step and include DC.
    H : channel transfer function.
    Tb: unit interval of the single-bit
    Tr: 0~100% risetime of the input NRZ pulse
    dt: output time step
    a: Kaiser window parameter
    
    Returns
    -------
    tsbr : time vector of sbr
    sbr : single-bit response

    '''
    import numpy as np
    from numpy import pi, sqrt, sinc, exp
    from 
    from 
    
    df =  # freq step
    M =  # total number of points for IDFT
    
    # Kaiiser window
    W = 
    W[0] 
    
    # spectrum of input trpezoidal pulse
    Vin = Vpp*
    
    sbr = 
    tsbr = np.arange(
    
    return tsbr, sbr

#%% define functions: get cursor values

def get_cursor(sbr,NpUI,Npre,Npost):
    '''
    Get the cursor values of the sbr, with Npre pre-cursors,
    Npost post-cursors

    Parameters
    ----------
    sbr : 1-dim np array
    NpUI : number of points in sbr per UI
    Npre : number of precursor taps
    Npost : number of postcursor taps (not counting main cursor)

    The main cursor is placed at the position of max(sbr).
    
    If Npre is larger than the available precursors, then
    zeros will be returned.


    Returns
    -------
    a: 1-dim np array of cursor values, with size Npre+Npost+1

    '''
    import numpy as np
    
    kmain = np.argmax(
    a = np.zeros(
    a[Npre:] = sbr[
    for i in range(1,Npre+1):
        if kmain-i*NpUI >= 0:
            a[Npre-i] = sbr[
        else: break
    
    return a
    
#%% define function: zero forcing solution for FFE taps

def ZFtap(a,Npre):
    '''
    Compute the zero forcing solution for FFE taps given
    the cursor values stored in the array 'a', in which
    Npre cursors are precursors.
    
    Note: the value of Npost is not needed, as all the 
    last (a.size-Npre-1) elements are postcursors.
    
    Parameters
    ----------
    a : 1-dim nparray of cursor values
    Npre: number of precursor taps
    
    Returns
    -------
    c: FFE taps

    '''
    import numpy as np
    N = a.size # total number of FFE taps
    if N == 1:
        c =    # If only one tap, meaning no FFE, what is c?
        print('c =',c)
        return c
    
    A = np.zeros( # construct the A matrix
    for i in range(Npre+1):
        A[:N-Npre+i,i] = a[
    for i in range(Npre+1,N):
        A[
    
    zzf = np.zeros(
    zzf[ ? ] = 1
    c = np.linalg.solve(
    c =    # normalize, such that sum(abs(c)) == 1
    
    
    print('a =',a)
    print('A =')
    print(A)
    print('c =',c)
    
    return c

#%% define function: FFE frequency response

def FFEfreqresp(c,f,Tb):
    '''
    Compute the equivalent frequency response of a given
    FFE taps given by variable 'c' at frequencies 'f'.

    Parameters
    ----------
    c : list of 1-dim np array of FFE taps
    f : 1-dim np array of frequency in Hertz
    Tb: unit interval in seconds
    
    Returns
    -------
    HFFE: equivalent frequency response of the FFE
    '''
    
    import numpy as np
    HFFE = np.ones(f.size,dtype=complex) *
    for i in range(1,len(c)):
        HFFE += 
        
    return HFFE
    


#%% define functions: sectional convolution with Tx FFE

def sctn_conv_random_NRZ_TxFFE(H,df,T2,Tb,Tr,NpUI,c=[1],Vpp=1,Nbp:int=2**16,a=0.5,seed=0):
    '''
    Perform sectional convolution to obtain output response of a given channel
    with transfer function H, assuming the input is purely random NRZ with the
    given parameters.
    
    The Tx FFE is specified by the variable 'c'. To consider the effects of 
    FFE, there are two ways: (1) construct the corresponding time-domain input
    waveform of FFE and than do sectional convolution. (2) compute the 
    frequency response of the FFE and multiply with the channel's H to obtain 
    the overall transfer function, and than assume the input is ideal NRZ.
    
    In this homework, we use method (1), because it is more flexible and can be
    extended to include other nonideal effects such as Tx jitter. Therefore,
    you need to think how to construct the post-FFE waveform in Python in an
    efficient way.
    
    The tap spacing is fixed to be Tb, which is the most common case.
    
    Parameters
    ----------
    H : Overall channel transfer function. Must have uniform df and include DC
    df: frequency step of H in Hz
    c: FFE taps, list or numpy 1-dim array,  [c0, c1, c2, ...]
       Note: c0 is not necessarily the main cursor tap. Instead,
       c0 is the "left-most" cursor tap. It may be the main cursor
       or it may be a precursor tap. If it is the 1st precursor, 
       than c1 is the main cursor, c2 is the 1st post cursor, and
       so on. If c=[1], then there is no FFE. The absolute values of
       the elements in c are not important; they will be normalized
       to meet the Vpp constraint. Only the relative ratios in c are
       relevant.
    T2: truncation time of channel impulse response in seconds
    Tb: unit interval in seconds
    Tr: risefall time in seconds
    NpUI: number of points per unit interval
    Vpp: peak-to-peak voltage of NRZ signal. The post-FFE waveform
         at Tx still has peak-to-peak voltage of Vpp. It may not be
         in the range [0,Vpp]. Instead, it may be in the range of
         [-q,Vpp-q] for some positive number q, depending on c.
    Nbp: minimum number of random bits to simulate
    a: Kaiser window parameter
    seed: seed value for random number generator

    Returns
    -------
    y: total output waveform
    ty: time vector of y

    '''
    import numpy as np
    from numpy import pi, sqrt
    from scipy.fft import fft, ifft
    from scipy.special import i0
    from time import time
    
    rng = np.random.default_rng(seed=seed)
    
    f =  # frequency vector of H
    dt = # output time step
    M =  # total number of points for DFT
    Nbps = # number of bits per section
    Ns = Nbp//Nbps+2 # number of sections
    NpTr =  # number of points per risetime
    NpT2 =  # number of points for duration T2
    Npos = Nbps*NpUI+NpT2 # number of points per output section
    unit_ramp = 1/NpTr * np.arange(NpTr,dtype=float)
    # Note: Vpp is not used in the above definition of ramp, because
    # with FFE, there are many different transitions (jumps). So we just
    # define the unit ramp array, and scale the magnitude on the fly.
    
    y = np.zeros(   # total output waveform
    # note: after convolution, y will be truncated to length (Ns-1)*Nbps*NpUI
    ty = dt*np.arange(y.size,dtype=float) # time vector for this y
    
    W =    # Kaiser window
    W[0] 
    
    bp = Vpp*rng.integers(0,2,Ns*Nbps).astype(float) # complete input bit pattern
    
    # FFE processing
    c = np.array(c) # in case c is a Python list, not ndarray
    NFFE = len(c) # number of FFE taps
    c = c/np.sum(abs(c)) # normalization
    
    bpFFE = np.convolve(?,?,'full')[:Nbps*Ns] # post-FFE input levels
    # the effect of risetime will be added to x in the for loop
    
    tic = time()
    
    print('Sectional convolution progress:     ',end='')
    
    for k in range(Ns): # for each section
        if abs((k%((Ns-1)/10))) < 1:
            # print(k,Ns,Ns/10)
            print('\b\b\b\b%3d%%'%(round(10*k/((Ns-1)/10))),end='')
        
        # construct input waveform of one section
    
        if k == 0:
            previous_bit = 0
        else:
            previous_bit = 
        
        bpk = bpFFE[
        
        xk = 
        
        if previous_bit == bpk[0]: # no transition
            xk[:NpUI] = 
        else: # rising or falling ramp
            xk[:NpTr] = 
            # hint: use previous_bit, bpk[0], and the unit_ramp
            xk[NpTr:NpUI] = 
        
        for j in range(1,Nbps):
            if bpk[?] == bpk[?]:
                xk[?:?] = 
            else:
                xk[?:?] = 
                xk[?:?] = 
        
        # convolution of one section. Same as before
        
        y[
    
    y = y[Nbps*NpUI:Nbps*Ns*NpUI] # truncate y
    ty = ty[Nbps*NpUI:Nbps*Ns*NpUI] # truncate ty
    
    txk = np.arange(xk.size)*dt
    print('\nTotal convolution time: %d sec'%(time()-tic))

    return ty, y, txk, xk
    # return xk for checking



#%% read S2P and plot impulse response
# use the s2p of previous homework

import 
import 
from 
def db()

f, sp, _ = snpreader('HW7_stripline_50cm.s2p', reciproc=True,maxsamplenum=10000)

H =  # channel transfer function

plt.plot(
plt.xlabel(
plt.ylabel(


df = 

# check impulse response
t, h = imp_resp(

#%% set parameters
T2 =  # truncation length (in seconds) of impulse response
Tb = 
Tr = 
NpUI = 
Vpp = 


#%% Case 1: no FFE
c = [1]
# sectional convolution and build stateye
ty, y, txk, xk = sctn_conv_random_NRZ_TxFFE(

# build statistical eye
stateye,  # same as previous homework

plt.figure(
plt.imshow(
plt.colorbar(
plt.xlabel(
plt.ylabel(
plt.title(

#%% check input waveform
# %matplotlib qt
plt.plot(
plt.xlim(
plt.grid(

#%% get SBR for FFE design

tsbr, sbr = get_sbr(
plt.plot(tsbr
plt.xlabel(
plt.ylabel(
plt.grid(

#%% get cursor values, design FFE taps by zero-forcing, 
# and compute equivalent FFE frequency response

Npre = 
Npost = 
a = get_cursor(
c = ZFtap(
# c = [0.65,-0.35] # You may manually set FFE taps if you want
HFFE = FFEfreqresp(


#%% FFE design by zero-forcing

# sectional convolution and build stateye
tz, z, txk, xk = sctn_conv_random_NRZ_TxFFE(

# build statistical eye
stateye,




# plot stateye
# Same as previous homework. Remember to change the trace label.

#%% check input waveform
plt.figure(
plt.plot(
plt.xlim(
plt.xlabel(
plt.ylabel(
plt.grid(

#%% plot time-domain waveform
# same as previous homework




