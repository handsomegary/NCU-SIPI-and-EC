

#%% HW 7: CTLE design and simulation


#%% define function: plot impulse response

def imp_resp(f,H,dt=1e-12,a=1):
    '''
    Compute the impulse response of a given channel response

    Parameters
    ----------
    f : frequency vector of H. Must be uniform step and include DC.
    H : channel transfer function.
    dt: output time step
    a: Kaiser window parameter
    
    Returns
    -------
    t : time vector of h
    h : impulse response

    '''
    import 
    from 
    from 
    from 
    import 
    
    df =   # frequency step
    M =    # total number of points for IDFT
    
    W =    # Kaiser window
    W[0] 
    
    h =    # impulse response of H
    t =    # time vector of h
    
    # %matplotlib qt
    
    plt.figure(
    plt.plot(
    plt.xlabel(
    plt.ylabel(
    plt.grid(
    plt.show()
    
    return t, h


#%% define functions: sectional convolution

def sctn_conv_random_NRZ(H,df,T2,Tb,Tr,NpUI,Vpp=1,Nbp:int=2**16,a=0.5,seed=0):
    '''
    Perform sectional convolution to obtain output response of a given channel
    with transfer function H, assuming the input is purely random NRZ with the
    given parameters.
    
    Parameters
    ----------
    H : Overall channel transfer function. Must have uniform df and include DC
    df: frequency step of H in Hz
    T2: truncation time of channel impulse response in seconds
    Tb: unit interval in seconds
    Tr: risefall time in seconds
    NpUI: number of points per unit interval
    Vpp: peak-to-peak voltage of NRZ signal
    Nbp: number of random bits to simulate
    a: Kaiser window parameter
    seed: seed valud for random number generator

    Returns
    -------
    y: total output waveform
    ty: time vector of y

    '''
    import 
    from 
    from 
    from 
    from 
    
    rng = np.random.default_rng(seed=seed)
    
    f =        # frequency vector of H
    dt =       # output time step
    M =        # total number of points for DFT
    Nbps =     # number of bits per section
    Ns =       # number of sections
    NpTr =     # number of points per risetime
    NpT2 =     # number of points for duration T2
    Npos =     # number of points per output section
    rising_ramp =   # use double precision (dtype=float), not float32
    falling_ramp = 
    
    y =        # total output waveform. Use dtype=float, not float32
    # note: after sectional convolution, y will be truncated to length Nbps*Ns*NpUI
    
    W =        # Kaiser window
    W[0] 
    
    tic = time()
    
    print('Sectional convolution progress:    ',end='')
    
    for k in range(Ns): # for each section
        if (k%(Ns//10)) == 0:
            print('\b\b\b%2d%%'%(10*k/(Ns//10)),end='')
        
        # construct input waveform of one section
    
        if k == 0:
            previous_bit 
        else:
        
        # set purely random bpk
        bpk = rng.integers(  
        
        
        xk = 
        
        # no hint here. Same as previous homework.
        
        
        
        
        
        # convolution of one section
        
        y[
    
    y = y[:Nbps*Ns*NpUI] # truncate y
    ty =       # time vector for output y
    
    print('\nTotal convolution time: %d sec'%(time()-tic))

    return ty, y

#%% define function: build fixed-rate stateye

def fixed_rate_unifdt_stateye(y,NpUI:int):
    '''
    Construct the fixed-rate statistical eye assuming y has uniform dt.

    Parameters
    ----------
    y : time domain waveform. ndarray with shape (y.size,)
    NpUI : int. Number of points per unit interval

    Returns
    -------
    stateye
    extent

    '''
    import numpy as np
    
    vmin =     # min voltage of stat eye (must be smaller than np.min(y))
    vmax =     # max voltage of stat eye (must be greater than np.max(y))
    dv = 1e-3  # size of each voltage bin
    Nvolt =    # number of Vbins
    
    # build 2D histogram
    eyehist = np.zeros((Nvolt,NpUI),dtype=
    
    # integer-array indexing ==> No for-loop needed ==> much faster!!
    idxt = np.arange(y.size) % NpUI # array of time-axis indices
    idxv =                          # array of volt-axis indices
    
    np.add.at(eyehist,              # unbuffered in place addition 
    # The above line is the key difference in this code. Previously, we
    # use a for-loop to build the histogram, which is slow. In this code,
    # the for-loop is removed, resulting in much faster construction.
              
    
    # vertical normalization
    total_count = 
    stateye =        # array broadcasting
    
    # construct 2-UI eye
    stateye = 1e-15 + np.hstack
    extent=[0,2,vmin,vmax]
    
    del eyehist, total_count
    
    return stateye, extent



#%% read S2P and plot impulse response

# import necessary modules you need
import ...


def db(x,power: bool=False): return 20/(int(power)+1)*np.log10(abs(x)+1e-15)

f, sp, _ = snpreader(
# remember to set the maxsamplenum parameters

S21 = sp[

plt.plot(
plt.xlabel(
plt.ylabel(


df = 

# check impulse response
t, h = imp_resp(

#%% set parameters
T2 =     # truncation length (in seconds) of impulse response
Tb =  
Tr = 
NpUI = 
Vpp = 

# Note: theoretically, we should set the truncation time T2 after
# CTLE is multiplied to S21, i.e., T2 should be determined based
# on the impulse response of S21*HCTLE, not just S21 alone.
# Nevertheless, in this problem because we already know the CTLE
# does not have time delay, so the T2 based on S21 is valid as 
# well for S21*HCTLE. If the equalizer has a significant time
# delay, then we need to determine T2 based on the impulse response
# of the overall channel, i.e., S21*Hequalizer.

#%% post-channel response and eye

# sectional convolution and build stateye
ty, y = sctn_conv_random_NRZ(

# build statistical eye
stateye, extent = fixed_rate_unifdt_stateye(

plt.figure(
plt.imshow(np.log10(
plt.colorbar(
plt.xlabel(
plt.ylabel(
plt.title(


#%% CTLE design

w = 2 * pi * f

# === CTLE === #
f1G =     # in GHz
f2G =     # in GHz
Adc_dB = 
Aac_dB = 
w1 = 2*pi*f1G*1e9
w2 = 
Adc = 10**(Adc_dB/20)
Aac = 
S21CTLE =   # CTLE's S21

H =      # overall transfer function


# sectional convolution and build stateye
tz, z = 

# build statistical eye
stateye, extent = 

# plot stateye
fig, axs = plt.subplots(1,2,figsize=(12,4),dpi=300,gridspec_kw={'width_ratios': [1, 2]})
# fig = plt.figure(dpi=300,figsize=(12,4))
axs[0].plot(f*1e-9,db(S21CTLE),label='CTLE')
axs[0].plot(f*1e-9,db(S21),label='post-channel')
axs[0].plot(f*1e-9,db(H),label='post-CTLE')
axs[0].grid(
axs[0].set_xlabel(
axs[0].set_ylabel(
axs[0].legend(loc=
axs[0].set_ylim(-30,10)

im = axs[1].imshow(
fig.colorbar(im,
axs[1].set_xlabel(
axs[1].set_ylabel(
axs[1].set_title(
axs[1].grid(

#%% plot time-domain waveform
tstart = 
tend = tstart + 
plt.figure(dpi=300,figsize=(9,2))
plt.plot(ty
plt.plot(tz
plt.xlabel(
plt.ylabel(
plt.legend(

