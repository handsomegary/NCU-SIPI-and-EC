
#%% HW 11: crosstalk simulation of eight coupled lines

#%% define functions

# copy necessary functions from previous homework to here


#%% define functions: sectional convolution with Tx FFE, for both SE and DP,

# SINGLE INPUT MULTIPLE OUTPUTS version, with specified initial time delay

def sctn_conv_random_NRZ_TxFFE_SEDP_SIMO(H, df, T2, Tb, Tr, Td, NpUI, 
            cFFE=[1], Vpp=1, Nbp:int=2**16, sig:str='SE', a=0.5, seed=0):
    '''
    Perform sectional convolution to obtain output responses of a SIMO channel
    with transfer functions given in H, assuming the input is purely random 
    NRZ with the given parameters. 
    
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
    H : channel transfer functions. Must have uniform df and include DC
        Allowable shapes of H:
            (Nfreq,) ==> single-input single-output
            (Nch, Nfreq) ==> single-input multiple-output (Nch outputs)
        where Nfreq = number of frequency points given in H
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
    
    Td: time delay of the driver (scalar)  <== NEW IN THIS VERSION
    
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
    ty: time vector of Y
    Y: total output waveform.
        Shape = (ty.size,) if H.ndim == 1
        Shape = (Nch, ty.size) if H.ndim == 2
    tx: time vector of x
    x: total input waveform
    bp: the overall input bit pattern for validation if needed
        It is always in {0,1}, regardless of 'sig' (not {1,-1}).

    '''
    import numpy as np
    from numpy import pi, sqrt
    from scipy.fft import fft, ifft
    from scipy.special import i0
    from time import time
    
    rng = 
    
    Nfreq = H.shape[ # number of frequency points in H
    if H.ndim == 1:
        Nch =          # number of output channels
    elif H.ndim == 2:
        Nch = 
    
    f =  # frequency vector of H
    dt =  # output time step
    M =  # total number of points for DFT
    Nbps =  # number of bits per section
    Ns =  # number of sections
    NpTr =  # number of points per risetime
    NpTd =  # number of points for Td
    NpT2 =  # number of points for duration T2
    Npos =  # number of points per output section
    unit_ramp = 
    # Note: Vpp is not used in the above definition of ramp, because
    # with FFE, there are many different transitions (jumps). So we just
    # define the unit ramp array, and scale the magnitude on the fly.
    
    
    Y = np.zeros((?,?),dtype=float) # total output waveform
    
    # Y will be truncated to length (Ns-1)*Nbps*NpUI after convolution
    ty =  # time vector for this y
    
    W =  # Kaiser window
    W[0] /=
    
    # FFE processing
    cFFE =  # in case cFFE is a Python list, not ndarray
    NFFE =  # number of FFE taps
    cFFE =  # normalization
    
    # Construct the complete input bit pattern in {0,1}.
    # It is the one that will be returned.
    bp = rng.
    
    # Construct the post-FFE input waveform bpFFE
    # The effect of risetime will be added to x later
    # The effect of Vpp is added to bpFFE now
    if sig=='SE':
        # use Vpp*bp for convolution with cFFE
        bpFFE = np.convolve(
    elif sig=='DP':
        # First convert the {0,1} NRZ signal to {-1,1} form then convolution
        # The overall p2p voltage of the resulting diff signal becomes 2Vpp.
        bpFFE = 
    else:
        raise ValueError(
    
    tic = time()
    
    # build x, the total input waveform, not just one section
    x = np.zeros(
    if bpFFE[0] != 0:
        x[NpTd:NpTd+NpTr] = 
        x[NpTd+NpTr:  ??] = 
    now = NpTd + NpUI
    for j in range(
        if ?? == ??
            x[now:  ?? ] = 
        else:
            x[      = 
            x[      = 
        now += 
    tx =   # time vector for x
    
    print('Sectional convolution progress:     ',end='')
    
    for k in range(Ns): # for each section
        if abs((k%((Ns-1)/10))) < 1:
            # print(k,Ns,Ns/10)
            print('\b\b\b\b%3d%%'%(round(10*k/((Ns-1)/10))),end='')
        
        # get input waveform of one section
        xk = x[k*
        
        # convolution of one section
        Y[:,
    
    Y =   # truncate y (掐頭去尾)
    ty =  # truncate ty
    # This is why the 'Ns' parameter was set to Nbp//Nbps+2. 
    # We add two more sections, then truncated.
    
    txk = np.arange(xk.size)*dt
    print('\nTotal convolution time: %d sec'%(time()-tic))
    
    if H.ndim == 1:
        Y.resize(   # make Y the same ndim as H
        
    return ty, Y, tx, x, bp
    # return x and bp for checking





#%% import modules and define functions
import



#%% read snp
fname = 'HW11_
freq, sp, Z0 = 
df = 

#%% check S parameters
plt.plot(


#%% define parameters
TxPorts = np.array([1,3,5,7,9,11,13,15])-1  # -1 because Python index starts from 0
RxPorts = 

RS = 
CS = 
RL = 
CL = 
Tb = 
Tr = 
Td = np.array([0,  0.5,  ...])*Tb # shape=(8,) array
Vpp = 
NpUI = 
cFFE = 
Nbp = 
sig = 

Nport = # determine number of ports based on sp.shape
Nch =   # number of channels, determined based on Nport
Id =    # identity matrix of size Nport x Nport



#%% Step 1: solve transfer function by MNA

# define the frequency-independent part of the A matrix
A0 = np.zeros((??,??),dtype=??)

# fill in G matrix
A0[TxPorts,TxPorts] = 
A0[??  ,  ??] = 

# fill in the Id in the upper-right block
A0[

# define right-hand side vectors (matrix)
B = np.zeros(
B[   = 1

# allocate transfer function array
H = np.zeros((?,?,?),dtype=?)

for idxf in range(  # frequency sweep
    A = A0.copy()
    
    # fill in C matrix
    A[   += 
    A[   += 
    
    # fill in S matrix
    A[    = Id - 
    A[    = -Z0*
    
    # solve voltage and currents
    X = 
    
    # store the required transfer functions
    H[?,?,?] = 
    
#%% compute impulse response

t, h = ss.imp_resp(freq,H)

#%% check impulse response to determine truncation time (T2)
# %matplotlib qt
plt.plot(t, (h.transpose((2,0,1)).reshape(t.size,Nch*Nch)))
plt.xlim(
    
#%% set truncation time
T2 = 


#%% sectional convolution to obtain output responses

ty, Y, tx, x = sctn_conv_random_NRZ_TxFFE_SEDP_SIMO(H[:,0,:],

#%
for j in range(
    Y += 

#%% build statistical eye

eyes = []
for i in range(Nch):
    stateye, extent = fixed_rate_unifdt_stateye(Y[
    eyes.append([stateye, extent])


#%% plot stateye

ch = 4 # 0 ~ Nch-1

%matplotlib inline
fig, axs = plt.subplots(1,1,figsize=(6,4),dpi=300)
im = axs.imshow(np.log10(eyes[ch][0])




#%% plot time-domain waveform
ch = 4

tstart = 
tend = tstart + 4e-9
tidx = (ty>=tstart) & (ty<=tend)
plt.figure(dpi=300,figsize=(11,3))
plt.plot(ty[tidx], 
         
         
         
         
plt.xlabel(
plt.ylabel(
plt.legend(
plt.grid(


