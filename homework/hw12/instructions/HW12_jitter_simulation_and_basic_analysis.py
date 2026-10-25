

#%% HW 12: jitter simulation and basic analysis

#%% define functions

# Copy necessary functions from previous homework to here

# In particular, in the StatEye construction function, do 
# not add 1e-15 to the eye array, because that will affect
# the subsequent jitter analysis. Instead, you can add
# 1e-15 to the array when plotting.



#%% define functions: sectional convolution with Tx FFE, for both SE and DP

def sctn_conv_random_NRZ_TxFFEjitter_SEDP(H, df, T2, Tb, Tr, Td, NpUI, 
            cFFE=[1],Vpp=1,SJamp=0,SJfreq=1e6,RJsigma=0,
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
    T2: truncation time of channel impulse response in seconds
    Tb: unit interval in seconds
    Tr: risefall time in seconds
    Td: initial time delay
    NpUI: number of points per unit interval
    Vpp: peak-to-peak voltage of NRZ signal. The post-FFE waveform
         at Tx still has peak-to-peak voltage of Vpp. It may not be
         in the range [0,Vpp]. Instead, it may be in the range of
         [-q,Vpp-q] for some positive number q, depending on cFFE.
         
         Important: for differential signaling (sig='DP'), this Vpp
         is the peak-to-peak voltage of each single-ended driver. 
         Therefore, the actual p2p voltage of the differential signal
         at Tx will be 2*Vpp.
    
    cFFE: FFE taps, list or numpy 1-dim array,  [c0, c1, c2, ...]
       Note: c0 is not necessarily the main cursor tap. Instead,
       c0 is the "left-most" cursor tap. It may be the main cursor
       or it may be a precursor tap. If it is the 1st precursor, 
       than c1 is the main cursor, c2 is the 1st post cursor, and
       so on. If cFFE=[1], then there is no FFE. The absolute values of
       the elements in cFFE are not important; they will be normalized
       to meet the Vpp constraint. Only the relative ratios in cFFE are
       relevant.
    
    SJamp: sinusoidal jitter amplitude (in seconds)
    SJfreq: sinusoidal jitter frequency (in Hertz)
    RJsigma: random jitter standard deviation (in seconds)
        RJ will follow the Gaussian distribution
    
    Nbp: minimum number of random bits to simulate
    sig: signaling scheme, a string, must be either 'SE' or 'DP'.
    a: Kaiser window parameter
    seed: seed value for random number generator

    Returns
    -------
    ty: time vector of y
    y: total output waveform
    tx: time vector of x
    x: total input waveform
    bp: the overall input bit pattern for validation if needed
        It is always in {0,1}, regardless of 'sig' (not {1,-1}).

    '''
    
    # import necessary modules
    import 
    from 
    from 
    from 
    from 
    
    rng = 
    
    f =  # frequency vector of H
    dt = # output time step
    M =  # total number of points for DFT
    NpT2 =  # number of points for duration T2
    c = 1 - NpT2/M
    if c < 0:
        raise ValueError('T2 too long. Use smaller dt or df (interpolation).')
    Nbps = round(0.8*c/df/Tb) # number of bits per section
    Ns =  # number of sections
    NpTr =  # number of points per risetime
    NpTd =  # number of points for delay Td
    Npos =  # number of points per output section
    unit_ramp = 
    
    print('(M,Nbps,Ns,Npos) = (%d,%d,%d,%d)'%(M,Nbps,Ns,Npos))
    
    y = np.zeros( # total output waveform
    # y will be truncated to length (Ns-1)*Nbps*NpUI after convolution
    ty =  # time vector for this y
    
    W =  # Kaiser window
    W[0
    
    # jitter processing
    wSJ =  # angular frequency of SJ
    SJ = 
    # array of SJ at the ideal transition points, size == Ns*Nbps
    
    RJ = rng.
    # array of RJ, size == Ns*Nbps
    
    TJ =  # TJ value at all transition times, size = Ns*Nbps
    
    # === hard limit the TJ to abs(Tb-Tr) === #
    TJ[TJ > (Tb-Tr)] = (Tb-Tr)
    TJ[TJ <-(Tb-Tr)] = 
    NpTJ =    # express TJ as number of points
    
    
    # === FFE processing === #
    cFFE =  # in case cFFE is a Python list, not ndarray
    cFFE =  # normalization
    
    # Construct the complete input bit pattern in {0,1}.
    # It is the one that will be returned.
    bp = rng.
    
    # Construct the post-FFE input waveform bpFFE
    # The effect of risetime will be added to x in the for loop
    # The effect of Vpp is added to bpFFE now
    if sig=='SE':
        # use Vpp*bp for convolution with cFFE
        bpFFE = 
    elif sig=='DP':
        # First convert the {0,1} NRZ signal to {-1,1} form then convolution
        # The overall p2p voltage of the resulting diff signal becomes 2Vpp.
        bpFFE = 
    else:
        raise ValueError('Unrecognized signaling scheme.')
    
    tic = time()
    
    
    # build x, the total input waveform, including jitter
    x = np.zeros(
    now = NpTd
    if bpFFE[0] != 0:
        if now + NpTJ[0] > 0:
            x[now+NpTJ[0]:??] = 
            x[?? : ??] = 
        else: # start transition at t = 0, not negative time
            x[ = 
            x[ = 
    now += 
    for j in range(
        if bpFFE[j] ==  # no transition, no jitter
            x[
        else:
            if NpTJ[j] > 0: # extend the waveform of previous bit
                x[now:now+NpTJ[j]] = 
            x[ = 
            x[ = 
        now 
    tx = 
    
    
    print('Sectional convolution progress:     ',end='')
    
    for k in range(Ns): # for each section
        if abs((k%((Ns-1)/10))) < 1:
            # print(k,Ns,Ns/10)
            print('\b\b\b\b%3d%%'%(round(10*k/((Ns-1)/10))),end='')
        
        xk = 
        
        # convolution of one section
        y[
    
    y = y[ # truncate y
    ty = ty[ # truncate ty
    # This is why the 'Ns' parameter was set to Nbp//Nbps+2. 
    # We add two more sections, then truncated.
    
    print('\nTotal convolution time: %d sec'%(time()-tic))

    # print(NpTJ[0:30])
    return ty, y, tx, x, bp
    # return xk and bp for checking



#%% HW12: jitter simulation
import 
import 
from 
from 
def 

 = snpreader(  # be careful of maxsamplenum

df = 

#%% check S parameters
plt.figure(
plt.plot(
plt.ylim(
plt.xlim(

#%% define parameters
# pretty much same as HW11










SJamp = 3e-12
SJfreq = 1e6
RJsigma = 1e-12

dv = 1e-3
dt = Tb / NpUI


#%% solve MNA to get transfer function (H)
# same as HW11














#%% compare H with the original S21
plt.figure(dpi=250)
plt.plot(freq*1e-9,db(sp[:,1]),label='S21')
plt.plot(freq*1e-9,db(H),label='H')
plt.ylim(
plt.xlim(
plt.xlabel(
plt.ylabel(
plt.legend(

#%% compute impulse response and find truncation time (T2)

t, h = 

#%% check impulse response to determine truncation time
# %matplotlib qt
plt.plot(t, h)
plt.xlim(
      
#%% set truncation time
T2 = 

#%% sectional convolution to obtain output responses

# Case 1: no Tx PJ & RJ
ty, y1, tx, x1 = sctn_conv_random_NRZ_TxFFEjitter_SEDP(
    
    )[:4]

# Case 2: Tx PJ on, no RJ
_, y2, _, x2 = 

# Case 3: Tx RJ on, no PJ
_, y3, _, x3 = 

# Case 4: Tx PJ & RJ on
_, y4, _, x4 = 


#%% build statistical eye

eyes_x1 = fixed_rate_unifdt_stateye(x1,NpUI,dv) # (stateye, extent)
eyes_y1 = 

eyes_x2 = 
eyes_y2 = 

eyes_x3 = 
eyes_y3 = 

eyes_x4 = 
eyes_y4 = 

#%% plot stateye

%matplotlib inline
fig, axs = plt.subplots(1,1,figsize=(7,4),dpi=300)
# im = axs.imshow(np.log10(eyes_x[0]),interpolation='None',cmap='gist_ncar_r',aspect='auto',origin='lower',vmax=-1,extent=eyes_x[1])
im = axs.imshow(np.log10(eyes_y4[0]+1e-30),cmap='gist_ncar_r',aspect='auto',origin='lower',
                vmax=-1,vmin=-15,extent=eyes_y4[1])
fig.colorbar(im,label=r'$\mathrm{log}_{10}$(probability)')
axs.set_xlabel('Time (UI)')
axs.set_ylabel('Voltage (V)')
axs.set_title('Fixed Rate Statistical Eye (UI = %.1f ps)'%(Tb*1e12))
axs.grid(color='lightgray')
axs.set_ylim(-0.2,1)



#%% plot time-domain waveform

tstart = 145e-9
tend = tstart + 4e-9
tidx = (ty>=tstart) & (ty<=tend)
plt.figure(dpi=300,figsize=(11,3))
plt.plot(ty[tidx], y[tidx],'b',linewidth=2,label='Rx waveform')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.legend(loc='upper right',ncol=1)
plt.grid(color='lightgray')


#%% check input waveform

tstart = 1000e-9
tend = tstart + 4e-9
tidx = (tx>=tstart) & (tx<=tend)
plt.figure(dpi=300,figsize=(11,3))
plt.plot(tx[tidx]-tstart, x[tidx],'b',linewidth=2,label='Tx waveform')
plt.xlabel('Time (s)')
plt.ylabel('Voltage (V)')
plt.legend(loc='upper right',ncol=1)
plt.grid(color='lightgray')



#%% jitter histogram
# For each case, you need to take out the right slice from
# stateye to get jitter distribution. Be careful that the 
# vertical axis of the stateye of each case is different.

# Tx jitter
Vtx =   # threshold voltage at Tx (hint: depend on Vpp)

# Case 1
vmin = eyes_x1[1][2]
Vtidx =    # vertical index of Vtx at the Case 1 eye
jtx1 = eyes_x1[0][Vtidx,

# Case 2
vmin = 
Vtidx = 
jtx2 = 

# Case 3
vmin = 
Vtidx = 
jtx3 = 

# Case 4
vmin = 
Vtidx = 
jtx4 = 


# Rx jitter
Vty =  # threshold voltage at Rx (hint: Vpp and H[0])

# Case 1
vmin = 
Vtidx = 
jty1 = 

# Case 2
vmin = 
Vtidx = 
jty2 = 

# Case 3
vmin = 
Vtidx = 
jty3 = 

# Case 4
vmin = 
Vtidx = 
jty4 = 



# note: these are not jitter PDF (PMF) yet
plt.figure(dpi=250)
plt.semilogy(jtx1,label='jtx1')
plt.semilogy(jtx2,label='jtx2')
plt.semilogy(jtx3,label='jtx3')
plt.semilogy(jtx4,label='jtx4')
plt.semilogy(jty1,label='jty1')
plt.semilogy(jty2,label='jty2')
plt.semilogy(jty3,label='jty3')
plt.semilogy(jty4,label='jty4')
plt.legend()

#%% take out only 1 UI of jitter
jtx1 = jtx1[?:?+NpUI]
jtx2 = jtx2[
jtx3 = jtx3[
jtx4 = jtx4[

jty1 = jty1[
jty2 = jty2[
jty3 = jty3[
jty4 = jty4[

# note: these are not jitter PDF (PMF) yet

plt.figure(dpi=250)
plt.semilogy(jtx1,label='jtx1')







plt.legend()


#%% normalize the jitter histogram to make it sum to 1
jtx1 = jtx1 / 
jtx2 = 
jtx3 = 
jtx4 = 
jty1 = 
jty2 = 
jty3 = 
jty4 = 

# Note: these are jitter PMF (probability mass function),
# not PDF (probability density function), because the value
# at each cell is "probability" of that time bin, not the
# "probability density". The PMF is convenient for computing
# the jitter mean, variance, and rms. However, keep in mind
# that in the next lecture we are going to fit the jitter
# PDF by the Gaussian PDF. There, we need to convert PMF
# to PDF. (Search online for the difference of PMF and PDF.)

# time vector for jitter PMF, in ps
tjp = np.linspace(?,?,?) * 1e12

# note: these are now the jitter PMF
plt.figure(dpi=250)
plt.semilogy(tjp,jtx1,







plt.xlabel
plt.ylabel(
plt.legend(


#%% compute jitter statistics
idx = np.arange(NpUI)

# Compute the mean of each PMF. Unit: none (index)
# Hints:
#   (1) A PMF is already normalized, so no need to normalize again
#   (2) Search online for the formula of mean given its PMF
#   (3) A simple np.sum function suffices
#   (4) The unit of these means are None, i.e., they are just indices.
#       These means are just for the computation of variance, so no need
#       to convert to the unit of time.
jtx1mean =
jtx2mean = 
jtx3mean = 
jtx4mean = 

# Compute the variance of each PMF.
# Search online for the formula of variance given its PMF
# The unit should be ps**2 or s**2, so dt is needed here.
# Also, the jtxmean computed above is also needed.
jtx1var = 
jtx2var = 
jtx3var = 
jtx4var = 


# Compute the same parameters at Rx
jty1mean = 
jty2mean = 
jty3mean = 
jty4mean = 
jty1var = 
jty2var = 
jty3var = 
jty4var = 


print('jtx1var =',jtx1var)
print('jtx2var =',
print('
print('

print('jty1var =',
print('jty2var =',
print('
print('

# print rms values
# Because should be centered, meaning that the mean is zero,
# so RMS value equals the square root of variance, i.e., rms
# value is just the standard deviation.
print('jtx1rms =',
print('jtx2rms =',
print('
print('

print('jty1rms =',
print('jty2rms =',
print('
print('

