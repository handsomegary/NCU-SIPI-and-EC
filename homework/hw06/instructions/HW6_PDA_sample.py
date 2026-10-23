#%% HW6 PDA

import 
import 
from 
from 
from 
from 
from 

# define pda function
def pda(sbr,NpUI,tol=0.0002):
    # sbr: single bit response (pulse response)
    # NpUI: number of time points per unit interval (integer)
    # tol: abs(voltage) smaller than tol will not be counted
    # It is assumed that the time vector for sbr is in uniform step!
    
    # Note: since NpUI is required to be an integer, we essentially
    # requires that dt (time step in t) divides Tb (unit interval).
    
    import numpy as np
    kmain = np.argmax(sbr) # main cursor center
    kstart = kmain - NpUI
    wc0 = np.zeros(                       # worst-case 0
    wc1 = sbr[kstart:                     # worst-case 1
    
    for k in range(2*NpUI):
        # pre-cursors
        precursors = sbr[kstart+k-NpUI::(-NpUI)]
        wc1[k] += # hint: use np.sum and Boolean indexing (conditional slicing)
        wc0[k] += 
        # wc1 post-cursors
        postcursors = 
        wc1[k] += 
        wc0[k] += 
        
    return wc1, wc0, kstart

#%% read S21 and compute single-bit response (sbr)

f, sp, _ = snpreader(
# remember to set the parameter 'maxsamplenum' for snpreader
S21 = 

# user defined parameters
Tb = 
NpUI = 
NpTr = 
Vpp = 

# dependent parameters
dt = 
df =    # frequency step in f
M =     # total number of points for IDFT
Tr =    # actual risetime in seconds

# spectrum of input trpezoidal
Vin = 

# Kaiser window. The parameter 'alpha' needs not be large 
# because Vin is already a lowpass response. If alpha is
# too large, the results will be further distorted.
a = 
W = 
W[0]    # halve the DC point in advance

sbr =   # IDFT of windowed Vout
tsbr =  # time vector to which sbr corresponds

#%% plot the sbr to determine a suitable truncation point


#%% truncate sbr
T2 = 10e-9       # truncation length (in seconds) of sbr
sbr = sbr[tsbr <= T2] # truncate
tsbr = tsbr[

#%% worst case eye
    = pda(

# normalized time vector for wc1 and wc0
twc = 

plt.figure(
plt.plot(twc,wc1,
plt.plot(

plt.
plt.
plt.


#%% Step 2: read input bit pattern
with open(
    bitpat = 

Nbp =  # number of total bits in bp
# convert string to bool
bp = np.array([int(x) for x in bitpat],dtype=bool)
del bitpat

#%% set parameters for sectional convolution
# Nbp = 2**16 
# Note: for bit pattern file, you already set Nbp above.
# For random NRZ, set an Nbp that is roughly the same as
# the Nbp for the bit file.


Nbps =  # number of bits per section

# dependent parameters
Ns =    # number of sections
NpT2 =  # number of points for duration T2
Npos =  # number of points per output section
rising_ramp = 
falling_ramp = 

#%% sectional convolution (almost the same as HW3)

tic = time()
y = np.zeros(Nbps*(Ns-1)*NpUI+Npos,dtype='float32') # total output waveform
ty =      # time vector for y

print('Sectional convolution completed %: ',end='')

for k in range(Ns): # for each section
    if (k%(Ns//10)) == 0:
        print(
    
    # construct input waveform of one section

    if k == 0:
        previous_bit = 
    else:
        previous_bit = 
    
    # set bpk (only activate one case below)
    # case 1: purely random
    bpk = np.random.  # hint: random integer array
    
    # case 2: from bit file
    bpk =           # the bit patterns of section k
    
    
    xk = 
    
    # No hint here. Same as HW 3
    
    
    
    y[

print('Total convolution time: %d sec'%(time()-tic))




#%% single-color eye
# Note: for this part, remember to plot the 1-UI eye first, and then
# duplicate the figure to make it 2-UI wide.

plt.figure(
# plot 1-UI single-colored eye
plt.plot(
# plot 1-UI single-colored eye and put it at the 2nd UI (i.e., shift the x axis)
plt.plot(

# plot PDA's worst-case boundaries
plt.plot(twc,wc1,
plt.plot(twc,

plt.xlim(
plt.ylim(
plt.title(
plt.xlabel(
plt.ylabel(


#%% build fixed-rate stateye (almost the same as HW5)

vmin =  # min voltage of stat eye (must be smaller than np.min(vo))
vmax =  # max voltage of stat eye (must be greater than np.max(vo))
dv =    # size of each voltage bin
Nvolt = # number of Vbins

# build 2D histogram
# Note: build 1-UI eye, not 2-UI eye
eyehist = np.zeros((
for i in range(y.size):
    idxt = 
    idxv = 
    eyehist[

# vertical normalization
total_count = 

stateye =  # array broadcasting

# duplicate the 1-UI eye to make it 2-UI wide
stateye = np.hstack((stateye

#%% plot statistical eye against PDA
# no hint here. See HW5




