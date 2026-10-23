#%% HW2: Kaiser and IDFT
#%% Part I (1): Kaiser parameter a vs 10~90% risetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0
from numpy import sin, sinh, sqrt, pi, sinc
def db(x): return 20*np.log10(abs(x))

F =      # Two times of truncation frequency. Could be any number for Part I.
a =      # Kaiser parameter alpha
t =      # Time vector. Remember to include negative time.
dt =     # Time step
wt = np.zeros((t.size,len(a)),dtype=float)
for i in range(len(a)):
    wt[:,i] =                # the time-domain waveform of Kaiser window
uW =         # integrate wt to obtain Kaiser step. You can use np.cumsum.
   # Remember to multiply the results by dt to obtain the correct values.

plt.figure()
plt.plot(
plt.ylim(
plt.xlim(
plt.grid(
plt.xlabel(
plt.ylabel(
plt.legend([r'$\alpha$ = %d'%a[i] for i in a],fontsize=12)

#%% Part I (2): find 10~90% risetime
Tr1090 = np.zeros((len(a)))
for i in range(len(a)):
    t10 = np.where(     # find the index where uW[:,i] crosses the 10% point
    t90 =
    Tr1090[i] = 

plt.figure(
plt.plot(
plt.grid(
plt.xlabel(
plt.ylabel(
plt.ylim(

#%% Part II: IDFT with Kaiser
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0
from numpy import sqrt, pi, sinc
from snpreader import snpreader
from scipy.fft import ifft
def db(x): return 20*np.log10(abs(x))

# load reference answer from Ngspice
data_ng = np.loadtxt(
tng =     # time vector
vng =     # voltage vector
dt =      # time step

# load S parameter
freq,sp,Z0 = snpreader('HW2_Chebyshev_LPF.s2p')

plt.figure(
plt.plot(freq*1e-9,db(sp[:,0]),label='S11')
plt.plot(
plt.xlabel(
plt.ylabel(
plt.grid(
plt.legend(
plt.ylim(
#%% multiply with the window function
S11 = sp[   # extract S11
df =     # frequency step (calculated from 'freq')

t1 =     # addition time delay
a =      # Kaiser parameter
WS11 =   # W(f) * S11
WS11[0]  # special treatment for the DC point

# IDFT
M =      # total number of points for IDFT
hW =     # windowed impulse response
sW =     # windowed step response


# plot results
plt.figure(
plt.plot(tng*1e9,
plt.plot(
plt.xlim(
plt.xlabel(
plt.ylabel(
plt.grid(
plt.legend(
