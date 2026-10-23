# %% Part II: IDFT with Kaiser
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0
from numpy import sqrt, pi, sinc
import snpreader
from scipy.fft import ifft

def db(x):
    return 20 * np.log10(abs(x))

# load S parameter
freq, sp, Z0 = snpreader.snpreader('HW2_Chebyshev_LPF.s2p', reciproc=True)

plt.figure()
plt.plot(freq * 1e-9, db(sp[:, 0]), label='S11')
#plt.plot(freq * 1e-9, db(sp[:, 1]), label='S21')  # Fill in S21 data if needed
plt.xlabel('Frequency (GHz)')  # Fill in the xlabel
plt.ylabel('Magnitude (dB)')  # Fill in the ylabel
plt.grid(True)  # Fill in grid parameters if needed
plt.legend()
plt.show()

plt.figure()
plt.plot(freq * 1e-9, db(sp[:, 1]), label='S21')
#plt.plot(freq * 1e-9, db(sp[:, 1]), label='S21')  # Fill in S21 data if needed
plt.xlabel('Frequency (GHz)')  # Fill in the xlabel
plt.ylabel('Magnitude (dB)')  # Fill in the ylabel
plt.grid(True)  # Fill in grid parameters if needed
plt.legend()
plt.show()
F = 40e9  # Replace with your desired value
#%% multiply with the window function
S11 = sp[:,0]   # extract S11
df = freq[1] - freq[0]    # frequency step (calculated from 'freq')
W = lambda f: i0(np.pi * a * np.sqrt(1 - (2 * f / F)**2)) / i0(np.pi * a)
t1 =  250e-15   # addition time delay
a =   2   # Kaiser parameter
WS11 = W(freq) * S11  # W(f) * S11
WS11[0] = 0 # special treatment for the DC point

# IDFT
M =  2**16  # total number of points for IDFT
WS11_padded = np.pad(WS11, (0, M - WS11.size), 'constant')
hW = 2 * df * M * np.real(ifft((WS11_padded)))  # windowed impulse response
dt = 1/ (df * M)
sW =  np.cumsum(hW) * dt   # windowed step response

t = np.arange(0, sW.size) * dt
# plot results
plt.figure()
plt.plot(t, sW.real, label='Windowed Step Response')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.xlim(0, 3e-9)
plt.grid(True)
plt.legend()
plt.show()
