# %% Part II: IDFT with Kaiser
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0
from snpreader import snpreader
from scipy.fft import ifft

def db(x):
    return 20 * np.log10(abs(x))

# load S parameter
freq, sp, Z0 = snpreader('HW2_Chebyshev_LPF.s2p', reciproc=True)

plt.figure(dpi=500,figsize=(6, 4))
plt.rcParams["font.family"] = "Arial"
plt.plot(freq * 1e-9, db(sp[:, 0]), color = (0, 0, 0.8), label='S11', linewidth = 2)
plt.plot(freq * 1e-9, db(sp[:, 1]), color = (1, 0, 0), label='S21', linewidth = 2)
plt.xlabel('Frequency (GHz)')  # Fill in the xlabel
plt.ylabel('Magnitude (dB)')  # Fill in the ylabel
plt.grid(True)  # Fill in grid parameters if needed
plt.legend()
plt.ylim(-90, 10)
plt.show()

F = 40e9  # Replace with your desired value
#%% multiply with the window function
S11 = sp[:,0]   # extract S11
df = freq[1] - freq[0]    # frequency step (calculated from 'freq')
W = lambda f: i0(np.pi * a * np.sqrt(1 - (2 * f / F)**2)) / i0(np.pi * a)
t1 = 25e-12  # addition time delay 
a = 1   # Kaiser parameter
WS11 = W(freq) * S11 * np.exp(-1j*2*np.pi*freq*t1)# W(f) * S11
WS11[0] = 0 # special treatment for the DC point

# IDFT
M =  2**20  # total number of points for IDFT
WS11_padded = np.pad(WS11, (0, M - WS11.size), 'constant')
hW = 2 * df * M * np.real(ifft((WS11_padded)))  # windowed impulse response
dt = 1/ (df * M)
sW =  np.cumsum(hW) * dt   # windowed step response

t = np.arange(0, sW.size) * dt

txt_array = np.loadtxt('ngspice.txt')
x = txt_array[:, 0]  # First column
y = txt_array[:, 1]  # Second column


plt.figure(dpi=500,figsize=(6, 4))
plt.rcParams["font.family"] = "Arial"
plt.plot(x * 1e9, y, color = (0, 0, 0.8), label='Ngspice')
plt.plot(t * 1e9, sW.real, color = (1, 0, 0), label='Kaiser IDFT')
plt.xlabel('Time (ns)')
plt.ylabel('Step response (V)')
plt.xlim(0, 3)
plt.grid(True)
plt.legend(prop={'weight': 'bold'})
plt.show()

