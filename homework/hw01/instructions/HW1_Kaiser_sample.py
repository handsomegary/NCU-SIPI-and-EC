# HW 1 sample code
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0
from scipy.fft import fft
from numpy import sqrt, pi
def db(x): return 20*np.log10(abs(x))

# 1 plot of time domain
data = np.loadtxt('D:\\HW1_Kaiser.txt',skiprows=1)
tn = data[:,0] * 1e-9 # time vector
xn = data[:,1]        # signal vector
dt = tn[1] - tn[0]    # time step (uniform)

plt.figure(dpi=200,figsize=(6,3))
plt.plot(tn*1e9,xn)
plt.xlabel('Time (ns)')
plt.ylabel('Signal (V)')
plt.tight_layout()
plt.show()

# 2 apply Kaiser window and perform spectral analysis
N = tn.size
n = np.arange(N)
a = [0.1, 0.2, 1, 2, 3, 4, 5] # the parameter 'alpha'

for alpha in a:
    wn = i0(pi * alpha * sqrt(1 - ((2*n/(N-1)) - 1)**2)) / i0(pi * alpha) # define Kaiser window
    G = np.sum(wn) / N       # window gain. Hint: use np.sum()
    wnxn = xn*wn             # multiply window to signal
    M = 2**20-1                # total number after zero padding
    df = 1 / (M * dt)        # the DFT bin (frequency step)
    f = np.arange(0, M)*df   # frequency vector (for plotting)
    Xk = dt*fft(wnxn, n=M)/G # estimate the continuous FT
    A = np.abs(Xk)/N*2/dt    # estimate the amplitude of frequency components
    # %matplotlib qt
    # %matplotlib inline
    plt.figure(dpi=200,figsize=(6,3))
    plt.plot(f*1e-9,db(A))
    plt.xlim(0,2)
    plt.ylim(-150,5)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude (dB)')
    plt.grid(color='lightgray')
    plt.tight_layout()
    plt.show()

# 3 The estimated Q, Ai, and fi
A = 1
f1 = 525*1e6
f2 = 601*1e6
x = np.arange(0, 20*1e-9, 0.01*1e-9)
y = A*np.sin(2*np.pi*f1*x)+A*np.sin(2*np.pi*f2*x)
plt.figure(dpi=200,figsize=(6,3))
plt.xlim(0,20*1e-9)
plt.ylim(-2,2)
plt.plot(x,y)
plt.xlabel('Time (ns)')
plt.ylabel('Amplitude (V)')
plt.tight_layout()
plt.show()
