
'''
HW4 spectrum of digital signal sample code
'''

#%% HW4: spectrum of digital (NRZ) signals
import numpy as np
from 
from 
from 
from 
import 
def db(x,power: bool): return 20/(int(power)+1)*np.log10(abs(x)+1e-15)

#%% read bit pattern (same as HW 3)
with open(


Nbp =  # number of total bits in bp
bp =   # bit pattern in numpy bool array

#%% averaged periodogram

Tb =     # unit interval
NbT =    # number of bits per truncation time (T)
NpUI =   # number of points per unit interval (UI)
NpTr =   # number of points per risetime (Tr)
NpT =    # number of points per T (in terms of above parameters)
K =      # averaging factor
df =     # desired frequency step
dt =     # time step based on Tb and NpUI
Tr =     # actual rise time in terms of above parameters
M =      # total number of points for DFT
f =      # frequency vector after DFT
Vpp =    # peak-to-peak voltage
rising_ramp = Vpp/NpTr * np.arange(NpTr,dtype='float')
falling_ramp = 
Z0 =     # transmission line characteristic impedance

PSD =    # theoretical PSD of purely random NRZ signals

# window functiion: use whatever window you like
n = 

wn =     # window function
G =      # window gain


S = np.zeros(M,dtype='float') # PSD initialize

tic = time()
for k in range(
    
    # Step 1: construct time-domain waveform (same as HW 3)
    # No hint here! Look at HW3's code sample
    
    
    # Step 2: perform fft (check periodogram method's equations)
    S += 

# outside the for loop, multiply necessary coefficients
S *= 

print('Total computation time: %d sec'%(time()-tic)) # Done!


#%% plot results
fG = f*1e-9
plt.figure(
plt.plot(
plt.plot(
plt.xlim(
plt.grid(
plt.ylim(
plt.legend()
plt.xlabel(
plt.ylabel(



