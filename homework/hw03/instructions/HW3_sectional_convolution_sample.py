
'''
HW3 sectional convolution sample code
'''

#%% Import necessary modules
import numpy as np
from 
from 
from 
from 
import 
from time import time
def db(x): return 20*np.log10(abs(x)+1e-15)

#%% Step 1: estimate length of impulse response
freq, sp, _ = snpreader(
S21 = sp[

dt =            # choose the dt you want
df =            # what is the df?

# multiply with the window function
a =             # Kaiser parameter
WS21 = S21 /    # S21 multiply with Kaiser 
WS21[0] /= 2    # divide DC point by 2

# IDFT
M =             # total number of points for IDFT
t =             # time vector
vimp =          # IDFT to get (windowed) impulse response

plt.plot(t*1e9,db(vimp))
plt.xlim(       # adjust the xlim to find out appropriate truncation time
plt.grid()

#%% Set T2
T2 = 

#%% Step 2: read input bit pattern
with open('HW3_bitpattern_1048576.txt','r') as f:
    bitpat = f.read()  # note: this is string

Nbp =            # number of total bits in bitpat

# convert string to bool array
bp = np.array([int(x) for x in bitpat],dtype=bool) # use dtype=bool to save memory

del bitpat

#%% Step 3: set parameters
# Some of the parameters below are free for you to define, while some are 
# related to others (e.g., dt is related to NpUI and Tb) and you need to 
# type in the correct expression rather than value.

BR = 10e9 # bit rate
Tb =      # unit interval = symbol period
Nbps =    # number of bits per section
Ns = Nbp//Nbps # number of sections
NpUI =    # number of points per UI
dt =      # sampling time, depend on NpUI (different from the dt in Step 1)
NpTr =    # number of points per risetime
NpT2 =    # number of points for duration T2
Npos =    # number of points per output section
Vpp = 1
rising_ramp = Vpp/NpTr * np.arange(NpTr,dtype='float32') 
          # rising waveform for constructing input waveform
          # use dtype='float32' instead float to save memory
falling_ramp = 
          # falling waveform for constructing input waveform

M =       # total number of points for DFT and also IDFT

tx = dt * np.arange(NpUI*Nbps,dtype='float32') 
          # time vector for one section of input waveform


a =       # Kaiser window parameter
W =       # Kaiser window for IDFT
W[0] /= 2 # halve the DC point in advance, so that we don't 
          # have to do this later in the for loop


#%% Step 4: sectional convolution

tic = time()

y = np.zeros(??,dtype='float32') # "total" output waveform
ty =                             # time vector for y

print('Sectional convolution completed %: ',end='')

for k in range(Ns): # for each section
    if (k%(Ns//10)) == 0:
        print('%d '%(10*k/(Ns//10)),end='')
    
    # Step 4-1: construct input waveform of one section

    if k == 0:
        previous_bit = 0 
    else:
        previous_bit =   # the last bit of the previous section
    bpk = bp[k*Nbps : ??] # the bit patterns of section k
    xk = np.zeros(NpUI*Nbps,dtype='float32') # input waveform
    
    if previous_bit:
        if bpk[0]: # 1-1, no transition
            xk[:NpUI] = Vpp
        else: # 1-0 transition
            xk[:NpTr] = 
    else: # previous_bit == 0
        if bpk[0]: # 0-1 transition
            xk[:NpTr] = 
            xk[NpTr:NpUI] = Vpp
        else: # 0-0, no transition
            pass # no action needed
    
    for j in range(1,Nbps):
        if bpk[j-1]:
            if bpk[j]: # 1-1, no transition
                xk[j*NpUI:??] = 
            else: # 1-0 transition
                xk[j*NpUI:??] = 
        else: # previous_bit == 0
            if bpk[j]: # 0-1 transition
                xk[ = 
                xk[??:??] = Vpp
            else: # 0-0, no transition
                
        
    # Step 4-2: convolution of one section
    
    y[k*Nbps*NpUI : ?? ] += 
        # complete all operations for each loop of sectional
        # convolution in just one line, including FFT, IFFT,
        # zero-padding, windowing, etc!

print('Total convolution time: %d sec'%(time()-tic))

#%% load ADS's waveform
data = np.loadtxt(
tADS = data[
yADS = data[

#%% plot results
tstart = 
tend = tstart + 10e-9
plt.figure(
plt.plot(ty
plt.plot(tADS
plt.xlabel(
plt.ylabel(
plt.legend(

