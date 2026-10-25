

#%% HW9: DFE simulation

# copy necessary functions from previous homework

#%% define function: DFE simulator

def fixed_rate_DFE(y,NpUI,c,Vt,sample_point):
    '''
    Simulate a decision feedback equalizer assuming input has fixed data rate.
    
    The slicer output of the DFE is assumed to be [0,1] (i.e., single-ended).
    
    The user has to supply a sample_point for this function, which is the point
    where the DFE makes decision (logic 0 or 1). The DFE will make decisions at
    points of the form (sample_point + i*NpUI) for all possible i (>0 and <0).
    
    This code assumes the sample_point is located at the center of the eye. In
    other words, the clock trigger points are assumed to be at (sample_point 
    + NpUI//2 + i*NpUI) for all possible i (>0 and <0).
    
    Parameters
    ----------
    y : 1-dim np array
        The input waveform at the front of DFE. Should be long enough.
        Should have uniform time step.
    NpUI : int
        Number of points per UI
    c : DFE tap values
        Python list or 1-dim np array
    Vt : float
        threshold voltage for decision making
    sample_point : int
        The sampling point at which the DFE makes decisions.
        The DFE will make decisions at (sample_point + i*NpUI).

    Returns
    -------
    z : 1-dim np array
        Post-DFE waveform, i.e., waveform at the summer output of the DFE.

    '''
    import numpy as np
    
    sample_start = sample_point % NpUI
    
    print('Threshold voltage Vt = %.3f V'%Vt)
    print('Sampling point relative to time axis = %d %%'%(sample_start/NpUI*100))
    
    z = y.copy() # output waveform
    
    for k in range(sample_start,?,?)
        # in each iteration, there are two things to do:
        # 1st: make decision for the current bit
        # 2nd: update the future waveform of z based on the current decision
        
        decision = 
        
        if decision:
            for i in 
                z[k+NpUI//2+i*NpUI : ?? ] += 
    
    
    return z





#%% read S2P and plot impulse response

import 
import 
from 
def 

f, sp, _ = snpreader(

H = # channel transfer function

plt.


df = 

# check impulse response
t, h = imp_resp(

#%% set parameters
T2 =  # truncation length (in seconds) of impulse response
Tb = 
Tr = 
NpUI = 
Vpp = 


#%% get SBR for DFE design

tsbr, sbr = get_sbr(
plt.

#%% get cursor values for DFE design

Npre = 
Npost = 
a = get_cursor(



#%% DFE design
# the three parameters below are the most important ones for a single-ended DFE
c =    # DFE taps. You can specify any number you wish, or based on 'a' (cursor value)
Vt =   # threshold voltage for DFE decision making. Think: how to choose it?
sample_pt =    # suggested sampling point. Think: how to choose it?


#%% compute time domain waveform

# Step 1: sectional convolution to construct PRBS response at Rx
# Below We use the function in HW 8, but please do not enable the FFE option.
# Alternatively, you can also use the sctn_conv_random_NRZ function of HW 7.
? = sctn_conv_random_NRZ_TxFFE(

# Step 2: compute dsum waveform of DFE
z = fixed_rate_DFE(

# build statistical eye
? = fixed_rate_unifdt_stateye(


#%% plot stateye

# No hint here. See previous homework

#%% plot time-domain waveform

# No hint here. See previous homework



