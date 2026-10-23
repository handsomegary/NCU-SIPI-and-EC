#%% HW 5: triggered eye
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt(
time = 
vo = 
clk = 


#%% single-color eye (just for reference, not part of homework)
plt.plot(np.mod(time*1e12,50),vo,'.',markersize=2)

#%% build fixed-rate stateye

# user defined parameters
UI =        # nominal unit interval (ps)
Ntime =     # number of time bins per UI (not per 2UI)
vmin =      # min voltage of stat eye (must be smaller than np.min(vo))
vmax =      # max voltage of stat eye (must be greater than np.max(vo))
dv =        # size of each voltage bin

# calculated parameters
dt =        # size of each time bin
Nvolt =     # number of voltage bins in stat eye

# build 2D histogram
eyehist = np.zeros((Nvolt,2*Ntime),dtype=int)
for i in range(time.size):
    idxt =    # horizontal index of time[i] in the 2-UI stateye
    idxv =    # vertical index of vo[i] in the stateye
    eyehist[idxv,idxt] += 

# vertical normalization
total_count = # total number of occurrence at each column of eyehist

stateye =     # final fixed-rate stat eye (use array broadcasting)


#%% plot statistical eye
plt.figure(dpi=300)
plt.imshow(stateye,cmap='??',origin='??',vmin=?,vmax=?,extent=[?,?,?,?])
plt.colorbar()
plt.xlabel(
plt.ylabel(
plt.title('Fixed Rate Statistical Eye (UI = %.1f ps)'



#%% build triggered stateye

# build 2D histogram
eyehist2 = np.zeros((
trigger_point = 0
for i in range(1,time.size):
    if     # trigger condition
        trigger_point = i
    if trigger_point == 0:
        continue
    
    idxt =    # horizontal index of time[i] in the 2-UI stateye
    if idxt >= (2*Ntime):
              # observation time is more than 2UI away from the previous trigger
    
    idxv =    # vertical index of vo[i] in the stateeye
    eyehist2[idxv,idxt] 

# vertical normalization
total_count =

stateye2 =    # final triggered stat eye


#%% plot statistical eye
plt.figure(
plt.imshow(
plt.colorbar()
plt.xlabel(
plt.ylabel(
plt.title(
