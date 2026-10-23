import numpy as np
import matplotlib.pyplot as plt
txt_array = np.loadtxt('HW3_Rx_waveform_ADS.txt',skiprows=1)
x = txt_array[:, 0]  # First column
y = txt_array[:, 1]  # Second column

# Plot the data
plt.plot(x, y)
dt = x[2] - x[1]
plt.xlim(4.058e-6, 4.068e-6)
plt.xlabel('X-axis label')
plt.ylabel('Y-axis label')
plt.title('Plot title')
plt.grid(True)
plt.show()
