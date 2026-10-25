import numpy as np
import matplotlib.pyplot as plt
txt_array = np.loadtxt('last.txt')
x = txt_array[:, 0]  # First column
y = txt_array[:, 1]  # Second column

# Plot the data
plt.plot(x, y)
plt.xlabel('X-axis label')
plt.ylabel('Y-axis label')
plt.title('Plot title')
plt.grid(True)
plt.show()
