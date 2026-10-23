import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0
from numpy import sinc
from scipy.integrate import quad

# Part I (1): Kaiser parameter a vs 10~90% risetime
F = 10  # Two times of truncation frequency. Could be any number for Part I.
a_values = [0, 1, 2, 3, 4]  # Kaiser parameter alpha
t = np.linspace(-3, 3, int(1e5))  # Time vector. Remember to include negative time.
dt = t[1] - t[0]  # Time step

def w(t, F, a):
    return (F / (i0(np.pi *a)) * np.real(np.sinc(np.sqrt((F * t)**2 - (a)**2, dtype=complex))))

uW = np.zeros((t.size, len(a_values)), dtype=float)

# Compute wt and uW for each alpha
for i, alpha in enumerate(a_values):
    wt = w(t, F, alpha)
    uW[:, i] = (np.cumsum(wt) * dt)

# Plot Part I (1) results
plt.figure(figsize=(12, 6))
for i, alpha in enumerate(a_values):
    plt.plot(F * t, uW[:, i], label=f'$\\alpha$={alpha}')

plt.xlabel('Time')
plt.xlim(-5, 5)
plt.ylabel('$u_{w}(t)$')
plt.legend()
plt.title('Part I (1): Kaiser parameter a vs 10~90% risetime')
plt.grid(True)
plt.show()