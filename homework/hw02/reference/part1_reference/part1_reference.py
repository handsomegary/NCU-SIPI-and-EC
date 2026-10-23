import numpy as np
import matplotlib.pyplot as plt
from scipy.special import i0
from numpy import sinc
from scipy.integrate import quad

# Part I (1): Kaiser parameter a vs 10~90% risetime
F = 40e9  # Two times of truncation frequency. Could be any number for Part I.
a_values = [0, 1, 2, 3, 4]  # Kaiser parameter alpha
t = np.linspace(-3e-9, 3e-9, int(1e5))  # Time vector. Remember to include negative time.
dt = t[1] - t[0]  # Time step

def w(t, F, a):
    return (F / (i0(np.pi *a)) * np.real(np.sinc(np.sqrt((F * t)**2 - (a)**2, dtype=complex))))

uW = np.zeros((t.size, len(a_values)), dtype=float)

# Compute wt and uW for each alpha
for i, alpha in enumerate(a_values):
    wt = w(t, F, alpha)
    uW[:, i] = (np.cumsum(wt) * dt)

# Plot Part I (1) results
plt.rcParams["font.family"] = "Times New Roman"
plt.figure(dpi=500,figsize=(5, 4))
for i, alpha in enumerate(a_values):
    plt.plot(F * t, uW[:, i], label=f'$\\alpha$={alpha}')

plt.xlabel('Normalized Time ($F \cdot t$)', fontweight='bold')
plt.xlim(-5, 5)
plt.ylabel('Kaiser Steps', fontweight='bold')
plt.legend()
plt.grid(True)
plt.show()

# Part 2: Plot the normalized 10~90% risetime of uW(t) with respect to α
risetime_values = []

for i, alpha in enumerate(a_values):
    t_10 = np.where(uW[:, i] >= 0.1)[0][0]
    t_90 = np.where(uW[:, i] >= 0.9)[0][0]
    risetime = (t[t_90] - t[t_10]) * F
    risetime_values.append(risetime)

# Plot Part 2 results
plt.rcParams["font.family"] = "Times New Roman"
plt.figure(dpi=500,figsize=(5, 3))
plt.plot(a_values, risetime_values, '--o', markerfacecolor='none', color = (0.8, 0, 0))
plt.xlabel('Kaiser parameter $\\alpha$', fontweight='bold')
plt.ylabel('Normalized 10~90% Risetime $F \cdot T_{r}$', fontweight='bold')
plt.ylim(0, 3)
plt.grid(True)
plt.show()