"""Independent numerical reference models for HW11 verification.

These helpers use terminal-wave boundary conditions and direct convolution,
not the production modified nodal analysis or FFT overlap-add algorithms.
"""

from __future__ import annotations

import numpy as np


def wave_terminated_transfer(frequency, scattering, z0, rs=30.0, cs=1.5e-12,
                             rl=65.0, cl=1e-12):
    """Voltage transfer from unit Thevenin sources using reflected waves.

    For voltage waves v=a+b and i=(a-b)/z0, termination KCL gives
    a=q+gamma*b. Network scattering gives b=S*a. Thus
    (I-diag(gamma)S)a=q and v=(I+S)a.
    """
    frequency = np.asarray(frequency)
    scattering = np.asarray(scattering)
    nports = scattering.shape[1]
    tx = np.arange(0, nports, 2)
    rx = tx + 1
    result = np.empty((len(frequency), len(rx), len(tx)), complex)
    identity = np.eye(nports)
    source = np.zeros((nports, len(tx)))
    source[tx, np.arange(len(tx))] = 1.0 / rs
    for k, f in enumerate(frequency):
        y = np.empty(nports, complex)
        y[tx] = 1.0 / rs + 2j * np.pi * f * cs
        y[rx] = 1.0 / rl + 2j * np.pi * f * cl
        denominator = 1.0 + z0 * y
        gamma = (1.0 - z0 * y) / denominator
        q = z0 * source / denominator[:, None]
        a = np.linalg.solve(identity - gamma[:, None] * scattering[k], q)
        result[k] = ((identity + scattering[k]) @ a)[rx]
    return result


def ideal_through_scattering(nchannels=1):
    """Zero-length uncoupled, lossless, reflectionless connections."""
    s = np.zeros((nchannels * 2, nchannels * 2))
    for j in range(nchannels):
        s[2*j, 2*j+1] = 1
        s[2*j+1, 2*j] = 1
    return s


def ideal_through_rc_transfer(frequency, rs=30.0, cs=1.5e-12,
                              rl=65.0, cl=1e-12):
    """Closed-form one-node RC divider, derived by scalar KCL."""
    return (1.0 / rs) / (1.0 / rs + 1.0 / rl
                         + 2j * np.pi * np.asarray(frequency) * (cs + cl))


def direct_mimo_convolution(x, taps):
    """Literal sum of 1-D linear convolutions, suitable for small oracles."""
    x = np.asarray(x)
    taps = np.asarray(taps)
    out = np.zeros((len(x) + len(taps) - 1, taps.shape[1]))
    for rx in range(taps.shape[1]):
        for tx in range(taps.shape[2]):
            out[:, rx] += np.convolve(x[:, tx], taps[:, rx, tx])
    return out


def scalar_nrz(bits, points_per_ui, rise_samples, delays, length):
    """Literal piecewise-linear source evaluation, one sample at a time."""
    bits = np.asarray(bits)
    output = np.zeros((length, bits.shape[1]))
    for tx in range(bits.shape[1]):
        for sample in range(length):
            position = sample - int(delays[tx])
            if position < 0 or position >= len(bits) * points_per_ui:
                continue
            bit_number, within_bit = divmod(position, points_per_ui)
            new = float(bits[bit_number, tx])
            old = float(bits[bit_number-1, tx]) if bit_number else 0.0
            if within_bit < rise_samples:
                output[sample, tx] = old + (new-old) * within_bit / rise_samples
            else:
                output[sample, tx] = new
    return output


def scalar_eye_counts(y, points_per_ui, voltage_edges):
    """Count each sample independently into voltage/time bins."""
    counts = np.zeros((len(voltage_edges) - 1, points_per_ui), np.int64)
    for k, voltage in enumerate(y):
        row = np.searchsorted(voltage_edges, voltage, side='right') - 1
        if voltage == voltage_edges[-1]:
            row -= 1
        if not 0 <= row < counts.shape[0]:
            raise ValueError('Sample outside independently supplied eye range')
        counts[row, k % points_per_ui] += 1
    return counts
