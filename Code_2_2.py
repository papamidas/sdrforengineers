# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 21:35:34 2020

@author: DM1CR
"""

# Code 2.2

import math
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 5]
import numpy as np

Fs = 1000  # Sample rate (Hz)
Fa = 905  # Input Frequency (Hz)
# Determine Nyquist zones
zone = 1 + math.floor(Fa / (Fs/2))
alias = Fa % Fs
if (zone % 2) == 0: # 2nd, 4th, 6th, ... Nyquist Zone
                    # Its not really a negative amplitude,
                    # but it is 180 degrees out of phase,
                    # which makes it harder to see on the time domain side,
                    # so we cheat to make the graphs look better.
    alias = -(Fs - alias)/Fs
else:               # 3rd, 5th, 7th, ... Nyquist Zone
    alias = (alias)/Fs

# Create the analog/time domain and digital sampling vectors
N = 2*1/abs(alias) + 1    # Number of Digital samples
points = 256              # Analog points between digital samples
analogIndexes = np.arange(0, N-1, 1/points)
samplingIndexes = np.arange(1, analogIndexes.size, points)
wave = np.sin(2*np.pi*Fa/Fs*analogIndexes)

fig, ax = plt.subplots()
ax.set(xlabel='Digital Samples', ylabel='Amplitude')
# Plot analog  input signal and sampled points
ax.bar(analogIndexes[samplingIndexes], wave[samplingIndexes], .1,
             label = 'Digital Sampling')
ax.plot(analogIndexes[samplingIndexes], wave[samplingIndexes], 'o',
              label = 'Digital Samples')
ax.plot(analogIndexes, wave, label = 'Analog Input')
# Plot digitally recreated signal
ax.plot(analogIndexes, np.sin(2*np.pi*alias*analogIndexes), '--',
              label = 'Digital Reconstruction')
ax.legend()
plt.title("Actual Input Frequency = %i\n" \
      "Measured Frequency = %i" % (Fa, abs(alias * Fs)))
plt.show()

