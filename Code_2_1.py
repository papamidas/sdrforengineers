# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 20:03:02 2020

@author: DM1CR
"""

# Code 2.1

import matplotlib.pyplot as plt
import numpy as np

samples = 5000
period = 1000
step = 2
harmonics = 9
t = np.arange(samples)


for i in range(1,harmonics,step):
    wave = 1/i * np.sin(2 * np.pi * i * t/period)
    if i == 1:
        wavesum = wave
    else:
        wavesum += wave

print("Summiert bis zur " + str(harmonics) + ". Harmonischen")
fig, ax = plt.subplots()
ax.plot(t, wavesum)
ax.set(xlabel='time (sample)', ylabel='amplitude (a.u.)',
       title='Code 2.1 auf Seite 20')
ax.grid()
plt.show()

