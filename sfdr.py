# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 12:56:34 2020

@author: DM1CR
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [5, 5]

def sfdr(x,fs): # attempt of a python version of MATLABs sfdr function 
    xw = x * np.kaiser(len(x),beta=38) /len(x)
    xw -= np.mean(xw)
    Y = np.fft.rfft(xw)
    freqs = np.fft.rfftfreq(len(xw), d=1.0/fs)
    mag = np.abs(Y)
    YdB = 20 * np.log10(mag)
    peakind = find_peaks(YdB, distance = 5)
    pksf=freqs[peakind[0]]
    pksY=YdB[peakind[0]]
    isorted = np.argsort(pksY)
    sfdrval = pksY[isorted[-1]] - pksY[isorted[-2]]
    #fig, ax = plt.subplots()
    ax = plt.gca()
    pkfa = pksf[isorted[-1]]
    pkYa = pksY[isorted[-1]]
    pkfb = pksf[isorted[-2]]
    pkYb = pksY[isorted[-2]]
    plt.fill_between((0,fs/2),(pkYb,pkYb),(pkYa,pkYa), label = 'SFDR',
                     color = "lightblue")    
    ax.plot(pkfa, pkYa, marker="s", label = 'fundamental')
    ax.plot(pkfb, pkYb, marker="s", label = 'spurs')
    ax.plot(freqs, YdB)
    ax.set(xlabel = 'Frequency (Hz)', ylabel = 'Power (dB)',
           title = "SFDR %.2f dB" % sfdrval)
    ax.set_xlim(0, fs / 2)
    ax.set_ylim(-400, 10)    
    ax.legend(loc = "upper right")
    return sfdrval

if __name__ == "__main__":
    
    deltat = 1e-8
    fs = 1/deltat
    t = np.arange(0, 1e-5-deltat, deltat)
    fundamental = 3959297  # Prime number
    x = 10e-3*np.sin(2*np.pi*fundamental*t)
    sfdr(x,fs)
    
