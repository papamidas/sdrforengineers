# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 17:14:53 2020

@author: DM1CR
"""

# additional calculations and plots for the article
# "Understanding the 'Phasing Method' of Single Sideband Demodulation"
# by Richard Lyons, 2012
# https://www.dsprelated.com/showarticle/176.php

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
plt.rcParams['figure.figsize'] = [18, 10]

def fftwin(s):
#    return (np.kaiser(len(s), beta=38))
    return (np.blackman(len(s)))

def calcdBFFT(s):
    s_fft = np.fft.fft(s * fftwin(s))
    s_dB = 20 * np.log10(np.abs(s_fft)+1e-200)
    return s_dB

def plotfft(ax, s, xlim1, xlim2,titlestr, xticklist, mutelist=[]):
    s_dB = calcdBFFT(s)
    freqs = np.fft.fftfreq(len(s), d=deltat)
    l, = ax.plot(freqs/1000, s_dB)
    ax.set_xlim(xlim1/1000, xlim2/1000)
    ax.set_ylim(-100,100)
    ax.set_xlabel('kHz')
    ax.xaxis.set_label_coords(1.02, -0.06)
    ax.set_xticks(xticklist)
    ax.set_title(titlestr)
    ax.grid(True)
    return l

usb = [[300,0.1],[1650,0.8],[3000,1]]   # USB frequencies and amplitudes
lsb = [[400,0.2],[1500,1.6],[2500,1.8]] # LSB frequencies and amplitudes
fc = 80000                              # Tx frequency, carrier frequency
c = 1                                   # carrier amplitude
fs = 2e6                                # sampling frequency
maxtim = 0.1                            # duration of the sampling block
deltat = 1/fs
t = np.arange(0, maxtim-deltat, deltat) # array of sampling times

ssb = np.zeros(len(t))                  # generating USB and LSB signals
for [freq, amp] in usb:
    ssb += amp*np.cos(2*np.pi*(freq+fc)*t)
for [freq, amp] in lsb:
    ssb += amp*np.cos(2*np.pi*(freq-fc)*t)

I = np.zeros(len(t))                    # calculating In-Phase Output
for [freq, amp] in usb:
    I += 0.5*amp*np.cos(2*np.pi*(freq+fc+fc)*t) + \
    0.5*amp*np.cos(2*np.pi*(freq+fc-fc)*t)
for [freq, amp] in lsb:
    I += 0.5*amp*np.cos(2*np.pi*(freq+fc+fc)*t) + \
    0.5*amp*np.cos(2*np.pi*(freq+fc-fc)*t)

fig, ax = plt.subplots()
plotfft(ax,ssb,-1.5*fc,1.5*fc,"receiver input signal",[-80,0,80])

fig, (ax1, ax2) = plt.subplots(2, 1)
fig.subplots_adjust(hspace=0.5, left=0.1, bottom=0.25)

axcolor = 'lightgoldenrodyellow'
axphs = plt.axes([0.1, 0.1, 0.75, 0.03], facecolor=axcolor)
axamp = plt.axes([0.1, 0.15, 0.75, 0.03], facecolor=axcolor)

sphs = Slider(axphs, 'phase err',-1.0,1.0,
              valinit=0,valfmt='%1.4f', valstep=0.001)
samp = Slider(axamp, 'amplitude err',0.99,1.01,
              valinit=1.0,valfmt='%1.4f',valstep=0.0001)

def calcQbar(ampQ, deltaPhiQ):          # calculating Quadrature Output
    Qbar = np.zeros(len(t))
    for [freq, amp] in usb:
        Qbar += ampQ*(-0.5*amp*np.cos(2*np.pi*(freq+fc+fc)*t+deltaPhiQ) - \
                0.5*amp*np.cos(2*np.pi*(freq+fc-fc)*t-deltaPhiQ))
    for [freq, amp] in lsb:
        Qbar += ampQ*(-0.5*amp*np.cos(2*np.pi*(freq+fc+fc)*t+deltaPhiQ) + \
                0.5*amp*np.cos(2*np.pi*(freq+fc-fc)*t-deltaPhiQ))
    return Qbar

ampQ = 1.0                              
deltaPhiQ = 0
Qbar = calcQbar(ampQ, deltaPhiQ)

l1=plotfft(ax1,I+Qbar,-5000,5000,
           "USB:I+Q(-90°) demodulated baseband signal",[-4,-2,0,2,4])  # USB
l2=plotfft(ax2,I-Qbar,-5000,5000,
           "LSB: I-Q(-90°) demodulated baseband signal",[-4,-2,0,2,4]) # LSB

def update(val):
    ampQ = samp.val
    deltaPhiQ = sphs.val / 180.0 * np.pi
    Qbar = calcQbar(ampQ, deltaPhiQ)
    l1.set_ydata(calcdBFFT(I+Qbar)) # re-claculate USB spectrum
    l2.set_ydata(calcdBFFT(I-Qbar)) # re-calculate LSB spectrum
    fig.canvas.draw_idle()

sphs.on_changed(update)
samp.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sphs.reset()
    samp.reset()
button.on_clicked(reset)

plt.show()
