# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:43:27 2020

@author: DM1CR
"""

# Simulation of DDS spurs
# use IPython console and "%matplotlib auto"
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
plt.rcParams['figure.figsize'] = [18, 10]

def fftwin(s):
    return (np.kaiser(len(s), beta=38))
    #return (np.blackman(len(s)))

def dds(clockfreq,simulationsteps,phaseinc,phaseregbitcount,dacbitcount):
    phaseregmax=2**phaseregbitcount
    dacmax=2**(dacbitcount-1)
    timeincrement=1/clockfreq
    steps=np.arange(0,simulationsteps)
    time=steps*timeincrement                        # time vector
    phase=np.mod(steps*phaseinc,phaseregmax)        # phase vector
    scaledphase=phase*2*np.pi/phaseregmax           # phase < 2Pi
    fout=(2**(dacbitcount-1)*np.sin(scaledphase)).astype(int)   # time series output
    return  time, fout

def ddsspectrum(time, fout):
    spectrum=np.fft.rfft(fout*fftwin(fout))
    mag = np.abs(spectrum)
    YdB = 20 * np.log10(mag + 1e-200)
    freqs = np.fft.rfftfreq(len(fout), d=time[1]-time[0])
    return freqs, YdB

# Global parameters
phaseregbitcount=16
clockfreq=1e6
dacbitcount=20
simulationsteps=20000
phaseincrement = 5169
nyquistbandsforspectrumcalculation=6

fig, (ax1, ax2) = plt.subplots(2, 1)
fig.subplots_adjust(hspace=0.5, left=0.1, bottom=0.25)
t,fout = dds(clockfreq, simulationsteps, phaseincrement, \
             phaseregbitcount, dacbitcount)
l1, = ax1.plot(t, fout)
f, YdB = ddsspectrum(t, fout)
l2, = ax2.plot(f, YdB)

axcolor = 'lightgoldenrodyellow'
axphinc = plt.axes([0.1, 0.1, 0.75, 0.03], facecolor=axcolor)
axdacb = plt.axes([0.1, 0.15, 0.75, 0.03], facecolor=axcolor)
sphinc = Slider(axphinc, 'phase increment',0,0.5*2**phaseregbitcount,
              valinit=1024,valfmt='%1.4f', valstep=1)
sdacb = Slider(axdacb, 'DAC bits',1,26,
              valinit=10,valfmt='%2.f', valstep=1)

def update(val):
    t,fout = dds(clockfreq, simulationsteps, sphinc.val, \
                phaseregbitcount, sdacb.val)
    l1.set_ydata(fout)              # re-plot dds time data
    f, YdB = ddsspectrum(t, fout)
    l2.set_ydata(YdB)               # re-plot dds spectrum
    fig.canvas.draw_idle()

sphinc.on_changed(update)
sdacb.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    sphinc.reset()
button.on_clicked(reset)

plt.show()
