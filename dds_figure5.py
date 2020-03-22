# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 14:15:39 2020

@author: DM1CR
"""

# python port of illustration skripts for
# Circuit cellar #217 - The darker side : DDS generator
# Robert Lacoste

import numpy as np
import matplotlib.pyplot as plt

# Global parameters
phaseregbitcount=16
clockfreq=1e6
dacbitcount=8
simulationsteps=1000
nyquistbandsforspectrumcalculation=6

#-----------------------------
# First simulation

phaseincrement=5169

# Deducted parameters
phaseregmax=2**phaseregbitcount
dacmax=2**(dacbitcount-1)
timeincrement=1/clockfreq

# steps vector
steps=np.arange(0,simulationsteps)

# Time vector
time=steps*timeincrement

# Phase vector
phase=np.mod(steps*phaseincrement,phaseregmax)

# Output sine lookup vector
scaledphase=phase*2*np.pi/phaseregmax
fout=(dacmax*np.sin(scaledphase)).astype(int)

# sine vector  with better time resolution
fout_ext=np.arange(0,simulationsteps*nyquistbandsforspectrumcalculation)

for i in range(0,simulationsteps*nyquistbandsforspectrumcalculation-1):
    fout_ext[i]=fout[int((i-1)/nyquistbandsforspectrumcalculation)]


# calculate spectrum
spectrum=np.fft.fft(fout_ext)
spectrum=2*spectrum/(simulationsteps*nyquistbandsforspectrumcalculation)
# However as the input signal is real then only the first half of the FFT
# is actually useful (the last half is the complex conjugate of the first half) :
usefulspectrum=spectrum[0:int(len(spectrum)/2)]
spectrumamplitudes=np.abs(usefulspectrum)
spectrumamplitudes[0]=dacmax #just to scale the graph

# sin(x)/x reference
x=np.arange(0, int(nyquistbandsforspectrumcalculation*simulationsteps/2))
x=x*np.pi/simulationsteps
sinxoverx=x
for i in range(1,len(x)-1):
    sinxoverx[i]=dacmax*abs( np.sin(x[i])/x[i])
    
# Plot signals
plt.subplot(3,2,1)
plt.plot(phase[0:49])
plt.title('phase (N=16, P=16, B=8, W=5169)')
plt.subplot(3,2,3)
plt.plot(fout_ext[0:50*nyquistbandsforspectrumcalculation-1])
plt.title('Fout')
plt.subplot(3,2,5)
plt.plot(np.arange(0,len(spectrumamplitudes)), spectrumamplitudes,'k')
plt.plot(np.arange(0,len(spectrumamplitudes)), sinxoverx,'g')
plt.title('spectrum')


#-----------------------------
# second simulation

phaseincrement=15673

# Deducted parameters
phaseregmax=2**phaseregbitcount
dacmax=2**(dacbitcount-1)
timeincrement=1/clockfreq

# steps vector
steps=np.arange(0,simulationsteps)

# Time vector
time=steps*timeincrement

# Phase vector
phase=np.mod(steps*phaseincrement,phaseregmax)

# Output sine lookup vector
scaledphase=phase*2*np.pi/phaseregmax
fout=(dacmax*np.sin(scaledphase)).astype(int)

# sine vector  with better time resolution
fout_ext=np.arange(0,simulationsteps*nyquistbandsforspectrumcalculation)
for i in range(0, simulationsteps*nyquistbandsforspectrumcalculation-1):
    fout_ext[i]=fout[int((i-1)/nyquistbandsforspectrumcalculation)]

# calculate spectrum
spectrum=np.fft.fft(fout_ext)
spectrum=2*spectrum/(simulationsteps*nyquistbandsforspectrumcalculation)
# However as the input signal is real then only the first half of the FFT
# is actually useful (the last half is the complex conjugate of the first half) :
usefulspectrum=spectrum[0:int(len(spectrum)/2)]
spectrumamplitudes=np.abs(usefulspectrum)
spectrumamplitudes[0]=dacmax #just to scale the graph

# sin(x)/x reference
x=np.arange(0,int(nyquistbandsforspectrumcalculation*simulationsteps/2))
x=x*np.pi/simulationsteps
sinxoverx=x
for i in range(1,len(x)-1):
    sinxoverx[i]=dacmax*abs(np.sin(x[i])/x[i])

# Plot signals
plt.subplot(3,2,2)
plt.plot(phase[0:49])
plt.title('phase (N=16, P=16, B=8, W=15673)')
plt.subplot(3,2,4)
plt.plot(fout_ext[0:50*nyquistbandsforspectrumcalculation-1])
plt.title('Fout')
plt.subplot(3,2,6)
plt.plot(np.arange(0,len(spectrumamplitudes)), \
         spectrumamplitudes,'k', \
         np.arange(0,len(spectrumamplitudes)), \
         sinxoverx,'g')
plt.title('spectrum')




