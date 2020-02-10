# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 12:56:34 2020

@author: DM1CR
"""

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [5, 5]
import numpy as np
from sfdr import sfdr
   

# Specify plot parameters
#txtsize=10
#ltxtsize=9
#pwidth=4
#pheight=4
#pxoffset=0.65
#pyoffset=0.5
#markersize=5

deltat = 1e-8
fs = 1/deltat
t = np.arange(0, 1e-5-deltat, deltat)
fundamental=3959297
x = 10e-3*np.sin(2*np.pi*fundamental*t)

f1=plt.figure(1)
sfdr(x,fs)
plt.ylim((-400, 10))
plt.savefig('fig1.png',bbox_inches='tight', pad_inches=0.5)

f2=plt.figure(2)
plt.plot(t*10e6,x)
plt.xlabel('time (us)')
plt.ylabel('Amplitude')
plt.ylim((-1.1, 1.1))
plt.savefig('fig2.png',bbox_inches='tight', pad_inches=0.5)


bits=2**11
x = np.round(10e-3*bits*np.sin(2*np.pi*fundamental*t))/bits

f3=plt.figure(3)
sfdr(x,fs)
plt.ylim((-400, 10))
plt.savefig('fig3.png',bbox_inches='tight', pad_inches=0.5)

f4=plt.figure(4)
plt.plot(t*10e6,x)
plt.xlabel('time (us)')
plt.ylabel('Amplitude')
plt.ylim((-1.1, 1.1))
plt.savefig('fig4.png',bbox_inches='tight', pad_inches=0.5)


bits=2**11
x = np.round(bits*np.sin(2*np.pi*fundamental*t))/bits

f5=plt.figure(5)
sfdr(x,fs)
plt.ylim((-400, 10))
plt.savefig('fig5.png',bbox_inches='tight', pad_inches=0.5)

f6=plt.figure(6)
plt.plot(t*10e6,x)
plt.xlabel('time (us)')
plt.ylabel('Amplitude')
plt.ylim(-1.1, 1.1)
plt.savefig('fig6.png',bbox_inches='tight', pad_inches=0.5)


t = np.arange(0, 1e-4-deltat, deltat)
x = np.round(bits*np.sin(2*np.pi*fundamental*t))/bits

f7=plt.figure(7)
sfdr(x,fs)
plt.ylim((-400, 10))
plt.savefig('fig7.png',bbox_inches='tight', pad_inches=0.5)

f8=plt.figure(8, figsize = [25.6, 4.8])
plt.plot(t*10e6,x)
plt.xlabel('time (us)')
plt.ylabel('Amplitude')
plt.ylim((-1.1, 1.1))
plt.xlim((-10, 1010))
plt.savefig('fig8.png',bbox_inches='tight', pad_inches=0.5)

fundamental=4000000
x = np.round(bits*np.sin(2*np.pi*fundamental*t))/bits

f9=plt.figure(9)
sfdr(x,fs)
plt.ylim((-400, 10))
plt.savefig('fig9.png',bbox_inches='tight', pad_inches=0.5)

f10=plt.figure(10, figsize = [25.6, 4.8])
plt.plot(t*10e6,x)
plt.xlabel('time (us)')
plt.ylabel('Amplitude')
plt.ylim((-1.1, 1.1))
plt.xlim((-10, 1010))
plt.savefig('fig10.png',bbox_inches='tight', pad_inches=0.5)

ran = np.random.rand(len(t)) - 0.5
x = np.round(bits*np.sin(2*np.pi*fundamental*t) + ran)/bits

f11=plt.figure(11)
sfdr(x,fs)
plt.ylim((-400, 10))
plt.savefig('fig11.png',bbox_inches='tight', pad_inches=0.5)

f12=plt.figure(12, figsize = [25.6, 4.8])
plt.plot(t*10e6,x)
plt.xlabel('time (us)')
plt.ylabel('Amplitude')
plt.ylim((-1.1, 1.1))
plt.xlim((-10, 1010))
plt.savefig('fig12.png',bbox_inches='tight', pad_inches=0.5)

fundamental=3959297
x = np.round(bits*np.sin(2*np.pi*fundamental*t) + ran)/bits

f13=plt.figure(13)
sfdr(x,fs)
plt.ylim((-400, 10))
plt.savefig('fig13.png',bbox_inches='tight', pad_inches=0.5)

f14=plt.figure(14, figsize = [25.6, 4.8])
plt.plot(t*10e6,x)
plt.xlabel('time (us)')
plt.ylabel('Amplitude')
plt.ylim((-1.1, 1.1))
plt.xlim((-10, 1010))
plt.savefig('fig14.png',bbox_inches='tight', pad_inches=0.5)
