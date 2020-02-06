# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 15:45:28 2020

@author: DM1CR
"""
# Code 2.3, 2.4, 2.5, 2.6, 2.7 (updownsample.m)

import matplotlib.pyplot as plt
plt.rc('figure', max_open_warning = 0)
import numpy as np
from scipy.signal import find_peaks_cwt, firls, lfilter, upfirdn

def length(a):
    return a.size

def time_plot(x1, x2, txtsize, ltxtsize, pwidth, pheight, pxoffset,
    pyoffset, markersize, titlestr):
    #persistent file
    plt.xlabel('Discrete Time (n)')
    plt.ylabel('Signal Amplitude')
    plt.ylim(-1.5, 1.5)
    plt.xlim(x1, x2)
    plt.title(titlestr)

def fft_plot(data, points, txtsize, ltxtsize, pwidth, pheight,
                      pxoffset, pyoffset, markersize, titlestr, units):
    #persistent file
    fft_data = do_fft(data)
    fft_axis = do_fs(data)
    #print("Länge data: %d" % length(data))
    #print("Länge fft_data: %d" % length(fft_data))
    #print("Länge fft_axis: %d" % length(fft_axis))
    ydata = fft_data[np.int32(length(data)/2):length(data)]
    xdata = fft_axis[np.int32(length(data)/2):length(data)]
    peakind = find_peaks_cwt(xdata, np.arange(1,10))
    
    psor=ydata[peakind]
    lsor=xdata[peakind]
    
    #'NPeaks', points, 'SortStr', 'descend')
    
    if np.mean(psor) > -23:
        (M,I) = min(lsor)
        plt.plot(fft_axis, fft_data, lsor(I), psor(I), 'o')
        plt.xlabel(r'$\frac{f_s}{2}$')
        tstr = '%2.1fdB @ %.3f' % (psor(I), lsor(I))
        plt.text(lsor[I]+.02, psor[I], tstr + r'$\frac{f_s}{2}$')
    else:
        plt.plot(fft_axis, fft_data)
    
    plt.xlabel( 'Frequency ' + units)
    plt.ylabel('Signal Amplitude (dB)')
    plt.ylim(-200, 0)
    plt.xlim(-0.5, 0.5)
    plt.xticks((-.5, -.25,  0,  .25,  .5))
    plt.grid(b=True)
    plt.title(titlestr)
    
def do_fft(indata):
    L = length (indata) # Window Length of FFT
    in_HannWnd = indata * np.hanning(L)
    out=20*(np.log10(np.abs(np.fft.fftshift(np.fft.fft(in_HannWnd,L))/(L/2))))
    return out

def do_fs(indata):
    # numpy.arange([start, ]stop, [step, ]dtype=None)
    # Stop-Wert ist nicht mehr Bestandteil des Arrays!
    out = (np.arange(0,length(indata),1)/(0.5*length(indata))-1)/2
    return out


# Specify plot parameters
txtsize=20
ltxtsize=9
pwidth=4 #4
pheight=4 #4
pxoffset=1 #0.65
pyoffset=1 #0.5
markersize=5 

# number filter taps
taps = 127
# where to plot from
start = np.int32((taps/100)+1)*100
inc = 50

Fs1 = 10*2*np.pi
# Create deterministic and stochastic digital data streams
n = np.arange(0, 100-(1/Fs1)+1, 1/Fs1)                # Time index vector
sin_wave = np.sin(5*n*2*np.pi)               # Generation of sinusoidal signal
random = 2*np.around(np.random.rand(n.size))-1  # Random string of +1/-1 values

plt.figure(1)
plt.plot(n[start:start+inc]*Fs1, sin_wave[start:start+inc], 'o',
    n[start:start+inc]*Fs1, sin_wave[start:start+inc])
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset,
          pyoffset, markersize,
          'Original Sinusoidal Signal : Time domain')

plt.figure(2)
plt.plot(n[start:start+inc]*Fs1, random[start:start+inc],'o',
    n[start:start+inc]*Fs1, random[start:start+inc])
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset,
          pyoffset, markersize,
          'Original Random Binary Signal : Time domain')

#print("start: %d" % start)
#print("Länge sin_wave: %d" % length(sin_wave))

plt.figure(3)
fft_plot(sin_wave[start:], 1, txtsize, ltxtsize, pwidth, pheight,
    pxoffset,  pyoffset, markersize,
    'Original Sinusoidal Signal : Fourier domain', r'$f_{s1}$')
         
plt.figure(4)
fft_plot(random[start:], 1, txtsize, ltxtsize, pwidth, pheight,
    pxoffset, pyoffset, markersize,
    'Original Random Binary Signal : Fourier domain', r'$f_{s1}$')


# Create lowpass filter and apply it to both data streams
# MATLAB:
#b = firls(n,f,a),
#     n is the FIR filter order
#     f is a vector of pairs of frequency points,
#     a is a vector containing the desired amplitude at the points in f
# SciPy:
#scipy.signal.firls(numtaps, bands, desired, weight=None, nyq=None, fs=None)
#   numtaps : int
#        The number of taps in the FIR filter. numtaps must be odd.
#   bands : array_like
#        A monotonic nondecreasing sequence containing the band edges in Hz.
#        All elements must be non-negative and less than or equal to the
#        Nyquist frequency given by nyq.
#   desired : array_like
#        A sequence the same size as bands containing the desired gain at
#        the start and end point of each band.
#   weight : array_like, optional
#        A relative weighting to give to each band region when solving
#       the least squares problem. weight has to be half the size of bands.
#   nyq : float, optional
#        Deprecated. Use `fs` instead. Nyquist frequency. Each frequency
#        in bands must be between 0 and nyq (inclusive). Default is 1.
#   fs : float, optional
#        The sampling frequency of the signal. Each frequency in bands
#        must be between 0 and fs/2 (inclusive). Default is 2.
coeffs1 = firls(taps,(0,0.2,0.22,1),(1,1,0,0)) # FIR filter coefficients
# MATLAB:
#y = filter(b,a,x) filters the input data x using a rational transfer
# function defined by the numerator and denominator coefficients b and a.
# SciPy:
#scipy.signal.lfilter(b, a, x, axis=-1, zi=None)
#    b : array_like
#        The numerator coefficient vector in a 1-D sequence.
#    a : array_like
#        The denominator coefficient vector in a 1-D sequence.
#        If a[0] is not 1, then both a and b are normalized by a[0].
#    x : array_like
#        An N-dimensional input array.
#    axis : int, optional
#        The axis of the input data array along which to apply the 
#        linear filter. The filter is applied to each subarray along
#        this axis. Default is -1.
#    zi : array_like, optional
#        Initial conditions for the filter delays. It is a vector
#        (or array of vectors for an N-dimensional input) of length
#        max(len(a), len(b)) - 1. If zi is None or is not given then
#        initial rest is assumed. See lfiltic for more information.
sin_bwlimited = lfilter(coeffs1,1,sin_wave)
random_bwlimited = lfilter(coeffs1,1,random)

plt.figure(5)
plt.plot(n[start:start+inc]*Fs1, sin_bwlimited[start:start+inc], '-o')
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset,
    pyoffset, markersize, 'Band Limited Sinusoidal Signal : Time domain')

plt.figure(6)
plt.plot(n[start:start+inc]*Fs1, random_bwlimited[start:start+inc],'o')
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset,
    pyoffset, markersize, 'Band Limited Random Binary Signal : Time domain')

plt.figure(7)
fft_plot(sin_bwlimited[start:], 1, txtsize, ltxtsize, pwidth, pheight,
    pxoffset, pyoffset, markersize,
    'Band Limited Sinusoidal Signal : Fourier domain', r'$f_{s1}$')

plt.figure(8)
fft_plot(random_bwlimited[start:], 1, txtsize, ltxtsize, pwidth,
    pheight, pxoffset, pyoffset, markersize,
    'Band Limited Random Binary Signal : Fourier domain', r'$f_{s1}$')

# MATLAB
# y = upsample(xn)
#     increases the sampling rate of x by inserting (n – 1) zeros
#     between samples.
# SciPy
#scipy.signal.upfirdn(h, x, up=1, down=1, axis=-1)
#    h : array_like
#        1-dimensional FIR (finite-impulse response) filter coefficients.
#    x : array_like
#        Input signal array.
#   up : int, optional
#        Upsampling rate. Default is 1.
# down : int, optional
#        Downsampling rate. Default is 1.
# axis : int, optional
#        The axis of the input data array along which to apply the
#        linear filter. The filter is applied to each subarray along
#        this axis. Default is -1.
N = 5
sin_up = upfirdn([1], sin_bwlimited, up=N)
random_up = upfirdn([1], random_bwlimited, up=N)
Fs2 = Fs1 * N

plt.figure(9)
plt.plot(n[start*N:(start+inc)*N]*Fs1, sin_up[start*N:(start+inc)*N], 'o')
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight,
    pxoffset, pyoffset, markersize,
    'Improperly Upsampled (%d) Sinusoidal Signal : Time domain' % N)

plt.figure(10)
plt.plot(n[start*N:(start+inc)*N]*Fs1, random_up[start*N:(start+inc)*N],'o')
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight,
    pxoffset, pyoffset, markersize,
    'Improperly Upsampled (%d) Random Binary Signal : Time domain' % N)

plt.figure(11)
fft_plot(sin_up[start:], N, txtsize, ltxtsize, pwidth, pheight,
    pxoffset, pyoffset, markersize,
    'Improperly Upsampled (%d) Sinusoidal Signal : Fourier domain' % N,
    r'$f_{s2}$')

plt.figure(12)
fft_plot(random_up[start:], N, txtsize, ltxtsize, pwidth, pheight,
    pxoffset, pyoffset, markersize,
    'Improperly Upsampled (%d) Random Binary Signal : Fourier domain' % N,
    r'$f_{s2}$')

# Attempt to downsampling by M without filtering
# This is incorrect, but is instructive to show what artifacts occur
M = 3
sin_up_down = upfirdn([1], sin_up, down=M)
random_up_down = upfirdn([1], random_up, down=M)
Fs3 = Fs2/M

plt.figure(13)
plt.plot(n[int(start*N/M):int((start+inc)*N/M)]*Fs1,
         sin_up_down[int(start*N/M):int((start+inc)*N/M)], 'o')
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight,
          pxoffset, pyoffset, markersize,
          'Improperly Upsampled (%d) then Downsampled (%d) '\
          'Sinusoidal Signal : Time domain' % (N, M))

plt.figure(14)
plt.plot(n[int(start*N/M):int((start+inc)*N/M)]*Fs1,
     random_up_down[int(start*N/M):int((start+inc)*N/M)],'o')
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight,
          pxoffset, pyoffset, markersize,
          'Improperly Upsampled (%d) then Downsampled (%d) '\
          'Random Binary Signal : Time domain' % (N, M))

plt.figure(15)
fft_plot(sin_up_down[start:], N, txtsize, ltxtsize, pwidth, pheight,
         pxoffset, pyoffset, markersize,
         'Improperly Upsampled (%d) then Downsampled (%d) '\
         'Sinusoidal Signal : Fourier domain' % (N, M), r'$f_{s3}$')

plt.figure(16)
fft_plot(random_up_down[start:], N, txtsize, ltxtsize, pwidth, pheight,
         pxoffset, pyoffset, markersize,
         'Improperly Upsampled (%d) then Downsampled (%d) '\
         'Random Binary Signal : Fourier domain' % (N, M), r'$f_{s3}$')

# Lowpass filtering of baseband periodic replica followed by downsampling
# (correct approach)
coeffs2 = firls(taps,(0,0.15,0.25,1),(N,N,0,0)) # FIR filter coefficients
sin_up_filtered = lfilter(coeffs2,1,sin_up)
sin_up_filtered_down = upfirdn([1], sin_up_filtered, down=M)
random_up_filtered = lfilter(coeffs2,1,random_up)
random_up_filtered_down = upfirdn([1], random_up_filtered, down=M)
start = start + np.int32(taps/(2*N))

plt.figure(17)
plt.plot(n[start*N:(start+inc)*N]*Fs1,
           sin_up_filtered[start*N:(start+inc)*N], 'o')
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight,
          pxoffset, pyoffset, markersize,
          'Properly Upsampled (%d) and Filtered '\
          'Sinusoidal Signal : Time domain' % N)

plt.figure(18)
plt.plot(n[start*N:(start+inc)*N]*Fs1,
           random_up_filtered[start*N:(start+inc)*N],'o')
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight,
          pxoffset, pyoffset, markersize,
          'Properly Upsampled (%d) and Filtered '\
          'Random Binary Signal : Time domain' % N)

plt.figure(19)
fft_plot(sin_up_filtered[start:], 1, txtsize, ltxtsize, pwidth, pheight,
         pxoffset, pyoffset, markersize,
         'Properly Upsampled (%d) and Filtered '\
         'Sinusoidal Signal : Fourier domain' % N, r'$f_{s2}$')

plt.figure(20)
fft_plot(random_up_filtered[start:], 1, txtsize, ltxtsize, pwidth, pheight,
         pxoffset, pyoffset, markersize,
         'Properly Upsampled (%d) and Filtered '\
         'Random Binary Signal : Fourier domain' % N, r'$f_{s2}$')

plt.figure(21)
plt.plot(n[int(start*N/M):int((start+inc)*N/M)]*Fs1,
         sin_up_filtered_down[int(start*N/M):int((start+inc)*N/M)], 'o')
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight,
          pxoffset, pyoffset, markersize,
          'Properly Upsampled (%d), Filtered, and Downsampled (%d) '\
          'Sinusoidal Signal : Time domain' % (N, M))

plt.figure(22)
plt.plot(n[int(start*N/M):int((start+inc)*N/M)]*Fs1,
         random_up_filtered_down[int(start*N/M):int((start+inc)*N/M)], 'o')
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight,
          pxoffset, pyoffset, markersize,
          'Properly Upsampled (%d), Filtered, and Downsampled (%d) '\
          'Random Binary Signal : Time domain' % (N, M))

plt.figure(23)
fft_plot(sin_up_filtered_down[start:], 1, txtsize, ltxtsize, pwidth, pheight,
         pxoffset, pyoffset, markersize,
         'Properly Upsampled (%d), Filtered and Downsampled (%d) '\
         'Sinusoidal Signal : Fourier domain' % (N, M), r'$f_{s3}$')

plt.figure(24)
fft_plot(random_up_filtered_down[start:], 1, txtsize, ltxtsize, pwidth,
         pheight, pxoffset, pyoffset, markersize,
         'Properly Upsampled (%d), Filtered and Downsampled (%d) '\
         'Random Binary Signal : Fourier domain' % (N,M), r'$f_{s3}$')
