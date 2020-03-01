# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 08:43:04 2020

@author: DM1CR
"""
import matplotlib.pyplot as plt
import numpy as np

def sample_plot(x1, x2, titlestr):
    plt.xlabel('Sample number')
    plt.ylabel('Amplitude')
    plt.xlim(x1, x2)
    plt.grid(b=True)
    plt.title(titlestr)


# Receive two vectors and return a vector resultant of
# convolution operation
def simple_conv(f,g):
    # Transform the vectors f and g in new vectors with the same length
    F = np.zeros( len(f) + len(g))
    G = np.zeros( len(f) + len(g))
    F[0:len(f)] = f
    G[0:len(g)] = g
    # Create a new vector C
    C = np.zeros(len(g)+len(f)-1)
    # FOR Loop to put the result of convolution between F and G vectors
    # in a new vector C. According to the convolution operation
    # characteristics, the length of a resultant vector of convolution
    # operation between two vector is the sum of vectors length minus 1
    for i in range(0, len(g)+len(f)-1): # range does not include last element
        # FOR Loop to walk through the vector F ang G
        for j in range(0, len(f)):
            if (i-j+1>0):
                C[i] = C[i] + F[j] * G[i-j+1]
    return C

samples1_len = 80
samples1_cnt = np.arange(0, samples1_len) #excluding stop
samples1_sin = 1.3*np.sin(samples1_cnt/samples1_len*2*4*np.pi)
samples1 = np.zeros(samples1_len)
samples1[8:69] = -samples1_sin[0:61] + 3.0/61*samples1_cnt[0:61]
m1 = samples1[68]/(68-72)
b1 = -m1 * 72
samples1[68:72] = m1 * samples1_cnt[68:72] + b1

# low-pass filter:

impulse1_len = 30
impulse1_phase = np.arange(0, np.pi*(impulse1_len+1)/impulse1_len,
                           np.pi/impulse1_len)
impulse1 = 0.055*np.sin(impulse1_phase)


plt.figure(num=1, figsize=(6, 4))
plt.plot(samples1,'o')
sample_plot(0, samples1_len, '')
plt.savefig('C2_16_fig1.png',bbox_inches='tight', pad_inches=0.5)

plt.figure(num=2, figsize=(2, 4))
plt.plot(impulse1,'o')
sample_plot(0, impulse1_len, 'Low-pass filter')
plt.savefig('C2_16_fig2.png',bbox_inches='tight', pad_inches=0.5)

filtered1 = simple_conv(samples1, impulse1)

plt.figure(num=3, figsize=(8, 4))
plt.plot(filtered1,'o')
plt.xticks(np.arange(0, 120, step=10))
plt.grid(b=True, which="minor")
plt.ylim(-2, 4)
sample_plot(0, len(filtered1), '')
plt.savefig('C2_16_fig3.png',bbox_inches='tight', pad_inches=0.5)

# high-pass filter:

impulse2 = -1 * impulse1
midindx = int(len(impulse2)/2)
impulse2[midindx] = 0
impulse2[midindx] = -np.sum(impulse2)


plt.figure(num=4, figsize=(6, 4))
plt.plot(samples1,'o')
sample_plot(0, samples1_len, '')
plt.savefig('C2_16_fig4.png',bbox_inches='tight', pad_inches=0.5)

plt.figure(num=5, figsize=(2, 4))
plt.plot(impulse2,'o')
sample_plot(0, impulse1_len, 'High-pass filter')
plt.savefig('C2_16_fig5.png',bbox_inches='tight', pad_inches=0.5)

filtered2 = simple_conv(samples1, impulse2)

plt.figure(num=6, figsize=(8, 4))
plt.plot(filtered2,'o')
plt.xticks(np.arange(0, 120, step=10))
plt.ylim(-2, 4)
sample_plot(0, len(filtered2), '')
plt.savefig('C2_16_fig6.png',bbox_inches='tight', pad_inches=0.5)

