#! /usr/bin/env python
from __future__ import absolute_import, print_function, division
import numpy as np
from matplotlib import pylab as p
import aipy as a

img = np.array([0,0,0,4,6,4,0,0,-2,-3,-2,0], dtype=np.float)
ker = np.array([3,2,0,0,0,0,0,0,0,0,0,2], dtype=np.float)

print('REAL TEST:')
print('img', img)
print('ker', ker)
cln, info = a.deconv.clean(img, ker)
print('cln', cln)
print(info)
print('-----------------------------------------------------------------')

img = np.array([0,0,0,4j,6j,4j,0,0,2+2j,3+3j,2+2j,0], dtype=np.complex)
ker = np.array([3,2,0,0,0,0,0,0,0,0,0,2], dtype=np.complex)

print('CMPLX TEST:')
print('img', img)
print('ker', ker)
mdl, info = a.deconv.clean(img, ker)
print('cln', cln)
print(info)
print('-----------------------------------------------------------------')

SIZE = 16
img = np.zeros((SIZE,), np.complex)
img[5] = 1+2j
img[10] = 3+4j
ker = np.zeros((SIZE,), np.complex)
ker[0] = 2j
ker[1] = 1+1j; ker[-1] = -1-1j
ker[2] = .5+.5j; ker[-2] = -.5-.5j
d1d = np.fft.ifft(np.fft.fft(img) * np.fft.fft(ker))
print('TEST:')
print('img', d1d)
print('ker', ker)
cln, info = a.deconv.clean(d1d, ker, maxiter=200)
print('cln', cln)
print(info)

p.plot(np.abs(cln), 'b')
p.plot(np.abs(img), 'k.')
p.plot(np.abs(img - cln), 'r')
p.show()
