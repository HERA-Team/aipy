#! /usr/bin/env python

from __future__ import print_function, division, absolute_import

import aipy as a, numpy as n
from matplotlib import pylab as p

img = n.array([0,0,0,4,6,4,0,0,-2,-3,-2,0], dtype=float)
ker = n.array([3,2,0,0,0,0,0,0,0,0,0,2], dtype=float)

print('REAL TEST:')
print('img', img)
print('ker', ker)
cln, info = a.deconv.clean(img, ker)
print('cln', cln)
print(info)
print('-----------------------------------------------------------------')

img = n.array([0,0,0,4j,6j,4j,0,0,2+2j,3+3j,2+2j,0], dtype=complex)
ker = n.array([3,2,0,0,0,0,0,0,0,0,0,2], dtype=complex)

print('CMPLX TEST:')
print('img', img)
print('ker', ker)
mdl, info = a.deconv.clean(img, ker)
print('cln', cln)
print(info)
print('-----------------------------------------------------------------')

SIZE = 16
img = n.zeros((SIZE,), complex)
img[5] = 1+2j
img[10] = 3+4j
ker = n.zeros((SIZE,), complex)
ker[0] = 2j
ker[1] = 1+1j; ker[-1] = -1-1j
ker[2] = .5+.5j; ker[-2] = -.5-.5j
d1d = n.fft.ifft(n.fft.fft(img) * n.fft.fft(ker))
print('TEST:')
print('img', d1d)
print('ker', ker)
cln, info = a.deconv.clean(d1d, ker, maxiter=200)
print('cln', cln)
print(info)

p.plot(n.abs(cln), 'b')
p.plot(n.abs(img), 'k.')
p.plot(n.abs(img - cln), 'r')
p.show()
