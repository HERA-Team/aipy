#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p

img = n.array([0,0,0,4,6,4,0,0,-2,-3,-2,0], dtype=n.float)
ker = n.array([3,2,0,0,0,0,0,0,0,0,0,2], dtype=n.float)

print 'REAL TEST:'
print 'img', img
print 'ker', ker
cln, info = a.deconv.clean(img, ker)
print 'cln', cln
print info
print '-----------------------------------------------------------------'

img = n.array([0,0,0,4j,6j,4j,0,0,2+2j,3+3j,2+2j,0], dtype=n.complex)
ker = n.array([3,2,0,0,0,0,0,0,0,0,0,2], dtype=n.complex)

print 'CMPLX TEST:'
print 'img', img
print 'ker', ker
mdl, info = a.deconv.clean(img, ker)
print 'cln', cln
print info
print '-----------------------------------------------------------------'

SIZE = 16
img = n.zeros((SIZE,), n.complex)
img[5] = 1+2j
img[10] = 3+4j
ker = n.zeros((SIZE,), n.complex)
ker[0] = 2j
ker[1] = 1+1j; ker[-1] = -1-1j
ker[2] = .5+.5j; ker[-2] = -.5-.5j
d1d = n.fft.ifft(n.fft.fft(img) * n.fft.fft(ker))
print 'TEST:'
print 'img', d1d
print 'ker', ker
cln, info = a.deconv.clean(d1d, ker, maxiter=200)
print 'cln', cln
print info

p.plot(n.abs(cln), 'b')
p.plot(n.abs(img), 'k.')
p.plot(n.abs(img - cln), 'r')
p.show()
