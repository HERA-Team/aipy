#! /usr/bin/env python

from __future__ import print_function, division, absolute_import

import aipy as a, numpy as n
from matplotlib import pylab as p

SIZE = 100
aim = n.zeros((SIZE,SIZE), dtype=float)
aim[10,10] = 10.
aim[20:25,20:25] = 1.
aim[30:40,30:40] = .1
dbm = a.img.gaussian_beam(2, shape=aim.shape).astype(float)
dbm[0,0] = 2
dbm /= dbm.sum()

dim = n.fft.ifft2(n.fft.fft2(aim) * n.fft.fft2(dbm)).astype(float)

print('REAL TEST:')
cim, info = a.deconv.clean(dim, dbm)
print(info)
print('-----------------------------------------------------------------')

aim = n.zeros((SIZE,SIZE), dtype=complex)
aim[10,10] = 10.
aim[20:25,20:25] = 1j
aim[30:40,30:40] = .1+.1j
dbm = a.img.gaussian_beam(2, shape=aim.shape).astype(complex)
dbm[0,0] = 2
dbm /= dbm.sum()

dim = n.fft.ifft2(n.fft.fft2(aim) * n.fft.fft2(dbm))

print('COMPLEX TEST:')
cim, info = a.deconv.clean(dim, dbm)
print(info)
print('-----------------------------------------------------------------')

p.subplot(221)
dat = n.log(n.abs(dim))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.subplot(222)
dat = n.log(n.abs(a.img.recenter(dbm, (SIZE/2,SIZE/2))))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.subplot(223)
dat = n.log(n.abs(cim))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.subplot(224)
dat = n.log(n.abs(aim-cim))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.show()
