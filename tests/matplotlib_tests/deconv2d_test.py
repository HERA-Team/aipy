#! /usr/bin/env python
from __future__ import absolute_import, print_function, division
import numpy as np
from matplotlib import pylab as p
import aipy as a

SIZE = 100
aim = np.zeros((SIZE,SIZE), dtype=np.float)
aim[10,10] = 10.
aim[20:25,20:25] = 1.
aim[30:40,30:40] = .1
dbm = a.img.gaussian_beam(2, shape=aim.shape).astype(np.float)
dbm[0,0] = 2
dbm /= dbm.sum()

dim = np.fft.ifft2(np.fft.fft2(aim) * np.fft.fft2(dbm)).astype(np.float)

print('REAL TEST:')
cim, info = a.deconv.clean(dim, dbm)
print(info)
print('-----------------------------------------------------------------')

aim = np.zeros((SIZE,SIZE), dtype=np.complex)
aim[10,10] = 10.
aim[20:25,20:25] = 1j
aim[30:40,30:40] = .1+.1j
dbm = a.img.gaussian_beam(2, shape=aim.shape).astype(np.complex)
dbm[0,0] = 2
dbm /= dbm.sum()

dim = np.fft.ifft2(np.fft.fft2(aim) * np.fft.fft2(dbm))

print('COMPLEX TEST:')
cim, info = a.deconv.clean(dim, dbm)
print(info)
print('-----------------------------------------------------------------')

p.subplot(221)
dat = np.log(np.abs(dim))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.subplot(222)
dat = np.log(np.abs(a.img.recenter(dbm, (SIZE/2,SIZE/2))))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.subplot(223)
dat = np.log(np.abs(cim))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.subplot(224)
dat = np.log(np.abs(aim-cim))
p.imshow(dat, vmin=dat.max()-6, vmax=dat.max())
p.colorbar(shrink=.5)
p.show()
