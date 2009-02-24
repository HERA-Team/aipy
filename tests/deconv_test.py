#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p

SIZE = 100
NOISE = .001
i = n.zeros((SIZE,SIZE), n.float)
i[10,10] = 10.
i[20:25,20:25] = 1.
i[30:40,30:40] = .1
b = a.img.gaussian_beam(2, shape=i.shape)
b[0,0] = .05

d = n.abs(n.fft.ifft2(n.fft.fft2(i) * n.fft.fft2(b)))
ns = n.random.normal(scale=NOISE, size=i.shape)
d = n.abs(d + ns)

print 'Clean'
p.subplot(221)
c,info = a.deconv.clean(d, b, verbose=True)
p.title('CLEAN')
p.imshow(n.log10(c), vmin=-5, vmax=1)

print 'LSQ'
p.subplot(222)
c,info = a.deconv.lsq(d, b, verbose=True)
p.title('LSQ')
p.imshow(n.log10(c), vmin=-5, vmax=1)

print 'MEM'
p.subplot(223)
c,info = a.deconv.maxent(d, b, n.var(d**2)*.5, verbose=True)
p.title('MEM')
p.imshow(n.log10(c), vmin=-5, vmax=1)

print 'Anneal'
p.subplot(224)
c,info = a.deconv.anneal(d, b, verbose=True)
p.title('Anneal')
p.imshow(n.log10(c), vmin=-5, vmax=1)
p.colorbar()

p.show()

