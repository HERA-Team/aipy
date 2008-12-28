#! /usr/bin/env python
import aipy._healpix, numpy

hp = aipy._healpix.HealpixBase()
print hp.order(), hp.nside(), hp.npix(), hp.scheme()
hp = aipy._healpix.HealpixBase(2**6)
print hp.order(), hp.nside(), hp.npix(), hp.scheme()
hp = aipy._healpix.HealpixBase(2**6, 'NEST')
print hp.order(), hp.nside(), hp.npix(), hp.scheme()
assert(hp.npix2nside(hp.npix()) == hp.nside())
px = numpy.arange(hp.npix())
px = hp.nest_ring_conv(px, 'RING')
hp.set_nside_scheme(hp.nside(), 'RING')
px = numpy.arange(10)
x,y,z = hp.px2crd(px, crdtype='vec')
th,phi = hp.px2crd(px, crdtype='pnt')
assert(numpy.all(px == hp.crd2px(x,y,z, interpolate=0)))
px2,wgt = hp.crd2px(x,y,z, interpolate=1)
z[-1] = numpy.Inf
try: hp.crd2px(x,y,z[:-1])
except(RuntimeError): pass
assert(numpy.all(px == hp.crd2px(th, phi, interpolate=0)))
px2,wgt = hp.crd2px(th, phi, interpolate=1)
