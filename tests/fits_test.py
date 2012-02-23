#! /usr/bin/env python
import numpy as n, aipy as a, pyfits

im = a.img.Img(size=200, res=.5)
L,M = im.get_LM(center=(200,200))
d_ra = n.arcsin(L[199,199]) * a.img.rad2deg
d_dec = n.arcsin(M[199,199]) * a.img.rad2deg

a.img.to_fits('test_img.fits', L.filled(0), clobber=True, 
    ra=10, d_ra=d_ra, dec=40, d_dec=d_dec)
d,h = a.img.from_fits('test_img.fits')
d = d.squeeze()
print d.shape
print d
print h
