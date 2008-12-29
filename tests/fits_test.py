#! /usr/bin/env python
import numpy as n, aipy as a, pyfits

data = n.arange(200**2, dtype=n.float); data.shape = (200,200)
a.img.to_fits('test_img.fits', data, clobber=True)
d,h = a.img.from_fits('test_img.fits')
d = d.squeeze()
print d.shape
print d
print h
