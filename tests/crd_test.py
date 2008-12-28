#! /usr/bin/env python

import aipy, ephem, numpy

ra,dec,epoch = ('19:59:28.3', '+40:44:02', ephem.J2000)
cyg = ephem.Equatorial(ra, dec, epoch=epoch)
print cyg.get()
ra,dec = cyg.get()
v_eq = aipy.coord.radec2eq((ra,dec))
print v_eq, 'J2000'
m = aipy.coord.convert_m('eq','eq', iepoch=epoch, oepoch='2007.5')
pv_eq = numpy.dot(m, v_eq)
print pv_eq, 'J2007.5'
ra,dec = aipy.coord.eq2radec(pv_eq)
pcyg = ephem.Equatorial(ra, dec, epoch='2007.5')
print pcyg.get()
pcyg = ephem.Equatorial(cyg, epoch='2007.5')
print 'Correct answer is: (19:59:43.88, 40:45:16.9)'
