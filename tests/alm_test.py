#! /usr/bin/env python
import aipy.healpix, sys, numpy

alm = aipy.healpix.Alm(1,1)
alm = aipy.healpix.Alm(2,1)
try:
    print "Checking that exception gets thrown:"
    alm = aipy.healpix.Alm(1,2)
    print "Should have failed (but didn't)"
    sys.exit(0)
except(RuntimeError): print "Ok."
alm = aipy.healpix.Alm(2,2)
print 'Lmax:', alm.lmax(), 'Mmax:', alm.mmax()
data = alm.get_data()
data[2] = 100
try: alm.set_data(data[:-1])
except(ValueError): pass
alm.set_data(data)
print alm.get_data(data)
for l,m in zip(*alm.lm_indices()):
    print l,m, '->', alm[l,m], '->',
    alm[l,m] = l+m*1j
    print alm[l,m]
print alm.get_data()
alm = aipy.healpix.Alm(2,2)
alm[1,0] = 1
alm[1,0] = 1.
alm[1,0] = 1j
try: alm[1,0] = "hi" ; print alm[1,0]
except(ValueError): pass
alm = aipy.healpix.Alm(2,2)
alm[0,0] = 4
hmap = aipy.healpix.HealpixMap(64, scheme="RING")
hmap.from_alm(alm)
#hmap.SetData(alm.to_map(hmap.Nside(), hmap.Scheme()))
print "Actual coeffiecients:"
print alm.get_data()
alm = hmap.to_alm(2, 2) ; print alm.get_data()
#alm.from_map(hmap.map.astype(numpy.int), 1); print alm.get_data()
print "long test:"
alm.from_map(hmap.map.astype(numpy.long), 1); print alm.get_data()
print "float test:"
alm.from_map(hmap.map.astype(numpy.float32), 1); print alm.get_data()
print "double test:"
for i in [1,2,3]:
    print "iter %d:" % i
    alm.from_map(hmap.map.astype(numpy.double), i); print alm.get_data()
#hmap.to_fits('alm_test.fits')
