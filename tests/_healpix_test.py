# -*- coding: utf-8 -*-
import unittest, aipy._healpix as h, numpy as n

class TestHealpix(unittest.TestCase):
    def setUp(self):
        self.hpb = h.HealpixBase()
    def test_order(self):
        """Test HEALpix order funtional attribute"""
        self.assertEqual(self.hpb.order(), -1)
    def test_nside(self):
        """Test HEALpix nside funtional attribute"""
        self.assertEqual(self.hpb.nside(), 0)
    def test_npix(self):
        """Test HEALpix npix funtional attribute"""
        self.assertEqual(self.hpb.npix(), 12*self.hpb.nside()**2)
    def test_scheme(self):
        """Test HEALpix scheme funtional attribute"""
        self.assertEqual(self.hpb.scheme(), 'RING')
    def test_npix2nside(self):
        """Test HEALpix npix2nside funtional attribute"""
        self.assertEqual(self.hpb.npix2nside(12*2**12), 2**6)

if False:
  class TestMemLeaks(unittest.TestCase):
    def setUp(self):
        self.hpb = h.HealpixBase(nside=256)
    def test_create(self):
        while True: hpb = h.HealpixBase(nside=256)
    def test_crd2pix(self):
        one = 0.99
        x = n.arange(-one, one, .01, dtype=n.double)
        y = n.zeros_like(x)
        z = n.sqrt(1 - x**2 - y**2)
        while True: crds = self.hpb.crd2px(x, y, z)

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy._healpix unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestHealpix))
        #self.addTests(loader.loadTestsFromTestCase(TestMemLeaks))

if __name__ == '__main__':
    unittest.main()

#hp = aipy._healpix.HealpixBase(2**6, 'NEST')
#print hp.order(), hp.nside(), hp.npix(), hp.scheme()
#px = numpy.arange(hp.npix())
#px = hp.nest_ring_conv(px, 'RING')
#hp.set_nside_scheme(hp.nside(), 'RING')
#px = numpy.arange(10)
#x,y,z = hp.px2crd(px, ncrd=3)
#th,phi = hp.px2crd(px, ncrd=2)
#assert(numpy.all(px == hp.crd2px(x,y,z, interpolate=0)))
#px2,wgt = hp.crd2px(x,y,z, interpolate=1)
#z[-1] = numpy.Inf
#try: hp.crd2px(x,y,z[:-1])
#except(RuntimeError): pass
#assert(numpy.all(px == hp.crd2px(th, phi, interpolate=0)))
#px2,wgt = hp.crd2px(th, phi, interpolate=1)
