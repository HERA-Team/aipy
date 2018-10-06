# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import unittest, aipy.healpix as h, numpy as n

class TestHealpixBase(unittest.TestCase):
    def setUp(self):
        self.hpb = h.HealpixBase(nside=32, scheme='RING')
    def test_order(self):
        """Test HEALpix order functional attribute"""
        self.assertEqual(self.hpb.order(), 5)
    def test_nside(self):
        """Test HEALpix nside functional attribute"""
        self.assertEqual(self.hpb.nside(), 32)
    def test_npix(self):
        """Test HEALpix npix functional attribute"""
        self.assertEqual(self.hpb.npix(), 12*self.hpb.nside()**2)
    def test_scheme(self):
        """Test HEALpix scheme functional attribute"""
        self.assertEqual(self.hpb.scheme(), 'RING')
    def test_npix2nside(self):
        """Test HEALpix npix2nside functional attribute"""
        self.assertEqual(self.hpb.npix2nside(12*2**12), 2**6)
    def test_set_nside_scheme(self):
        self.hpb.set_nside_scheme(64,'NEST')
        self.assertEqual(self.hpb.nside(), 64)
        self.assertEqual(self.hpb.scheme(), 'NEST')
        self.hpb.set_nside_scheme(32,'RING')
        self.assertEqual(self.hpb.nside(), 32)
        self.assertEqual(self.hpb.scheme(), 'RING')
    def test_nest_ring_conv(self):
        ipx = n.array([0])
        px = self.hpb.nest_ring_conv(ipx, 'NEST')
        self.assertEqual(px, 1023)
        self.hpb.set_nside_scheme(32, 'NEST')
        ipx = n.array([0])
        px = self.hpb.nest_ring_conv(ipx, 'RING')
        self.assertEqual(px, 5968)
    def test_ang2px(self):
        th,ph = n.linspace(0,n.pi,10), n.linspace(-n.pi,n.pi,10)
        px = self.hpb.crd2px(th,ph)
        self.assertEqual(len(px), len(th))
        n.testing.assert_allclose(px, n.array([    2,   398,  1375,  3114,  5049,  7239,  9173, 10912, 11889, 12286]))
    def test_vec2px(self):
        x,y = n.linspace(-.5,.5,10), n.linspace(-.5,.5,10)
        z = 1 - n.sqrt(x**2 + y**2)
        px = self.hpb.crd2px(x,y,z)
        self.assertEqual(len(px), len(x))
        n.testing.assert_allclose(px, n.array([3728, 2192, 1069,  247,   19,   13,  225, 1023, 2128, 3664]))
    def test_px2vec(self):
        x,y = n.linspace(-.5,.5,10), n.linspace(-.5,.5,10)
        z = 1 - n.sqrt(x**2 + y**2)
        px = n.array([0])
        x, y, z = self.hpb.px2crd(px)
        self.assertEqual(len(px), len(x))
        n.testing.assert_allclose(x, n.array([0.018041]), rtol=1e-4)
        n.testing.assert_allclose(y, n.array([0.018041]), rtol=1e-4)
    def test_ang2px_interp(self):
        th,ph = n.linspace(0,n.pi,3), n.linspace(-n.pi,n.pi,3)
        px,wgt = self.hpb.crd2px(th,ph,interpolate=True)
        self.assertEqual(len(px), len(th))
        self.assertEqual(wgt.shape, (len(th),4))
        n.testing.assert_allclose(px, n.array(
            # The ordering of pixels can be healpix version dependent
              #[[    2,     3,     1,     0],
              # [ 6208,  6080,  6207,  5952],
              # [12286, 12285, 12287, 12284]]))
              [[    3,     0,     1,     2],
               [ 6207,  6080,  6208,  6209],
               [12285, 12286, 12287, 12284]]))
        n.testing.assert_allclose(wgt, n.array(
            # The ordering of pixels can be healpix version dependent
              #[[ 0.25,  0.25,  0.25,  0.25],
              # [ 0.25,  0.25,  0.25,  0.25],
              # [ 0.25,  0.25,  0.25,  0.25]]))
              [[ 0.25,  0.25,  0.25,  0.25],
               [ 0.5 ,  0.5 ,  0.  ,  0.  ],
               [ 0.25,  0.25,  0.25,  0.25]]))
    def test_vec2px_interp(self):
        x,y = n.linspace(-.5,.5,3), n.linspace(-.5,.5,3)
        z = 1 - n.sqrt(x**2 + y**2)
        px,wgt = self.hpb.crd2px(x,y,z,interpolate=True)
        self.assertEqual(len(px), len(x))
        self.assertEqual(wgt.shape, (len(x),4))
        n.testing.assert_allclose(px, n.array(
            # The ordering of pixels can be healpix version dependent
              #[[3984, 3856, 3855, 3728],
              # [   0,    1,    3,    2],
              # [3920, 3792, 3791, 3664]]))
              [[3728, 3729, 3855, 3856],
               [   1,    2,    3,    0],
               [3664, 3665, 3791, 3792]]))
        n.testing.assert_allclose(wgt, n.array(
            # The ordering of pixels can be healpix version dependent
              #[[ 0.099602,  0.215996,  0.215996,  0.468407],
              # [ 0.25    ,  0.25    ,  0.25    ,  0.25    ],
              # [ 0.099602,  0.215996,  0.215996,  0.468407]]), rtol=1e-5)
              [[ 0.367711,  0.      ,  0.316145,  0.316145],
               [ 0.25    ,  0.25    ,  0.25    ,  0.25    ],
               [ 0.367711,  0.      ,  0.316145,  0.316145]]), rtol=1e-5)

class TestHealpixMap(unittest.TestCase):
    def setUp(self):
        self.hpm = h.HealpixMap(32, 'RING')
    def test_set(self):
        self.hpm[0] = 1
        self.assertEqual(self.hpm[0], 1)
        self.hpm[0,0] = 1
        self.assertEqual(self.hpm[0,0], 1)
        self.hpm[0,1,0] = 1
        self.assertEqual(self.hpm[0,1,0], 1)
    def test_set_px(self):
        c = n.arange(15)
        self.assertEqual(len(self.hpm[c]), len(c))
        self.hpm[c] += c
        n.testing.assert_allclose(self.hpm[c], c)
        self.hpm[c] += c
        n.testing.assert_allclose(self.hpm[c], 2*c)
        self.hpm[c] = c
        n.testing.assert_allclose(self.hpm[c], c)
    def test_set_ang(self):
        th,ph = n.linspace(0,n.pi,10), n.linspace(-n.pi,n.pi,10)
        self.assertEqual(len(self.hpm[th,ph]), len(th))
        self.hpm[th,ph] += th
        n.testing.assert_allclose(self.hpm[th,ph], th)
        self.hpm[th,ph] += th
        n.testing.assert_allclose(self.hpm[th,ph], 2*th)
        self.hpm[th,ph] = th
        n.testing.assert_allclose(self.hpm[th,ph], th)
    def test_set_vec(self):
        x,y = n.linspace(-.5,.5,10), n.linspace(-.5,.5,10)
        z = 1 - n.sqrt(x**2 + y**2)
        self.assertEqual(len(self.hpm[x,y,z]), len(x))
        self.hpm[x,y,z] += x
        n.testing.assert_allclose(self.hpm[x,y,z], x)
        self.hpm[x,y,z] += x
        n.testing.assert_allclose(self.hpm[x,y,z], 2*x)
        self.hpm[x,y,z] = x
        n.testing.assert_allclose(self.hpm[x,y,z], x)
        

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
