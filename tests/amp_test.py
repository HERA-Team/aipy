# -*- coding: utf-8 -*-
import unittest
import aipy.amp as amp, numpy as n

class TestRadioBody(unittest.TestCase):
    def setUp(self):
        self.fqs = n.arange(.1,.2,.01)
        bm = amp.Beam(self.fqs)
        ant0 = amp.Antenna(0,0,0,bm)
        self.aa = amp.AntennaArray(('0','0'), [ant0])
    def test_attributes(self):
        """Test aipy.amp.RadioFixedBody attributes"""
        s = amp.RadioFixedBody('0:00', '0:00',
            jys=100, index=-2, mfreq=.1, name='src1')
        self.assertEqual(s._jys, 100)
        self.assertEqual(s.index, -2)
        s.compute(self.aa)
        self.assertTrue(n.all(s.get_jys() == 100 * (self.fqs / .1)**-2))

class TestBeam(unittest.TestCase):
    def setUp(self):
        self.fqs = n.arange(.1,.2,.01)
        self.bm = amp.Beam(self.fqs)
    def test_response(self):
        """Test retrieving aipy.amp.Beam response"""
        xyz = (0,0,1)
        self.assertTrue(n.all(self.bm.response(xyz) == n.ones_like(self.fqs)))
        x = n.array([0, .1, .2])
        y = n.array([.1, .2, 0])
        z = n.array([.2, 0, .1])
        xyz = (x,y,z)
        self.assertTrue(n.all(self.bm.response(xyz) == \
            n.ones((self.fqs.size,3))))
        self.bm.select_chans([0,1,2])
        self.assertTrue(n.all(self.bm.response(xyz) == n.ones((3,3))))

class TestBeam2DGaussian(unittest.TestCase):
    def setUp(self):
        self.fqs = n.arange(.1,.2,.01)
        self.bm = amp.Beam2DGaussian(self.fqs, .05, .025)
    def test_response(self):
        """Test retrieving a 2D Gaussian beam response"""
        xyz = (0,0,1)
        self.assertTrue(n.all(self.bm.response(xyz) == n.ones_like(self.fqs)))
        x = n.array([0, .05,   0])
        y = n.array([0,   0, .05])
        z = n.array([1,   1,   1])
        xyz = (x,y,z)
        resp = self.bm.response(xyz)
        self.assertEqual(resp.shape, (self.fqs.size,3))
        ans = n.sqrt(n.array([1., n.exp(-0.5), n.exp(-2)]))
        ans.shape = (1,3)
        self.bm.select_chans([0])
        resp = self.bm.response(xyz)
        self.assertTrue(n.all(n.round(resp - ans, 3) == 0))

class TestBeamAlm(unittest.TestCase):
    def setUp(self):
        self.fqs = n.arange(.1,.2,.01)
        self.coeffs = {0:n.array([1.,0,0], dtype=n.complex)}
        self.bm = amp.BeamAlm(self.fqs, lmax=1, mmax=1, deg=0, coeffs=self.coeffs)
    def test_init_update(self):
        self.assertTrue(n.allclose(self.bm.hmap[0].map, 1./n.sqrt(4*n.pi)))
    def test_response(self):
        x = n.array([0, .05,   0])
        y = n.array([0,   0, .05])
        z = n.array([1,   1,   1])
        xyz = (x,y,z)
        resp = self.bm.response(xyz)
        self.assertEqual(resp.shape, (self.fqs.size,3))
        ans = n.ones(3, dtype=n.complex) / n.sqrt(4*n.pi)
        ans.shape = (1,3)
        self.bm.select_chans([0])
        resp = self.bm.response(xyz)
        self.assertTrue(n.all(n.round(resp - ans, 3) == 0))

class TestAntenna(unittest.TestCase):
    def setUp(self):
        self.fqs = n.arange(.1,.2,.01)
        bm = amp.Beam2DGaussian(self.fqs, .05, .025)
        self.ant = amp.Antenna(0,0,0, beam=bm)
    def test_passband(self):
        """Test the Antenna passband"""
        pb = self.ant.passband()
        self.assertTrue(n.all(pb == n.ones_like(self.fqs)))
        self.ant.select_chans([0,1,2])
        pb = self.ant.passband()
        self.assertEqual(pb.shape, (3,))
    def test_bm_response(self):
        """Test the Antenna beam response"""
        xyz = (.05,0,1)
        self.ant.select_chans([0])
        # resp = self.ant.bm_response(xyz, pol='x')
        # self.assertAlmostEqual(resp, n.sqrt(n.exp(-0.5)), 3)
        # resp = self.ant.bm_response(xyz, pol='y')
        # self.assertAlmostEqual(resp, n.sqrt(n.exp(-2)), 3)

#class TestMemLeaks(unittest.TestCase):
#    def test_antenna_create(self):
#        freqs = n.arange(.1,.2,.001)
#        beam = amp.Beam(freqs)
#        while True: ant = amp.Antenna(0,0,0,beam, pointing=(0,n.pi/2,.1))
#    def test_aa_create(self):
#        freqs = n.arange(.1,.2,.001)
#        beam = amp.Beam(freqs)
#        ants = [amp.Antenna(0,0,0,beam) for i in range(100)]
#        while True: aa = amp.AntennaArray(('0:00','0:00'), ants)
        
class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.amp unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestRadioBody))
        self.addTests(loader.loadTestsFromTestCase(TestBeam))
        self.addTests(loader.loadTestsFromTestCase(TestBeam2DGaussian))
        self.addTests(loader.loadTestsFromTestCase(TestAntenna))
        #self.addTests(loader.loadTestsFromTestCase(TestMemLeaks))

if __name__ == '__main__':
    unittest.main()
