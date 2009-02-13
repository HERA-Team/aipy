import unittest, ephem, random
import aipy as a, numpy as n

class TestPointingError(unittest.TestCase):
    def setUp(self):
        self.pnterr = a.ant.PointingError('Error String')
    def teststring(self):
        self.assertEqual(str(self.pnterr), 'Error String')
    def testraise(self):
        def raise_pnterr(): raise self.pnterr
        self.assertRaises(a.ant.PointingError, raise_pnterr)

class TestJulDates(unittest.TestCase):
    def setUp(self):
        self.ephemzero = ephem.date('1899/12/31 12:00')
        self.jdzero = 2415020.
    def testephemzero(self):
        self.assertEqual(a.ant.juldate2ephem(self.jdzero), self.ephemzero)
        self.assertEqual(a.ant.ephem2juldate(self.ephemzero), self.jdzero)
    def testrandom(self):
        for i in range(10):
            d1 = random.random() * ephem.now()
            d2 = a.ant.juldate2ephem(a.ant.ephem2juldate(d1))
            self.assertAlmostEqual(d1, d2)

class TestRadioBody(unittest.TestCase):
    def test_attributes(self):
        epoch = ephem.B1950
        s = a.ant.RadioFixedBody('0:00', '0:00', mfreq=.200, name='src1',
            epoch=epoch, ionref=(.0,.0), srcshape=(.003, .005, .6))
        self.assertEqual(s._ra, ephem.hours('0:00'))
        self.assertEqual(s._dec, ephem.degrees('0:00'))
        self.assertEqual(s.mfreq, .200)
        self.assertEqual(s.src_name, 'src1')
        #self.assertEqual(s._epoch, epoch)
        self.assertEqual(s.ionref, [.0, .0])
        self.assertEqual(s.srcshape, [.003, .005, .6])
        self.assertRaises(RuntimeError, lambda: s.ra)
        self.assertRaises(RuntimeError, lambda: s.dec)
        self.assertRaises(AttributeError, lambda: s.map)
        o = a.ant.ArrayLocation(('0:00','0:00'))
        o.set_ephemtime(epoch)
        s.compute(o)
        self.assertEqual(len(s.get_crd('eq', ncrd=2)), 2)
        self.assertEqual(len(s.get_crd('top', ncrd=2)), 2)
        self.assertEqual(len(s.get_crd('eq', ncrd=3)), 3)
        self.assertEqual(len(s.get_crd('top', ncrd=3)), 3)
        self.assertEqual(s.map.shape, (3,3))
    #def testephem(self):
    #    o = a.ant.ArrayLocation(('0:00','0:00'))
    #    for epoch in [ephem.B1900, ephem.B1950, ephem.J2000]:
    #      for ra in n.arange(0, 2*n.pi, n.pi/4):
    #        for dec in n.arange(-n.pi/2, n.pi/2, n.pi/4):
    #            s = a.ant.RadioFixedBody(ra, dec, epoch=epoch)
    #            o.set_ephemtime(epoch)
    #            s.compute(o)
    #            self.assertEqual(s.ra, s._ra)
    #            self.assertEqual(s.dec, s._dec)
    def test_compute(self):
        epoch = ephem.J2000
        o = a.ant.ArrayLocation(('0:00','0:00'))
        o.set_ephemtime(epoch)
        diagonal = n.array([[1.,0,0],[0,1,0],[0,0,1]])
        for ra in n.arange(1*n.pi/8, 2*n.pi, n.pi/8):
          for dec in n.arange(-3*n.pi/8, n.pi/2, n.pi/8):
            s = a.ant.RadioFixedBody(ra, dec, epoch=epoch)
            s.compute(o)
            m = a.coord.top2eq_m(o.sidereal_time()-s.ra, s.dec)
            err = n.abs(diagonal  - n.dot(m, s.map)).sum()
            self.assertAlmostEqual(err, 0, 10)
    def test_get_crd(self):
        epoch = ephem.J2000
        o = a.ant.ArrayLocation(('0:00','0:00'))
        o.set_ephemtime(epoch)
        for ra in n.arange(1*n.pi/8, 2*n.pi, n.pi/8):
          for dec in n.arange(-3*n.pi/8, n.pi/2, n.pi/8):
            s = a.ant.RadioFixedBody(ra, dec, epoch=epoch)
            s.compute(o)
            ra, dec = s.get_crd('eq', ncrd=2)
            self.assertEqual(s.ra, ra); self.assertEqual(s.dec, dec)
            az, alt = s.get_crd('top', ncrd=2)
            self.assertEqual(s.az, az); self.assertEqual(s.alt, alt)
            eq = s.get_crd('eq', ncrd=3)
            ra, dec = a.coord.eq2radec(eq)
            self.assertAlmostEqual(s.ra, ra, 10)
            self.assertAlmostEqual(s.dec, dec, 10)
            top = s.get_crd('top', ncrd=3)
            az, alt = a.coord.top2azalt(top)
            self.assertAlmostEqual(s.az, az, 10)
            self.assertAlmostEqual(s.alt, alt, 10)
    

if __name__ == '__main__':
    unittest.main()
