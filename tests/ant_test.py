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

'''
#! /usr/bin/env python
import aipy as a, numpy as n, ephem as e

JULDATE = 2454000
epoch = a.ant.juldate2ephem(JULDATE)
#epoch = e.J2000

src1 = a.ant.RadioFixedBody('12:00','45:00', epoch=epoch, name='src1')
src2 = a.ant.RadioFixedBody('12:00','90:00', epoch=epoch, name='src2')
sun = a.ant.RadioSpecial('Sun')
cat = a.ant.SrcCatalog([src1,src2])

bm = a.ant.Beam(n.array([.1,.2,.3]), active_chans=n.array([0,1]))
ant1 = a.ant.Antenna(0, 0, 0, bm)
ant2 = a.ant.Antenna(100, 0, 0, bm)
aa = a.ant.AntennaArray(('0:00','0:00'), [ant1,ant2])

aa.set_jultime(JULDATE)
sun.compute(aa)
cat.compute(aa)

print '---------------------------------------------------------------'
print '(0,1) eq raw:', aa.get_baseline(0,1,src='r')
print '(0,1) eq now:', aa.get_baseline(0,1,src='e')
print '(0,1) top:', aa.get_baseline(0,1,src='z')
print '---------------------------------------------------------------'
print 'Sun at RA=%s, DEC=%s' % (sun.ra, sun.dec)
print '       AZ=%s, ALT=%s' % (sun.az, sun.alt)
print '(0,1) Sun:', aa.get_baseline(0,1,src=sun)
print '---------------------------------------------------------------'
print 'src1 at RA=%s, DEC=%s' % (cat['src1'].ra, cat['src1'].dec)
print '        AZ=%s, ALT=%s' % (cat['src1'].az, cat['src1'].alt)
print '(0,1) src1:', aa.get_baseline(0,1,src=cat['src1'])
print '---------------------------------------------------------------'
print 'src2 at RA=%s, DEC=%s' % (cat['src2'].ra, cat['src2'].dec)
print '        AZ=%s, ALT=%s' % (cat['src2'].az, cat['src2'].alt)
print '(0,1) src2:'
x,y,z = aa.get_baseline(0,1,src=cat.get_crds('eq', srcs=['src2']))
print 'x=%s\ny=%s\nz=%s' % (x,y,z)
print '---------------------------------------------------------------'
print '(0,1) %s:' % cat.keys()
x,y,z = aa.get_baseline(0,1,src=cat.get_crds('eq'))
print 'x=%s\ny=%s\nz=%s' % (x,y,z)
print '---------------------------------------------------------------'
print 'uvw Sun for freqs %s:' % (bm.afreqs)
u,v,w = aa.gen_uvw(0,1, src=sun)
print 'u=%s\nv=%s\nw=%s' % (u,v,w)
print '---------------------------------------------------------------'
print 'uvw %s:'% cat.keys()
u,v,w = aa.gen_uvw(0,1, src=cat.get_crds('eq'))
print 'u=%s\nv=%s\nw=%s' % (u,v,w)
print '---------------------------------------------------------------'
print aa.resolve_src(cat.get_crds('eq', srcs=['src2']), 0,1, 
    srcshape=(.02,.01,0.))
print aa.resolve_src(cat['src2'], 0, 1, srcshape=(.02,.01,n.pi/2))
print aa.resolve_src(cat['src2'], 0, 1, srcshape=(.01,.02,0.))
print aa.resolve_src(cat.get_crds('eq'), 0, 1, srcshape=(.01,.02,0.))
'''
