# -*- coding: utf-8 -*-
import unittest, ephem, random
import aipy as a, numpy as n

class TestPointingError(unittest.TestCase):
    def setUp(self):
        self.pnterr = a.phs.PointingError('Error String')
    def teststring(self):
        """Test setting a PointingError message string"""
        self.assertEqual(str(self.pnterr), 'Error String')
    def testraise(self):
        """Test raising a PointingError"""
        def raise_pnterr(): raise self.pnterr
        self.assertRaises(a.phs.PointingError, raise_pnterr)

class TestJulDates(unittest.TestCase):
    def setUp(self):
        self.ephemzero = ephem.date('1899/12/31 12:00')
        self.jdzero = 2415020.
    def testephemzero(self):
        """Test converting between ephem dates and JD - zero point"""
        self.assertEqual(a.phs.juldate2ephem(self.jdzero), self.ephemzero)
        self.assertEqual(a.phs.ephem2juldate(self.ephemzero), self.jdzero)
    def testrandom(self):
        """Test converting between ephem dates and JD - various"""
        for i in range(10):
            d1 = random.random() * ephem.now()
            d2 = a.phs.juldate2ephem(a.phs.ephem2juldate(d1))
            self.assertAlmostEqual(d1, d2)

class TestRadioBody(unittest.TestCase):
    def test_attributes(self):
        """Test aipy.phs.RadioFixedBody attributes"""
        epoch = ephem.B1950
        s = a.phs.RadioFixedBody('0:00', '0:00', mfreq=.200, name='src1',
            epoch=epoch, ionref=(.0,.0), srcshape=(.003, .005, .6))
        self.assertEqual(s._ra, ephem.hours('0:00'))
        self.assertEqual(s._dec, ephem.degrees('0:00'))
        self.assertEqual(s.mfreq, .200)
        self.assertEqual(s.src_name, 'src1')
        self.assertAlmostEqual(s._epoch, epoch, 3)
        self.assertEqual(s.ionref, [.0, .0])
        self.assertEqual(s.srcshape, [.003, .005, .6])
        self.assertRaises(RuntimeError, lambda: s.ra)
        self.assertRaises(RuntimeError, lambda: s.dec)
        self.assertRaises(AttributeError, lambda: s.map)
        o = a.phs.ArrayLocation(('0:00','0:00'))
        o.set_ephemtime(epoch)
        s.compute(o)
        self.assertEqual(len(s.get_crds('eq', ncrd=2)), 2)
        self.assertEqual(len(s.get_crds('top', ncrd=2)), 2)
        self.assertEqual(len(s.get_crds('eq', ncrd=3)), 3)
        self.assertEqual(len(s.get_crds('top', ncrd=3)), 3)
        self.assertEqual(s.map.shape, (3,3))
    def testephem(self):
        """Test the aipy.phs.RadioFixedBody ephem interface"""
        o = a.phs.ArrayLocation(('0:00','0:00'))
        for epoch in [ephem.B1900, ephem.B1950, ephem.J2000]:
          for ra in n.arange(1*n.pi/8, 2*n.pi, n.pi/8):
            for dec in n.arange(-3*n.pi/8, n.pi/2, n.pi/8):
                s = a.phs.RadioFixedBody(ra, dec, epoch=epoch)
                o.set_ephemtime(epoch)
                s.compute(o)
                self.assertAlmostEqual(s.a_ra, s._ra, 9)
                self.assertAlmostEqual(s.a_dec, s._dec, 9)
    def test_compute(self):
        """Test the aipy.phs.RadioFixedBody ephem calculations"""
        epoch = ephem.J2000
        o = a.phs.ArrayLocation(('0:00','0:00'))
        o.set_ephemtime(epoch)
        diagonal = n.array([[1.,0,0],[0,1,0],[0,0,1]])
        for ra in n.arange(1*n.pi/8, 2*n.pi, n.pi/8):
          for dec in n.arange(-3*n.pi/8, n.pi/2, n.pi/8):
            s = a.phs.RadioFixedBody(ra, dec, epoch=epoch)
            s.compute(o)
            m = a.coord.top2eq_m(o.sidereal_time()-s.ra, s.dec)
            err = n.abs(diagonal  - n.dot(m, s.map)).sum()
            self.assertAlmostEqual(err, 0, 10)
    def test_get_crds(self):
        """Test the aipy.phs.RadioFixedBody calculated coordinates"""
        epoch = ephem.J2000
        o = a.phs.ArrayLocation(('0:00','0:00'))
        o.set_ephemtime(epoch)
        for ra in n.arange(1*n.pi/8, 2*n.pi, n.pi/8):
          for dec in n.arange(-3*n.pi/8, n.pi/2, n.pi/8):
            s = a.phs.RadioFixedBody(ra, dec, epoch=epoch)
            s.compute(o)
            ra, dec = s.get_crds('eq', ncrd=2)
            self.assertEqual(s.ra, ra); self.assertEqual(s.dec, dec)
            az, alt = s.get_crds('top', ncrd=2)
            self.assertEqual(s.az, az); self.assertEqual(s.alt, alt)
            eq = s.get_crds('eq', ncrd=3)
            ra, dec = a.coord.eq2radec(eq)
            self.assertAlmostEqual(s.ra, ra, 10)
            self.assertAlmostEqual(s.dec, dec, 10)
            top = s.get_crds('top', ncrd=3)
            az, alt = a.coord.top2azalt(top)
            self.assertAlmostEqual(s.az, az, 10)
            self.assertAlmostEqual(s.alt, alt, 10)

class TestSrcCatalog(unittest.TestCase):
    def setUp(self):
        src1 = a.phs.RadioFixedBody('1:00', '1:00', name='src1')
        src2 = a.phs.RadioFixedBody('2:00', '2:00', name='src2')
        src3 = a.phs.RadioFixedBody('3:00', '3:00', name='src3')
        self.srcs = [src1, src2, src3]
        self.cat = a.phs.SrcCatalog(self.srcs)
    def test_add_srcs(self):
        """Test adding sources to a aipy.phs.SrcCatalog() catalog"""
        cat = a.phs.SrcCatalog()
        src1b = a.phs.RadioFixedBody('0:00', '0:00', name='src1')
        src4 = a.phs.RadioFixedBody('4:00', '4:00', name='src4')
        cat.add_srcs(self.srcs)
        self.assertEqual(len(cat), 3)
        cat.add_srcs(src1b)
        self.assertEqual(len(cat), 3)
        cat.add_srcs(src1b, src4)
        self.assertEqual(len(cat), 4)
        srclist = [src1b] + self.srcs[1:] + [src4]
        for name, src in zip([s.src_name for s in srclist], srclist):
            self.assertEqual(src, cat[name])
    def test_get_srcs(self):
        """Test retrieving sources from a aipy.phs.SrcCatalog() catalog"""
        self.assertEqual(self.cat.get_srcs('src1','src2','src3'), self.srcs)
        self.assertEqual(self.cat.get_srcs(['src1','src2']), self.srcs[:2])
        self.assertRaises(KeyError, lambda: self.cat.get_srcs('bad'))
    def test_compute(self):
        """Test the ephem interfaces for a aipy.phs.SrcCatalog() catalog"""
        o = ephem.Observer()
        self.cat.compute(o)
        for src in self.cat.values():
            self.assertNotEqual(src.ra, None)
            self.assertNotEqual(src.dec, None)
    def test_get_crds(self):
        """Test coordinates calculated from a aipy.phs.SrcCatalog() catalog"""
        o = ephem.Observer()
        self.cat.compute(o)
        crd1 = self.cat.get_crds('eq', srcs=['src1'])
        self.assertEqual(crd1.shape, (3,1))
        self.assertTrue(n.all(crd1[:,0] == self.srcs[0].get_crds('eq')))
        crd2 = self.cat.get_crds('top', srcs=['src1','src2'])
        self.assertEqual(crd2.shape, (3,2))
    def test_get(self):
        """Test retrieving source attributes from a aipy.phs.SrcCatalog() catalog"""
        mfreq = self.cat.get('mfreq',srcs=['src1'])
        self.assertEqual(mfreq.shape, (1,))
        mfreq = self.cat.get('mfreq',srcs=['src1','src2'])
        self.assertEqual(mfreq.shape, (2,))
        ionrefs = self.cat.get('ionref',srcs=['src1','src2'])
        self.assertEqual(ionrefs.shape, (2,2))
        srcshapes = self.cat.get('srcshape',srcs=['src1','src2'])
        self.assertEqual(srcshapes.shape, (3,2))

class TestBeam(unittest.TestCase):
    def setUp(self):
        self.fq = n.arange(0,1,.1)
        self.bm = a.phs.Beam(self.fq)
    def test_attributes(self):
        """Test accessing aipy.phs.Beam attributes"""
        self.assertTrue(n.all(self.bm.freqs == self.fq))
        self.assertTrue(n.all(self.bm.chans == n.arange(self.fq.size)))
        self.assertTrue(n.all(self.bm.afreqs == self.fq))
    def test_select_chans(self):
        """Test selecting various aipy.phs.Beam channels"""
        chans = n.array([1,2,3])
        self.bm.select_chans(chans)
        self.assertTrue(n.all(self.bm.chans == chans))
        self.assertTrue(n.all(self.bm.afreqs == self.fq.take(chans)))

class TestAntenna(unittest.TestCase):
    def setUp(self):
        fq = n.arange(0,1,.1)
        self.bm = a.phs.Beam(fq)
        self.ant = a.phs.Antenna(1, 2, 3, self.bm, phsoff=[0,1])
    def test_attributes(self):
        """Test accessing aipy.phs.Antenna attributes"""
        self.assertEqual(self.ant.beam, self.bm)
        pos = n.array([1,2,3], n.float64)
        self.assertTrue(n.all(self.ant.pos == pos))
        x,y,z = self.ant
        self.assertEqual(x, pos[0])
        self.assertEqual(y, pos[1])
        self.assertEqual(z, pos[2])
        self.assertTrue(n.all(self.ant + self.ant == pos*2))
        self.assertTrue(n.all(self.ant - self.ant == 0))
        self.assertTrue(n.all(-self.ant == -pos))
    def test_select_chans(self):
        """Test selecting various aipy.phs.Antenna channels"""
        chans = n.array([1,2,3])
        self.ant.select_chans(chans)
        self.assertTrue(n.all(self.bm.chans == chans))
        self.assertTrue(n.all(self.ant.phsoff == 1))
        
class TestArrayLocation(unittest.TestCase):
    def setUp(self):
        self.aa = a.phs.ArrayLocation(('0','0'))
    def test_attributes(self):
        self.assertEqual(self.aa.pressure, 0)
        self.assertEqual(self.aa.lat, 0)
        self.assertEqual(self.aa.long, 0)
        self.assertEqual(self.aa.elev, 0)
        self.assertTrue(n.all(self.aa._eq2zen == a.coord.eq2top_m(0., 0.)))
    def test_set_jultime(self):
        jd = 2454554.9798841
        self.aa.set_jultime(jd)
        self.assertEqual(self.aa.date, self.aa.epoch)
        self.assertEqual(self.aa.date, a.phs.juldate2ephem(jd))
        self.assertAlmostEqual(self.aa.sidereal_time(), 0, 7)
        eq2now_rnd = n.round(self.aa._eq2now, 7)
        self.assertTrue(n.all(eq2now_rnd ==
            n.array([[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]])))
    def test_get_jultime(self):
        self.aa.set_jultime(2454555)
        self.assertEqual(self.aa.get_jultime(), 2454555)

class TestAntennaArray(unittest.TestCase):
    def setUp(self):
        bm = a.phs.Beam(n.arange(0,1,.1))
        a1 = a.phs.Antenna(0,0,0,bm, [0,0])
        a2 = a.phs.Antenna(1,0,0,bm, [0,1])
        a3 = a.phs.Antenna(0,1,0,bm, [0,2])
        a4 = a.phs.Antenna(0,0,1,bm, [0,3])
        self.ants = [a1,a2,a3,a4]
        self.aa = a.phs.AntennaArray(('0','0'), self.ants)
    def test_attributes(self):
        self.assertEqual(len(self.aa), 4)
        for ai, aj in zip(self.aa, self.ants): self.assertEqual(ai, aj)
        for i in range(4): self.assertEqual(len(self.aa[:i]), i)
    def test_select_chans(self):
        chans = n.array([1,2,3])
        self.aa.select_chans(chans)
        for ant in self.aa.ants:
            self.assertTrue(n.all(ant.beam.chans == chans))
    def test_ij2bl(self):
        self.assertEqual(self.aa.ij2bl(0,1), 258)
        self.assertEqual(self.aa.ij2bl(0,2), 259)
        self.assertEqual(self.aa.ij2bl(1,2), 515)
    def test_bl2ij(self):
        self.assertEqual(self.aa.bl2ij(258), (0,1))
        self.assertEqual(self.aa.bl2ij(259), (0,2))
        self.assertEqual(self.aa.bl2ij(515), (1,2))
    def test_get_baseline(self):
        for j, ant in enumerate(self.aa):
            if j in [0,3]:
                self.aa.set_jultime(2454554.9)
                self.assertTrue(n.all(self.aa.get_baseline(0,j,'r') == ant.pos))
                self.assertTrue(n.all(self.aa.get_baseline(0,j,'e') == ant.pos))
            else:
                self.assertTrue(n.all(self.aa.get_baseline(0,j,'r') == ant.pos))
                self.aa.set_jultime(2454554.9798841)
                bl_rnd = n.round(self.aa.get_baseline(0,j,'e'), 7)
                self.assertTrue(n.all(bl_rnd == ant.pos))
                self.aa.set_jultime(2454554.9)
                bl_rnd = n.round(self.aa.get_baseline(0,j,'e'), 7)
                self.assertFalse(n.all(bl_rnd == ant.pos))
        src = a.phs.RadioFixedBody('12:00', '0:00')
        src.compute(self.aa)
        self.assertRaises(a.phs.PointingError, 
            lambda: self.aa.get_baseline(0,1,src))
        for t in n.random.random((10,)):
            self.aa.set_jultime(2454555. + t)
            src = a.phs.RadioFixedBody(self.aa.sidereal_time(), self.aa.lat,
                epoch=self.aa.epoch)
            src.compute(self.aa)
            zbl_rnd = n.round(self.aa.get_baseline(0,1,'z'), 3)
            sbl_rnd = n.round(self.aa.get_baseline(0,1,src), 3)
            self.assertTrue(n.all(zbl_rnd == sbl_rnd))
    def test_get_phs_offset(self):
        self.assertTrue(n.all(self.aa.get_phs_offset(0,0) == 0))
        self.assertTrue(n.all(self.aa.get_phs_offset(0,1) == 1))
        self.assertTrue(n.all(self.aa.get_phs_offset(0,2) == 2))
        self.assertTrue(n.all(self.aa.get_phs_offset(0,3) == 3))
    def test_gen_uvw(self):
        self.aa.select_chans()
        afreqs = self.aa[0].beam.afreqs
        u,v,w = self.aa.gen_uvw(0,1,'z')
        self.assertEqual(u.shape, (1,afreqs.size))
        self.assertTrue(n.all(u == 0*afreqs))
        self.assertTrue(n.all(v == 0*afreqs))
        self.assertTrue(n.all(w == 1*afreqs))
        u,v,w = self.aa.gen_uvw(0,2,'z')
        self.assertTrue(n.all(u == 1*afreqs))
        self.assertTrue(n.all(v == 0*afreqs))
        self.assertTrue(n.all(w == 0*afreqs))
        u,v,w = self.aa.gen_uvw(0,3,'z')
        self.assertTrue(n.all(u == 0*afreqs))
        self.assertTrue(n.all(v == 1*afreqs))
        self.assertTrue(n.all(w == 0*afreqs))
    def test_gen_phs(self):
        self.aa.select_chans([1,2,3])
        afreqs = self.aa[0].beam.afreqs
        for t in n.random.random((40,)):
            self.aa.set_jultime(2454555. + t)
            src = a.phs.RadioFixedBody(self.aa.sidereal_time(), self.aa.lat,
                epoch=self.aa.epoch)
            src.compute(self.aa)
            if t > .5: 
                seq = src.get_crds('eq', ncrd=3)
                if t > .75: seq = n.array([seq,seq]).transpose()
            else: seq = src
            phs = n.round(self.aa.gen_phs(seq, 0, 1, mfreq=.1), 6)
            ans = n.round(n.exp(-1j*2*n.pi*afreqs), 6)
            if t > .75: self.assertEqual(phs.shape, (2,3))
            else: self.assertEqual(phs.shape, (3,))
            self.assertTrue(n.all(phs == ans))
            phs = n.round(self.aa.gen_phs(src, 0, 2, mfreq=.1), 3)
            self.assertTrue(n.all(phs == 1+0j))
            phs = n.round(self.aa.gen_phs(src, 0, 3, mfreq=.1), 3)
            self.assertTrue(n.all(phs == 1+0j))
        phs1 = self.aa.gen_phs(src, 0, 2, mfreq=.1, ionref=(.001,.001))
        phs2 = self.aa.gen_phs(src, 0, 2, mfreq=.1, srcshape=(.01,.01,0), 
            resolve_src=True)
        self.assertTrue(n.all(phs1 != 1+0j))
        self.assertTrue(n.all(phs2 != 1+0j))
    def test_resolve_src(self):
        amp = self.aa.resolve_src(100., 100., srcshape=(0,0,0))
        self.assertEqual(amp, 1)
        amp1 = self.aa.resolve_src(100., 50., srcshape=(.01,0,n.pi/2))
        amp2 = self.aa.resolve_src(100., 50., srcshape=(0,.01,0))
        self.assertAlmostEqual(amp1, amp2, 15)
        amp1 = self.aa.resolve_src(100., 100., srcshape=(.02,0,0))
        amp2 = self.aa.resolve_src(100., 100., srcshape=(.01414,.01414,0))
        self.assertAlmostEqual(amp1, amp2, 3)
        amp = self.aa.resolve_src(100., 0., srcshape=(0.001,0,0))
        x = 2*n.pi * .1
        self.assertEqual(amp, 2*a._cephes.j1(x)/x)
    def test_refract(self):
        self.aa.select_chans([1,2,3])
        afreqs = self.aa[0].beam.afreqs
        zeros = n.zeros((1,3), dtype=n.float)
        ones = n.ones((1,3), dtype=n.float)
        # Test non-vectors, dra->u association
        dw = self.aa.refract(ones, zeros, mfreq=.1, ionref=(.001,0))
        ans = .001 / (afreqs/.1)**2 ; ans.shape = (1, ans.size)
        self.assertEqual(len(dw.shape), 2)
        self.assertTrue(n.all(n.round(dw,10) == n.round(ans,10)))
        # Test non-vectors, no dra-> v association
        dw = self.aa.refract(zeros, ones, mfreq=.1,ionref=(.001,0))
        self.assertTrue(n.all(dw == 0))
        # Test non-vectors, no ddec->u association
        dw = self.aa.refract(ones, zeros, mfreq=.1, ionref=(0,.001))
        self.assertTrue(n.all(dw == 0))
        # Test non-vectors, ddec->v association, v scaling
        dw = self.aa.refract(zeros, 2*ones, mfreq=.1, ionref=(0,.001))
        self.assertTrue(n.all(n.round(dw,10) == n.round(2*ans,10)))
        # Test vectors, mfreq scaling
        ones = n.ones((2,3), dtype=n.float)
        ionref = (n.array([0, .001]), n.array([.001, 0]))
        mfreq = n.array([.1, .2])
        ans = n.array([.002/(afreqs / .1)**2, .002/(afreqs / .2)**2])
        dw = self.aa.refract(2*ones, 2*ones, mfreq=mfreq, ionref=ionref)
        self.assertTrue(n.all(n.round(dw,10) == n.round(ans,10)))
    def test_phs2src(self):
        self.aa.select_chans([1,2,3])
        self.aa.set_jultime(2454555.)
        src = a.phs.RadioFixedBody('0:00', '20:00')
        src.compute(self.aa)
        self.assertTrue(n.all(self.aa.phs2src(1.,src,0,1) == \
            self.aa.gen_phs(src,0,1)))
    def test_unphs2src(self):
        self.aa.select_chans([1,2,3])
        self.aa.set_jultime(2454555.)
        src = a.phs.RadioFixedBody('0:00', '20:00')
        src.compute(self.aa)
        self.assertTrue(n.all(
            self.aa.unphs2src(self.aa.gen_phs(src,0,1),src,0,1) == 1.))

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.phs unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestPointingError))
        self.addTests(loader.loadTestsFromTestCase(TestJulDates))
        self.addTests(loader.loadTestsFromTestCase(TestRadioBody))
        self.addTests(loader.loadTestsFromTestCase(TestSrcCatalog))
        self.addTests(loader.loadTestsFromTestCase(TestBeam))
        self.addTests(loader.loadTestsFromTestCase(TestAntenna))
        self.addTests(loader.loadTestsFromTestCase(TestArrayLocation))
        self.addTests(loader.loadTestsFromTestCase(TestAntennaArray))

if __name__ == '__main__':
    unittest.main()
