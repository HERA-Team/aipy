import unittest, aipy._src.helm as h, aipy as a, numpy as n

class TestHelmboldtFixedBody(unittest.TestCase):
    def setUp(self):
        self.s1 = h.HelmboldtFixedBody('0:00','0:00',name='test1',
            jys=100.,index=[-1],mfreq=.074)
        self.s2 = h.HelmboldtFixedBody('0:00','0:00',name='test2',
            jys=100.,index=[-1,-4,-1],mfreq=.074)
    def test_init(self):
        self.assertEqual(self.s1._jys, 100.)
        self.assertEqual(self.s1.index, [-1])
        self.assertEqual(self.s2._jys, 100.)
        self.assertEqual(self.s2.index, [-1,-4,-1])
    def test_compute(self):
        fq = n.array([.074, .150, .300])
        X = n.log10(fq / .074)
        bm = a.phs.Beam(fq)
        ant = a.phs.Antenna(0,0,0,bm)
        aa = a.phs.AntennaArray(('0:00','0:00'),[ant])
        self.s1.compute(aa)
        B = self.s1.index[0]
        self.assertTrue(n.all(n.round(self.s1.get_jys() - 100*10**(B*X), 10) == 0))
        self.s2.compute(aa)
        B,C,D = self.s2.index
        self.assertTrue(n.all(n.round(self.s2.get_jys() - 100*10**(B*X+C*n.exp(D*X)), 10) == 0))
    def test_get_params(self):
        prms = self.s1.get_params(['jys', 'index'])
        self.assertEqual(prms['jys'], 100.)
        self.assertEqual(prms['index'], [-1])
        prms = self.s2.get_params(['jys', 'index'])
        self.assertEqual(prms['jys'], 100.)
        self.assertEqual(prms['index'], [-1,-4,-1])
    def test_set_params(self):
        s1 = h.HelmboldtFixedBody('0:00','0:00',name='test1',
            jys=100.,index=[-1],mfreq=.074)
        s1.set_params({'index':[-2]})
        self.assertEqual(s1.index, [-2])
        s1.set_params({'index':[-1,-4,-1]})
        self.assertEqual(s1.index, [-1,-4,-1])

class TestHelmboldtCatalog(unittest.TestCase):
    def setUp(self):
        self.cat = h.HelmboldtCatalog()
        self.cat.fromfile(h.POSFILE, h.FITFILE)
    def test_crosscheck(self):
        srcs = {}
        srclines = [L for L in open(h.POSFILE).readlines() if L.startswith('J')]
        for line in srclines:
            s = line[:9]
            srcs[s] = srcs.get(s, []) + [map(float, line[58:73].split())]
        for s in srcs:
            if not self.cat.has_key(s): continue
            d = n.array(srcs[s])
            fq = d[:,0] / 1e3
            flx = d[:,1]
            bm = a.phs.Beam(fq)
            ant = a.phs.Antenna(0,0,0,bm)
            aa = a.phs.AntennaArray(('0:00','0:00'),[ant])
            self.cat[s].compute(aa)
            jys = self.cat[s].get_jys()
            ratio = jys / flx
            bad = n.logical_or(ratio > 1.5, ratio < .5).sum()
            if isinstance(self.cat[s], h.HelmboldtFixedBody) and s != 'J0320+413':
                self.assertTrue(bad < len(ratio)/2.)
    def test_get_metadata(self):
        md = self.cat.get_metadata()
        self.assertEqual(md['J0000+554'][0], (.074, 15.92, 0.37))    
        self.assertEqual(md['J0000+554'][1], (.325,  5.78, 0.01))    
        self.assertEqual(md['J2359+440'][-1], (14.900,  0.12, 0.01))    
        


if __name__ == '__main__':
    unittest.main()
