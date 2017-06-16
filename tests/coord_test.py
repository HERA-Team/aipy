# -*- coding: utf-8 -*-
import unittest, ephem as e, random
import aipy as a, numpy as n

class TestConvert(unittest.TestCase):
    def testtypecheck(self):
        """Test coordinate conversion types"""
        crd2, crd3 = (0, 0), (1, 0, 0)
        for sys in a.coord.sys_dict:
            ans2 = a.coord.convert(crd2, sys, sys)
            self.assertEqual(crd2[0], ans2[0])
            self.assertEqual(crd2[1], ans2[1])
            ans3 = a.coord.convert(crd3, sys, sys)
            self.assertEqual(crd2[0], ans3[0])
            self.assertEqual(crd2[1], ans3[1])
        self.assertRaises(KeyError, a.coord.convert, crd2, 'bad','bad')
    def testprecession(self):
        """Test coordinate precessions for accuracy"""
        crdpairs = [[('0:00','0:00'), ('00:02:33.77','00:16:42.1')],
            [('6:00','0:00'), ('06:02:33.75','-00:00:05.6')],
            [('0:00','89:00'), ('00:03:03.75','+89:16:41.7')]]
        for b1950,j2000 in crdpairs:
            c1 = a.coord.convert(b1950,'eq','eq',iepoch=e.B1950,oepoch=e.J2000)
            c2 = a.coord.convert(j2000,'eq','eq',iepoch=e.J2000,oepoch=e.B1950)
            c1_ck = e.Equatorial(j2000[0],j2000[1],epoch=e.J2000).get()
            c2_ck = e.Equatorial(b1950[0],b1950[1],epoch=e.B1950).get()
            self.assertAlmostEqual(c1[0], c1_ck[0], 4)
            self.assertAlmostEqual(c1[1], c1_ck[1], 4)
            self.assertAlmostEqual(c2[0], c2_ck[0], 4)
            self.assertAlmostEqual(c2[1], c2_ck[1], 4)
    def testcrdsys(self):
        """Test coordinate conversions for accuracy"""
        eq_ec_ga = [
            [   ('19:59:28.3566','40:44:02.096'),
                ('317.7054323','59.3254895'),
                ('76.1898379','5.7554756')],
            [   ('12:30:49.4233','12:23:28.043'),
                ('182.0592608','14.4166861'),
                ('283.7777978','74.4911308')],
            [   ('13:25:27.6152','-43:01:08.805'),
                ('217.1433477','-31.3319020'),
                ('309.5158743','19.4173247')]]
        for eq,ec,ga in eq_ec_ga:
            eq_ec = a.coord.convert(eq,'eq','ec')
            eq_ga = a.coord.convert(eq,'eq','ga')
            ec_eq = a.coord.convert(ec,'ec','eq')
            ec_ga = a.coord.convert(ec,'ec','ga')
            ga_eq = a.coord.convert(ga,'ga','eq')
            ga_ec = a.coord.convert(ga,'ga','ec')
            eq_ck = e.Equatorial(eq[0], eq[1], epoch=e.J2000).get()
            ec_ck = e.Ecliptic(ec[0], ec[1], epoch=e.J2000).get()
            ga_ck = e.Galactic(ga[0], ga[1], epoch=e.J2000).get()
            self.assertAlmostEqual(eq_ec[0], ec_ck[0], 4)
            self.assertAlmostEqual(eq_ec[1], ec_ck[1], 4)
            self.assertAlmostEqual(eq_ga[0], ga_ck[0], 4)
            self.assertAlmostEqual(eq_ga[1], ga_ck[1], 4)
            self.assertAlmostEqual(ec_eq[0], eq_ck[0], 4)
            self.assertAlmostEqual(ec_eq[1], eq_ck[1], 4)
            self.assertAlmostEqual(ec_ga[0], ga_ck[0], 4)
            self.assertAlmostEqual(ec_ga[1], ga_ck[1], 4)
            self.assertAlmostEqual(ga_eq[0], eq_ck[0], 4)
            self.assertAlmostEqual(ga_eq[1], eq_ck[1], 4)
            self.assertAlmostEqual(ga_ec[0], ec_ck[0], 4)
            self.assertAlmostEqual(ga_ec[1], ec_ck[1], 4)

class TestConvert_m(unittest.TestCase):
    def testshape(self):
        """Test conversion matrices for shape preservation"""
        for c1 in a.coord.sys_dict:
            for c2 in a.coord.sys_dict:
                m = a.coord.convert_m(c1,c2)
                self.assertEqual(len(m.shape), 2)
                self.assertEqual(m.shape[0], 3)
                self.assertEqual(m.shape[1], 3)
    def testdiagonal(self):
        """Test conversion matrices for normalcy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        for c1 in a.coord.sys_dict:
            for c2 in a.coord.sys_dict:
                m1 = a.coord.convert_m(c1,c2)
                m2 = a.coord.convert_m(c2,c1)
                m = n.round(n.dot(m1, m2), 10)
                self.assertTrue(n.all(m == diag))

class TestRot_m(unittest.TestCase):
    def testshape(self):
        """Test rotation matrices for shape preservation"""
        m = a.coord.rot_m(0, n.array([1.,0,0]))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 3)
        self.assertEqual(m.shape[1], 3)
        m = a.coord.rot_m(n.array([0., 0.]), n.array([[1.,0,0],[1.,0,0]]))
        self.assertEqual(len(m.shape), 3)
        self.assertEqual(m.shape[0], 2)
        self.assertEqual(m.shape[1], 3)
        self.assertEqual(m.shape[2], 3)
    def testaccuracy(self):
        """Test rotation matricies for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.rot_m(2*n.pi, n.array([1.,0,0])), 10)
        self.assertTrue(n.all(m == diag))
        m = n.round(a.coord.rot_m(n.pi/2, n.array([0.,0,1])), 10)
        self.assertTrue(n.all(m == n.array([-e2,e1,e3])))
        m = n.round(a.coord.rot_m(n.pi/2, n.array([0,1,0])), 10)
        self.assertTrue(n.all(m == n.array([e3,e2,-e1])))
        m = n.round(a.coord.rot_m(n.pi/2, n.array([1,0,0])), 10)
        self.assertTrue(n.all(m == n.array([e1,-e3,e2])))
        
class TestXyz2thphi(unittest.TestCase):
    def testshape(self):
        """Test the x,y,z to theta,phi conversion for shape preservation"""
        m = a.coord.xyz2thphi((0,0,1))
        self.assertEqual(len(m.shape), 1)
        self.assertEqual(m.shape[0], 2)
        m = a.coord.xyz2thphi((n.array([0.,0,0]),
            n.array([0.,0,0]),n.array([1.,1,1])))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 2)
        self.assertEqual(m.shape[1], 3)
    def testaccuracy(self):
        """Test the x,y,z to theta,phi conversion for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.xyz2thphi(e3), 10)
        self.assertTrue(n.all(m == n.array([0.,0])))
        m = n.round(a.coord.xyz2thphi(e1), 10)
        self.assertTrue(n.all(m == n.round(n.array([n.pi/2,0]), 10)))
        m = n.round(a.coord.xyz2thphi(e2), 10)
        self.assertTrue(n.all(m == n.round(n.array([n.pi/2,n.pi/2]), 10)))
        
class TestThphi2xyz(unittest.TestCase):
    def testshape(self):
        """Test the theta,phi to x,y,z conversion for shape preservation"""
        m = a.coord.thphi2xyz((0,0))
        self.assertEqual(len(m.shape), 1)
        self.assertEqual(m.shape[0], 3)
        m = a.coord.thphi2xyz((n.array([0.,0]),n.array([0.,0])))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 3)
        self.assertEqual(m.shape[1], 2)
    def testaccuracy(self):
        """Test the theta,phi to x,y,z conversion for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.thphi2xyz((0,0)), 10)
        self.assertTrue(n.all(m == e3))
        m = n.round(a.coord.thphi2xyz((n.pi/2,0)), 10)
        self.assertTrue(n.all(m == e1))
        m = n.round(a.coord.thphi2xyz((n.pi/2,n.pi/2)), 10)
        self.assertTrue(n.all(m == e2))
        
class TestEq2radec(unittest.TestCase):
    def testshape(self):
        """Test the equatorial to ra,dec conversion for shape preservation"""
        m = a.coord.eq2radec((0,0,1))
        self.assertEqual(len(m.shape), 1)
        self.assertEqual(m.shape[0], 2)
        m = a.coord.eq2radec((n.array([0.,0,0]),
            n.array([0.,0,0]),n.array([1.,1,1])))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 2)
        self.assertEqual(m.shape[1], 3)
    def testaccuracy(self):
        """Test the equatorial to ra,dec conversion for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.eq2radec(e3), 10)
        self.assertTrue(n.all(m == n.round(n.array([0.,n.pi/2]), 10)))
        m = n.round(a.coord.eq2radec(e1), 10)
        self.assertTrue(n.all(m == n.round(n.array([0,0]), 10)))
        m = n.round(a.coord.eq2radec(e2), 10)
        self.assertTrue(n.all(m == n.round(n.array([n.pi/2,0]), 10)))

class TestRadec2eq(unittest.TestCase):
    def testshape(self):
        """Test the ra,dec to equatorial conversion for shape preservation"""
        m = a.coord.radec2eq((0,0))
        self.assertEqual(len(m.shape), 1)
        self.assertEqual(m.shape[0], 3)
        m = a.coord.radec2eq((n.array([0.,0]),n.array([0.,0])))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 3)
        self.assertEqual(m.shape[1], 2)
    def testaccuracy(self):
        """Test the ra,dec to equatorial conversion for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.radec2eq((0,0)), 10)
        self.assertTrue(n.all(m == e1))
        m = n.round(a.coord.radec2eq((n.pi/2,0)), 10)
        self.assertTrue(n.all(m == e2))
        m = n.round(a.coord.radec2eq((n.pi/2,n.pi/2)), 10)
        self.assertTrue(n.all(m == e3))

class TestLatlong2xyz(unittest.TestCase):
    def testshape(self):
        """Test the lat,long to x,y,z conversion for shape preservation"""
        m = a.coord.latlong2xyz((0,0))
        self.assertEqual(len(m.shape), 1)
        self.assertEqual(m.shape[0], 3)
        m = a.coord.latlong2xyz((n.array([0.,0]),n.array([0.,0])))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 3)
        self.assertEqual(m.shape[1], 2)
    def testaccuracy(self):
        """Test the lat,long to x,y,z conversion for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.latlong2xyz((0,0)), 10)
        self.assertTrue(n.all(m == e1))
        m = n.round(a.coord.latlong2xyz((n.pi/2,0)), 10)
        self.assertTrue(n.all(m == e3))
        m = n.round(a.coord.latlong2xyz((0,n.pi/2)), 10)
        self.assertTrue(n.all(m == e2))

class TestTop2azalt(unittest.TestCase):
    def testshape(self):
        """Test the x,y,z to az,alt conversion for shape preservation"""
        m = a.coord.top2azalt((0,0,1))
        self.assertEqual(len(m.shape), 1)
        self.assertEqual(m.shape[0], 2)
        m = a.coord.top2azalt((n.array([0.,0,0]),
            n.array([0.,0,0]),n.array([1.,1,1])))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 2)
        self.assertEqual(m.shape[1], 3)
    def testaccuracy(self):
        """Test the x,y,z to az,alt conversion for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = a.coord.top2azalt(e3)
        self.assertAlmostEqual(m[1], n.pi/2, 10)
        m = n.round(a.coord.top2azalt(e1), 10)
        self.assertTrue(n.all(m == n.round(n.array([n.pi/2,0]), 10)))
        m = n.round(a.coord.top2azalt(e2), 10)
        self.assertTrue(n.all(m == n.round(n.array([0,0]), 10)))

class TestAzalt2top(unittest.TestCase):
    def testshape(self):
        """Test the az,alt to x,y,z conversion for shape preservation"""
        m = a.coord.azalt2top((0,0))
        self.assertEqual(len(m.shape), 1)
        self.assertEqual(m.shape[0], 3)
        m = a.coord.azalt2top((n.array([0.,0]),n.array([0.,0])))
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 3)
        self.assertEqual(m.shape[1], 2)
    def testaccuracy(self):
        """Test the az,alt to x,y,z conversion for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.azalt2top((0,0)), 10)
        self.assertTrue(n.all(m == e2))
        m = n.round(a.coord.azalt2top((n.pi/2,0)), 10)
        self.assertTrue(n.all(m == e1))
        m = n.round(a.coord.azalt2top((0,n.pi/2)), 10)
        self.assertTrue(n.all(m == e3))

class TestEq2top_m(unittest.TestCase):
    def testshape(self):
        """Test the equatorial/x,y,z rotation matrix for shape preservation"""
        m = a.coord.eq2top_m(0, 0)
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 3)
        self.assertEqual(m.shape[1], 3)
        m = a.coord.eq2top_m(n.array([0., 0]), n.array([0.,0]))
        self.assertEqual(len(m.shape), 3)
        self.assertEqual(m.shape[0], 2)
        self.assertEqual(m.shape[1], 3)
        self.assertEqual(m.shape[2], 3)
    def testaccuracy(self):
        """Test the equatorial/x,y,z rotation matrix for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.eq2top_m(0, 0), 10)
        self.assertTrue(n.all(m == n.array([e2,e3,e1])))
        m = n.round(a.coord.eq2top_m(-n.pi/2, 0.), 10)
        self.assertTrue(n.all(m == n.array([-e1,e3,e2])))
        m = n.round(a.coord.eq2top_m(0, n.pi/2), 10)
        self.assertTrue(n.all(m == n.array([e2,-e1,e3])))
        
class TestTop2eq_m(unittest.TestCase):
    def testshape(self):
        """Test the x,y,z/equatorial rotation matrix for shape preservation"""
        m = a.coord.top2eq_m(0, 0)
        self.assertEqual(len(m.shape), 2)
        self.assertEqual(m.shape[0], 3)
        self.assertEqual(m.shape[1], 3)
        m = a.coord.top2eq_m(n.array([0., 0]), n.array([0.,0]))
        self.assertEqual(len(m.shape), 3)
        self.assertEqual(m.shape[0], 2)
        self.assertEqual(m.shape[1], 3)
        self.assertEqual(m.shape[2], 3)
    def testaccuracy(self):
        """Test the x,y,z/equatorial rotation matrix for accuracy"""
        diag = n.array([[1,0,0],[0,1,0],[0,0,1]], dtype=n.double)
        e1,e2,e3 = diag
        m = n.round(a.coord.top2eq_m(0, 0), 10)
        self.assertTrue(n.all(m == n.array([e3,e1,e2])))
        m = n.round(a.coord.top2eq_m(-n.pi/2, 0.), 10)
        self.assertTrue(n.all(m == n.array([-e1,e3,e2])))
        m = n.round(a.coord.top2eq_m(0, n.pi/2), 10)
        self.assertTrue(n.all(m == n.array([-e2,e1,e3])))

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.coord unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestConvert))
        self.addTests(loader.loadTestsFromTestCase(TestConvert_m))
        self.addTests(loader.loadTestsFromTestCase(TestRot_m))
        self.addTests(loader.loadTestsFromTestCase(TestXyz2thphi))
        self.addTests(loader.loadTestsFromTestCase(TestThphi2xyz))
        self.addTests(loader.loadTestsFromTestCase(TestEq2radec))
        self.addTests(loader.loadTestsFromTestCase(TestRadec2eq))
        self.addTests(loader.loadTestsFromTestCase(TestLatlong2xyz))
        self.addTests(loader.loadTestsFromTestCase(TestTop2azalt))
        self.addTests(loader.loadTestsFromTestCase(TestAzalt2top))
        self.addTests(loader.loadTestsFromTestCase(TestEq2top_m))
        self.addTests(loader.loadTestsFromTestCase(TestTop2eq_m))

if __name__ == '__main__':
    unittest.main()
