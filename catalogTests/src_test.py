# -*- coding: utf-8 -*-
import unittest, aipy.src as src

class TestSrc(unittest.TestCase):
    def test_get_catalog(self):
        """Test setting up a catalog"""
        cat = src.get_catalog(srcs=['J0535+220'])
        self.assertTrue(cat.has_key('J0535+220'))
        s = cat['J0535+220']
        self.assertEqual(s._jys, 1883.)
        self.assertEqual(s.index, -99.)
        cat = src.get_catalog(srcs=['J0535+220'], catalogs=[])
        self.assertEqual(len(cat.keys()), 0)

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.src unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTest(loader.loadTestsFromTestCase(TestSrc))

if __name__ == '__main__':
    unittest.main()        
