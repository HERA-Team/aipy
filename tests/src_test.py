import unittest, aipy.src as src

class TestSrc(unittest.TestCase):
    def test_get_catalog(self):
        cat = src.get_catalog(srcs=['J0535+220'])
        self.assertTrue(cat.has_key('J0535+220'))
        s = cat['J0535+220']
        self.assertEqual(s._jys, 1883.)
        self.assertEqual(s.index, -99.)
        cat = src.get_catalog(srcs=['J0535+220'], cats=[])
        self.assertEqual(len(cat.keys()), 0)

if __name__ == '__main__':
    unittest.main()        
