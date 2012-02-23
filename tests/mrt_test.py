import unittest, aipy._src.mrt as h, aipy as a, numpy as n

class TestMRTCatalog(unittest.TestCase):
    def setUp(self):
        self.cat = h.MRTCatalog()
        self.cat.fromfile(h.MRTFILE)
    def test_spotcheck(self):
        for srcname in self.cat:
            src = self.cat[srcname]
            self.assertEqual(src.index, 0)
            self.assertEqual(len(src.srcshape), 3)
            self.assertEqual(src.srcshape[0], 0)
            self.assertEqual(src.srcshape[1], 0)
            self.assertEqual(src.srcshape[2], 0)
            self.assertEqual(src.mfreq, .150)

if __name__ == '__main__':
    unittest.main()
