import unittest, aipy._src.three_cr as h, aipy as a, numpy as n

class TestThreeCRCatalog(unittest.TestCase):
    def setUp(self):
        self.cat = h.ThreeCRCatalog()
        self.cat.fromfile(h.THREECRFILE)
    def test_spotcheck(self):
        for srcname in self.cat:
            src = self.cat[srcname]
            self.assertEqual(src.index, 0)
            self.assertEqual(len(src.srcshape), 3)
            self.assertEqual(src.srcshape[0], 0)
            self.assertEqual(src.srcshape[1], 0)
            self.assertEqual(src.srcshape[2], 0)
            self.assertEqual(src.mfreq, .178)

if __name__ == '__main__':
    unittest.main()
