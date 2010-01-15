import unittest, aipy._src.vlss as h, aipy as a, numpy as n

class TestVLSSCatalog(unittest.TestCase):
    def setUp(self):
        self.cat = h.VLSSCatalog()
        self.cat.fromfile(h.VLSSFILE)
    def test_spotcheck(self):
        for srcname in self.cat:
            src = self.cat[srcname]
            self.assertEqual(src.index, 0)
            self.assertEqual(len(src.srcshape), 3)
            self.assertEqual(src.srcshape[0], 0)
            self.assertEqual(src.srcshape[1], 0)
            self.assertEqual(src.srcshape[2], 0)
            self.assertEqual(src.mfreq, 0.074)

if __name__ == '__main__':
    unittest.main()
