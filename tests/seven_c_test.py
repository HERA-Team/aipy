import unittest, aipy._src.seven_c as h, aipy as a, numpy as n

class TestSevenCCatalog(unittest.TestCase):
    def setUp(self):
        self.cat = h.SevenCCatalog()
        self.cat.fromfile(h.SEVENCFILE)
    def test_spotcheck(self):
        for srcname in self.cat:
            src = self.cat[srcname]
            self.assertEqual(src.index, 0)
            self.assertEqual(len(src.srcshape), 3)
            self.assertEqual(src.srcshape[0], 0)
            self.assertEqual(src.srcshape[1], 0)
            self.assertEqual(src.srcshape[2], 0)
            self.assertEqual(src.mfreq, .151)

if __name__ == '__main__':
    unittest.main()
