import unittest, aipy._src.culgoora as h, aipy as a, numpy as n

class TestCulgooraCatalog(unittest.TestCase):
    def setUp(self):
        self.cat = h.CulgooraCatalog()
        self.cat.fromfile(h.CULGOORAFILE)
    def test_spotcheck(self):
        for srcname in self.cat:
            src = self.cat[srcname]
            self.assertEqual(len(src.srcshape), 3)
            self.assertEqual(src.srcshape[0], 0)
            self.assertEqual(src.srcshape[1], 0)
            self.assertEqual(src.srcshape[2], 0)
            self.assertTrue(src.mfreq in [.080, .160])
    def test_metadata(self):
        md = self.cat.get_metadata()
        for src in md:
            self.assertEqual(len(md[src]), 3)
            jy080, jy160, index = md[src]
            if index is None:
                self.assertTrue(jy080 is None or jy160 is None)
            else:
                self.assertNotEqual(jy080, None)
                self.assertNotEqual(jy160, None)
                #self.assertAlmostEqual(index, n.log2(jy160/jy080), 1)
                

if __name__ == '__main__':
    unittest.main()
