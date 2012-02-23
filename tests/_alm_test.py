import unittest, aipy._alm as a, numpy as n

class TestAlm(unittest.TestCase):
    def test_init(self):
        alm = a.Alm(10, 5)
        self.assertEqual(alm.lmax(), 10)
        self.assertEqual(alm.mmax(), 5)
        self.assertRaises(RuntimeError, a.Alm, 5, 10)

if False:
  class TestMemLeaks(unittest.TestCase):
    def setUp(self):
        self.alm = a.Alm(10,10)
    def test_alm_create(self):
        while True: alm = a.Alm(20,20)
    def test_alm_to_map(self):
        while True: d = self.alm.to_map(256, 'RING')
    def test_alm_from_map(self):
        d = n.zeros(12*256**2, dtype=n.float)
        while True: self.alm.from_map(d, 2)
    def test_alm_get_data(self):
        while True: c = self.alm.get_data()
    def test_alm_set_data(self):
        c = self.alm.get_data()
        while True: self.alm.set_data(c)

if __name__ == '__main__':
    unittest.main()
