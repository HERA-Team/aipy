# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import unittest
import aipy.healpix as a
import numpy as np

class TestAlm(unittest.TestCase):
    def test_init(self):
        """Test initializing a aipy._alm object"""
        alm = a.Alm(10, 5)
        self.assertEqual(alm.lmax(), 10)
        self.assertEqual(alm.mmax(), 5)
        self.assertRaises(Exception, a.Alm, 5, 10)
    def test_get_data(self):
        alm = a.Alm(1,1)
        data = alm.get_data()
        self.assertEqual(data.size, 3)
    def test_set_data(self):
        alm = a.Alm(1,1)
        data = np.array([1+1j, 2+2j, 3+3j], dtype=np.complex128)
        alm.set_data(data)
        np.testing.assert_equal(alm.get_data(), data)
    def test_set_to_zero(self):
        alm = a.Alm(3,3)
        np.testing.assert_equal(alm.get_data(), 0)
    def test_get_set(self):
        alm = a.Alm(1,1)
        data = np.array([1+1j, 2+2j, 3+3j], dtype=np.complex128)
        alm.set_data(data)
        self.assertEqual(alm[0,0], data[0])
        self.assertEqual(alm[1,0], data[1])
        self.assertEqual(alm[1,1], data[2])
        alm[1,1] = 4+4j
        self.assertEqual(alm[1,1], 4+4j)
    def test_to_map(self):
        alm = a.Alm(1,1)
        alm[0,0] = 1
        d = alm.to_map(32, 'RING')
        np.testing.assert_almost_equal(d, np.ones(12*32**2, dtype=np.float) * 0.2820948)
    def test_from_map(self):
        alm = a.Alm(1,1)
        d = np.ones(12*32**2, dtype=np.float32)
        alm.from_map(d, 3)
        self.assertAlmostEqual(alm[0,0], 3.5449077018110313)
        self.assertAlmostEqual(alm[1,0], 0)
        self.assertAlmostEqual(alm[1,1], 0)

if False:
  class TestMemLeaks(unittest.TestCase):
    def setUp(self):
        self.alm = a.Alm(10,10)
    def test_alm_create(self):
        while True: alm = a.Alm(20,20)
    def test_alm_to_map(self):
        while True: d = self.alm.to_map(256, 'RING')
    def test_alm_from_map(self):
        d = np.zeros(12*256**2, dtype=n.float)
        while True: self.alm.from_map(d, 2)
    def test_alm_get_data(self):
        while True: c = self.alm.get_data()
    def test_alm_set_data(self):
        c = self.alm.get_data()
        while True: self.alm.set_data(c)

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy._alm unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestAlm))
        #self.addTests(loader.loadTestsFromTestCase(TestMemLeaks))

if __name__ == '__main__':
    unittest.main()
