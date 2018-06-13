# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import unittest, numpy as n
import aipy.twodgauss as g

class TestTwodgauss(unittest.TestCase):
    def test_moments(self):
        """Test twodgauss.moments"""

        a0 = 2.0
        x0, y0 = 5.0, 7.0
        wx, wy = 1.0, 1.0
        x, y = n.arange(11, dtype=n.float64), n.arange(11, dtype=n.float64)
        x, y = n.meshgrid(x, y)
        data = a0*n.exp(-((x-x0)/wx)**2 / 2.0 - ((y-y0)/wy)**2 / 2.0)
        params = g.moments(data)
        self.assertAlmostEqual(a0, params[1], 3)
        self.assertAlmostEqual(x0, params[2], 3)
        self.assertAlmostEqual(y0, params[3], 3)
        #self.assertAlmostEqual(wx, params[4], 3)	# Why not?
        #self.assertAlmostEqual(wy, params[5], 3)	# Why not?

    def test_twodgaussian(self):
        """Test twodgauss.twodgaussian"""

        a0 = 2.0
        x0, y0 = 5.0, 7.0
        wx, wy = 1.0, 1.0
        x, y = n.arange(11, dtype=n.float64), n.arange(11, dtype=n.float64)
        x, y = n.meshgrid(x, y)
        data0 = a0*n.exp(-((x-x0)/wx)**2 / 2.0 - ((y-y0)/wy)**2 / 2.0)
        data1 = g.twodgaussian([0.0, a0, x0, y0, wx, wy], shape=y.shape)
        for i in range(data0.shape[0]):
            for j in range(data0.shape[1]):
                self.assertAlmostEqual(data0[i,j], data1[i,j], 6)

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.phs unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestTwodgauss))

if __name__ == '__main__':
    unittest.main()
