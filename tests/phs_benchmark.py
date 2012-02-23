# -*- coding: utf-8 -*-
import sys
import unittest
import timeit

class TestSpeed(unittest.TestCase):
    def test_phsspeed(self):
        """Test the speed of the aipy.phs module"""
        setup = '''
import numpy as n, aipy as a
freqs = n.arange(.1,.2,.0001)
beam = a.fit.Beam(freqs)
ants = [
    a.fit.Antenna(  0,   0,  0, beam),
    a.fit.Antenna(  0, 100,  0, beam),
    a.fit.Antenna(100,   0,  0, beam),
    a.fit.Antenna(100, 100,  0, beam),
]
aa = a.fit.AntennaArray(('45:00','90:00'), ants)
aa.set_jultime(2455400.1)
s_eqs = n.array([[0,1,0]]*100).transpose()
'''
        expr = '''
aa.gen_phs(s_eqs,0,1,)
'''
        t = timeit.Timer(expr, setup=setup)
        sys.stderr.write("%.3f ms ... " % ((t.timeit(number=10)/ 10) * 1e3))

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.phs unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestSpeed))

if __name__ == '__main__':
    unittest.main()
