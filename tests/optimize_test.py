# -*- coding: utf-8 -*-

# Python3 compatibility
from __future__ import print_function, division, absolute_import

import unittest, numpy as n
import aipy.optimize as o

class TestOptimize(unittest.TestCase):
    def test_fmin(self):
        """Test optimize.fmin"""
        x0 = [0.8, 1.2, 0.7]
        fit = o.fmin(o.rosen, x0, disp=0)
        self.assertAlmostEqual(1.0, fit[0], 4)
        self.assertAlmostEqual(1.0, fit[1], 4)
        self.assertAlmostEqual(1.0, fit[2], 4)
        
    def test_fmin_powell(self):
        """Test optimize.fmin_powell"""
        x0 = [0.8, 1.2, 0.7]
        fit = o.fmin_powell(o.rosen, x0, disp=0)
        self.assertAlmostEqual(1.0, fit[0], 6)
        self.assertAlmostEqual(1.0, fit[1], 6)
        self.assertAlmostEqual(1.0, fit[2], 6)
        
    def test_fmin_ncg(self):
        """Test optimize.fmin_ncg"""
        x0 = [0.8, 1.2, 0.7]
        fit = o.fmin_ncg(o.rosen, x0, o.rosen_der, disp=0, maxiter=80)
        self.assertAlmostEqual(1.0, fit[0], 6)
        self.assertAlmostEqual(1.0, fit[1], 6)
        self.assertAlmostEqual(1.0, fit[2], 6)
        
    def test_fmin_ncg_hess_prod(self):
        """Test optimize.fmin_ncg with Hessian product"""
        x0 = [0.8, 1.2, 0.7]
        fit = o.fmin_ncg(o.rosen, x0, o.rosen_der, fhess_p=o.rosen_hess_prod, disp=0, maxiter=80)
        self.assertAlmostEqual(1.0, fit[0], 6)
        self.assertAlmostEqual(1.0, fit[1], 6)
        self.assertAlmostEqual(1.0, fit[2], 6)
        
    def test_fmin_ncg_hess_full(self):
        """Test optimize.fmin_ncg with full Hessian"""
        x0 = [0.8, 1.2, 0.7]
        fit = o.fmin_ncg(o.rosen, x0, o.rosen_der, fhess=o.rosen_hess, disp=0, maxiter=80)
        self.assertAlmostEqual(1.0, fit[0], 6)
        self.assertAlmostEqual(1.0, fit[1], 6)
        self.assertAlmostEqual(1.0, fit[2], 6)

class TestAnneal(unittest.TestCase):
    def test_cauchy(self):
        """Test optimize.anneal with the cauchy schedule"""
        func = lambda x: n.cos(14.5*x-0.3) + (x+0.2)*x
        fit = o.anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='cauchy')
        self.assertAlmostEqual(-0.195, fit[0], 3)
        
    def test_fast(self):
        """Test optimize.anneal with the fast schedule"""
        func = lambda x: n.cos(14.5*x-0.3) + (x+0.2)*x
        fit = o.anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='fast')
        self.assertAlmostEqual(-0.195, fit[0], 3)
        
    def test_boltzmann(self):
        """Test optimize.anneal with the boltzmann schedule"""
        func = lambda x: n.cos(14.5*x-0.3) + (x+0.2)*x
        fit = o.anneal(func,1.0,full_output=1,upper=3.0,lower=-3.0,feps=1e-4,maxiter=2000,schedule='boltzmann')
        self.assertAlmostEqual(-0.195, fit[0], 3)
        
    #def test_cauchy_twovar(self):
    #    func = lambda x: n.cos(14.5*x[0]-0.3) + (x[1]+0.2)*x[1] + (x[0]+0.2)*x[0]
    #    fit = o.anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='cauchy')
    #    self.assertAlmostEqual(-0.195, fit[0][0], 3)
    #    self.assertAlmostEqual(-0.100, fit[0][1], 3)
        
    #def test_fast_twovar(self):
    #    func = lambda x: n.cos(14.5*x[0]-0.3) + (x[1]+0.2)*x[1] + (x[0]+0.2)*x[0]
    #    fit = o.anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='fast')
    #    self.assertAlmostEqual(-0.195, fit[0][0], 3)
    #    self.assertAlmostEqual(-0.100, fit[0][1], 3)
        
    #def test_boltzmann_twovar(self):
    #    func = lambda x: n.cos(14.5*x[0]-0.3) + (x[1]+0.2)*x[1] + (x[0]+0.2)*x[0]
    #    fit = o.anneal(func,[1.0, 1.0],full_output=1,upper=[3.0, 3.0],lower=[-3.0, -3.0],feps=1e-4,maxiter=2000,schedule='boltzmann')
    #    self.assertAlmostEqual(-0.195, fit[0][0], 3)
    #    self.assertAlmostEqual(-0.100, fit[0][1], 3)

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.phs unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestOptimize))
        self.addTests(loader.loadTestsFromTestCase(TestAnneal))

if __name__ == '__main__':
    unittest.main()
