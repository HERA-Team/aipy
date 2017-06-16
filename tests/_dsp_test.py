import unittest
import aipy._dsp as _dsp
import numpy as n

class Testgrid1D_c(unittest.TestCase):
    def test_sanity(self):
        buf = n.zeros(32, dtype=n.complex64)
        ind = n.array([5, 10.1, 14.9], dtype=n.float32)
        dat = n.array([1, 1, 1], dtype=n.complex64)
        _dsp.grid1D_c(buf, ind, dat)
        x = n.arange(32)
        ans = n.exp(-(x-5)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2) + \
            n.exp(-(x-10.1)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2) + \
            n.exp(-(x-14.9)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2)
        if False:
            from matplotlib import pylab as P
            P.plot(ind, n.abs(dat), '.')
            P.plot(n.abs(buf))
            P.semilogy(x, n.exp(-(x-5)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            P.semilogy(x, n.exp(-(x-10.1)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            P.semilogy(x, n.exp(-(x-14.9)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            P.semilogy(x, n.abs(ans))
            P.ylim(1e-10, 1)
            P.show()
        self.assertAlmostEqual(n.max(n.abs(buf - ans)), 0, 6)

class Testgrid2D_c(unittest.TestCase):
    def test_sanity(self):
        buf = n.zeros((32,32), dtype=n.complex64)
        ind = n.array([[5,5], [10.1,10.1], [14.5, 15.5]], dtype=n.float32)
        ind1 = ind[:,0].copy()  # copy necessary b/c not contigous memory otherwise
        ind2 = ind[:,1].copy()
        dat = n.array([1, 1, 1], dtype=n.complex64)
        _dsp.grid2D_c(buf, ind1, ind2, dat)
        self.assertAlmostEqual(buf[5,5], 0.63661977236758149, 6)
        self.assertAlmostEqual(buf[10,10], 0.63661977236758149 * n.exp(-2*(.1**2+.1**2)), 6)
        self.assertAlmostEqual(buf[14,15], 0.63661977236758149 * n.exp(-2*(.5**2+.5**2)), 6)
        if False:
            from matplotlib import pylab as P
            P.imshow(n.log10(n.abs(buf)), vmax=0, vmin=-6, interpolation='nearest')
            #P.plot(n.abs(buf))
            #P.semilogy(x, n.exp(-(x-5)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            #P.semilogy(x, n.exp(-(x-10.1)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            #P.semilogy(x, n.exp(-(x-14.9)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            #P.semilogy(x, n.abs(ans))
            #P.ylim(1e-10, 1)
            P.show()

class Testdegrid2D_c(unittest.TestCase):
    def test_sanity(self):
        buf = n.ones((32,32), dtype=n.complex64)
        ind = n.array([[5,5], [10.1,10.1], [14.5, 15.5]], dtype=n.float32)
        ind1 = ind[:,0].copy()  # copy necessary b/c not contigous memory otherwise
        ind2 = ind[:,1].copy()
        dat = n.zeros(ind1.shape, dtype=n.complex64)
        _dsp.degrid2D_c(buf, ind1, ind2, dat)
        self.assertTrue(n.all(dat == 1))
        if False:
            from matplotlib import pylab as P
            P.imshow(n.log10(n.abs(buf)), vmax=0, vmin=-6, interpolation='nearest')
            #P.plot(n.abs(buf))
            #P.semilogy(x, n.exp(-(x-5)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            #P.semilogy(x, n.exp(-(x-10.1)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            #P.semilogy(x, n.exp(-(x-14.9)**2 / (2*.5**2))/n.sqrt(2*n.pi*.5**2))
            #P.semilogy(x, n.abs(ans))
            #P.ylim(1e-10, 1)
            P.show()

if __name__ == '__main__':
    unittest.main()
