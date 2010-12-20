import unittest, numpy as np, os
import aipy.miriad as m, aipy._miriad as _m

class TestMiriadUV(unittest.TestCase):
    def setUp(self):
        self.filename = '/tmp/test.uv'
        if os.path.exists(self.filename): os.system('rm -rf %s' % self.filename)
        uv = m.UV(self.filename, status='new', corrmode='r')
        uv['history'] = 'Made this file from scratch.\n'
        uv.add_var('nchan', 'i')
        uv.add_var('pol', 'i')
        uv['nchan'] = 4
        uv['pol'] = -5
        uvw = np.array([1,2,3], dtype=np.double)
        preamble = (uvw, 12345.6789, (0,1))
        self.data = np.ma.array([1j,2,3j,4], mask=[0,0,1,0], dtype=np.complex64)
        uv.write(preamble,self.data)
        uv['pol'] = -6
        uv.write(preamble,self.data)
    def test_vartable(self):
        uv = m.UV(self.filename)
        self.assertEqual(uv.vartable['corr'], 'r')
        self.assertEqual(uv.vartable['nchan'], 'i')
        self.assertEqual(uv.vartable['pol'], 'i')
    def test_data(self):
        uv = m.UV(self.filename)
        self.assertEqual(uv['history'], 'Made this file from scratch.\n')
        (uvw,t,bl),d = uv.read()
        self.assertEqual(uv['nchan'], 4)
        self.assertEqual(uv['pol'], -5)
        self.assertEqual(bl, (0,1))
        self.assertEqual(t, 12345.6789)
        self.assertTrue(np.all(uvw == np.array([1,2,3], dtype=np.double)))
        self.assertTrue(np.all(d == self.data))
        (uvw,t,bl),d = uv.read()
        self.assertEqual(uv['nchan'], 4)
        self.assertEqual(uv['pol'], -6)
        self.assertEqual(bl, (0,1))
        self.assertEqual(t, 12345.6789)
        self.assertTrue(np.all(uvw == np.array([1,2,3], dtype=np.double)))
        self.assertTrue(np.all(d == self.data))
    def tearDown(self):
        os.system('rm -rf %s' % self.filename)

if __name__ == '__main__':
    unittest.main()
