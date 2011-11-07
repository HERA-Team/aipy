import unittest, numpy as np, os
import aipy.miriad as m, aipy._miriad as _m

class TestMiriadUV(unittest.TestCase):
    def setUp(self):
        self.filename1 = '/tmp/test1.uv'
        self.filename2 = '/tmp/test2.uv'
        if os.path.exists(self.filename1): os.system('rm -rf %s' % self.filename1)
        if os.path.exists(self.filename2): os.system('rm -rf %s' % self.filename2)
        uv = m.UV(self.filename1, status='new', corrmode='r')
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
    def test_immediate_corr(self):
        uv = m.UV(self.filename2, 'new')
        self.assertEqual(uv.vartable['corr'], 'r')
    def test_vartable(self):
        uv = m.UV(self.filename1)
        self.assertEqual(uv.vartable['corr'], 'r')
        self.assertEqual(uv.vartable['nchan'], 'i')
        self.assertEqual(uv.vartable['pol'], 'i')
    def test_data(self):
        uv = m.UV(self.filename1)
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
        os.system('rm -rf %s' % self.filename1)
        os.system('rm -rf %s' % self.filename2)

if __name__ == '__main__':
    unittest.main()
