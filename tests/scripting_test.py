import unittest, aipy as a, numpy as n
from aipy.miriad import ij2bl

class TestParseAnts(unittest.TestCase):
    'Tests aipy.scripting.parse_ants()'
    def testant4(self):
        nants = 4
        cases = {
            'all': [],
            'auto': [('auto',1)],
            'cross': [('auto',0)],
            '0_1': [(ij2bl(0,1),1)],
            '0_1,1_2': [(ij2bl(0,1),1), (ij2bl(1,2),1)],
            '(0,1)_2': [(ij2bl(0,2),1), (ij2bl(1,2),1)],
            '0_(1,2)': [(ij2bl(0,1),1), (ij2bl(0,2),1)],
            '(0,1)_(2,3)': [(ij2bl(0,2),1), (ij2bl(0,3),1),
                            (ij2bl(1,2),1), (ij2bl(1,3),1)],
            '0_(1,-2)': [(ij2bl(0,1),1), (ij2bl(0,2),0)],
            '0,1,all': [],
        }
        for i in range(nants):
            cases[str(i)] = map(lambda x: (ij2bl(x,i),1), range(nants))
            cases['-'+str(i)] = map(lambda x: (ij2bl(x,i),0), range(nants))
        for ant_str in cases:
            self.assertEqual(a.scripting.parse_ants(ant_str, nants), 
                cases[ant_str])
        self.assertRaises(ValueError, a.scripting.parse_ants, '(0_1)_2', nants)

class TestParseChans(unittest.TestCase):
    'Tests aipy.scripting.parse_chans()'
    def setUp(self):
        self.nchan = 256
        self.cases = {
            'all': [n.arange(self.nchan)],
            '0,1,2,3': [n.array([0]),n.array([1]),n.array([2]),n.array([3])],
            '0_20': [n.arange(21)],
            '5,7,40_50': [n.array([5]),n.array([7]),n.arange(40,51)],
        }
    def testchan_concat(self):
        for case in self.cases:
            chans = a.scripting.parse_chans(case, self.nchan, concat=True)
            self.assertTrue(n.all(chans == n.concatenate(self.cases[case])))
        self.assertRaises(AssertionError,
            a.scripting.parse_chans, '0_1_2', self.nchan)
    def testchan_noconcat(self):
        for case in self.cases:
            chans = a.scripting.parse_chans(case, self.nchan, concat=False)
            for ch1, ch2 in zip(chans, self.cases[case]):
                self.assertTrue(n.all(ch1 == ch2))

if __name__ == '__main__':
    unittest.main()

