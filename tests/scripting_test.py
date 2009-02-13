import unittest
import aipy as a
from aipy.miriad import ij2bl

class TestParseAnts(unittest.TestCase):
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

if __name__ == '__main__':
    unittest.main()

