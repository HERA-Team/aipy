# -*- coding: utf-8 -*-
import unittest, aipy as a, numpy as n, re
from aipy.miriad import ij2bl

class TestParseAnts(unittest.TestCase):
    'Tests aipy.scripting.parse_ants()'
    def testant4(self):
        """Test aipy.scripting.parse_ants()"""
        nants = 4
        cases = {
            'all': [],
            'auto': [('auto',1)],
            'cross': [('auto',0)],
            '0_1': [(ij2bl(0,1),1)],
            '0_1,1_2': [(ij2bl(0,1),1), (ij2bl(1,2),1)],
            '0x_1x': [(ij2bl(0,1),1,'xx')],
            '(0x,0y)_1x': [(ij2bl(0,1),1,'xx'), (ij2bl(0,1),1,'yx')],
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
        # inelegantly paste on the new pol parsing flag on the above tests
        # XXX really should add some new tests for the new pol parsing
        for k in cases:
            cases[k] = [(v+(-1,))[:3] for v in cases[k]]
        for ant_str in cases:
            self.assertEqual(a.scripting.parse_ants(ant_str, nants), 
                cases[ant_str])
        self.assertRaises(ValueError, a.scripting.parse_ants, '(0_1)_2', nants)

class TestParseSrcs(unittest.TestCase):
    'Tests aipy.scripting.parse_srcs()'
    def testsrc(self):
        """Test aipy.scripting.parse_srcs()"""
        cat = 'misc,helm'
        src_opt = 'cyg,cas,12_40,12:00_40:00,12:00:00.0_40:00:00.0,bug_man'
        srclist,cutoff,catalogs = a.scripting.parse_srcs(src_opt, cat)
        self.assertEqual(srclist[0], 'cyg')
        self.assertEqual(srclist[1], 'cas')
        self.assertEqual(srclist[2].src_name, '12_40')
        self.assertEqual(srclist[3].src_name, '12:00_40:00')
        self.assertEqual(srclist[4].src_name, '12:00:00.0_40:00:00.0')
        self.assertEqual(srclist[5], 'bug_man')

class TestParseChans(unittest.TestCase):
    'Tests aipy.scripting.parse_chans()'
    def setUp(self):
        self.nchan = 256
        self.cases = {
            'all': [n.arange(self.nchan)],
            '0,1,2,3': [n.array([0]),n.array([1]),n.array([2]),n.array([3])],
            '0_20': [n.arange(21)],
            '0_20_4': [n.arange(0,21,4)],
            '5,7,40_50': [n.array([5]),n.array([7]),n.arange(40,51)],
        }
    def testchan_concat(self):
        """Test aipy.scripting.parse_chans() - concatenate"""
        for case in self.cases:
            chans = a.scripting.parse_chans(case, self.nchan, concat=True)
            self.assertTrue(n.all(chans == n.concatenate(self.cases[case])))
        self.assertRaises(AssertionError,
            a.scripting.parse_chans, '0_1_2_3', self.nchan)
    def testchan_noconcat(self):
        """Test aipy.scripting.parse_chans() - without concatenate"""
        for case in self.cases:
            chans = a.scripting.parse_chans(case, self.nchan, concat=False)
            for ch1, ch2 in zip(chans, self.cases[case]):
                self.assertTrue(n.all(ch1 == ch2))

class TestParsePrms(unittest.TestCase):
    def setUp(self):
        self.test_strs = [
            'src1=jys',
            'src2=jys/100',
            'src3=jys/100/1',
            'src4=jys//1',
            '(src5/src6)=jys',
            '(src7/src8)=jys/(100/200)',
            '(src9/src10)=jys/(100/200)/(1/2)',
            '(src11/src12)=jys/(100/200)/1',
            '(src13/src14/src15)=(jys/index)//1',
        ]
    def testname(self):
        """Test aipy.scripting.name matching"""
        r = re.compile(a.scripting.name)
        self.assertEqual(r.match('src1').groups(), ('src1',))
        self.assertEqual(r.match('1').groups(), ('1',))
        self.assertEqual(r.match('1,1').groups(), ('1',))
        self.assertEqual(r.match('1=1').groups(), ('1',))
        self.assertEqual(r.match('1/1').groups(), ('1',))
        self.assertEqual(r.match('1(1').groups(), ('1',))
        self.assertEqual(r.match('1)1').groups(), ('1',))
        self.assertEqual(r.match(')1'), None)
    def testgrp(self):
        """Test aipy.scripting.grp matching"""
        r = re.compile(a.scripting.grp)
        self.assertEqual(r.match('src1').groups()[0], 'src1')
        self.assertEqual(r.match('src1/src2').groups()[0], 'src1')
        self.assertEqual(r.match('(src1/src2)').groups()[0], '(src1/src2)')
        self.assertEqual(r.match('(1/2/3)').groups()[0], '(1/2/3)')
        self.assertEqual(r.match('(src1/)'), None)
        self.assertEqual(r.match('(src1/src2'), None)
    def testprm(self):
        """Test aipy.scripting.prm_rgx matching"""
        r = a.scripting.prm_rgx
        for t in self.test_strs:
            self.assertEqual(r.match(t).groups()[0], t)
    def testprmlist(self):
        """Test aipy.scripting.parse_prms()"""
        t = ','.join(self.test_strs)
        prms = a.scripting.parse_prms(t)
        self.assertEqual(prms['src1']['jys'], (None,None))
        self.assertEqual(prms['src2']['jys'], (100.,None))
        self.assertEqual(prms['src3']['jys'], (100.,1.))
        self.assertEqual(prms['src4']['jys'], (None,1.))
        self.assertEqual(prms['src5']['jys'], (None,None))
        self.assertEqual(prms['src6']['jys'], (None,None))
        self.assertEqual(prms['src7']['jys'], (100.,None))
        self.assertEqual(prms['src8']['jys'], (200.,None))
        self.assertEqual(prms['src9']['jys'], (100.,1.))
        self.assertEqual(prms['src10']['jys'], (200.,2.))
        self.assertEqual(prms['src11']['jys'], (100.,1.))
        self.assertEqual(prms['src12']['jys'], (200.,1.))
        self.assertEqual(prms['src13']['jys'], (None,1.))
        self.assertEqual(prms['src13']['index'], (None,1.))
        self.assertEqual(prms['src14']['jys'], (None,1.))
        self.assertEqual(prms['src14']['index'], (None,1.))
        self.assertEqual(prms['src15']['jys'], (None,1.))
        self.assertEqual(prms['src15']['index'], (None,1.))
        self.assertRaises(AssertionError, 
            a.scripting.parse_prms,'(a/b)=(c/d)/(1/2)/(3/4)')
        t = '(1/2/3)=jys,(2/3)=index'
        prms = a.scripting.parse_prms(t)
        self.assertEqual(len(prms['1']), 1)
        self.assertEqual(len(prms['2']), 2)
        self.assertEqual(len(prms['3']), 2)
        prms = a.scripting.parse_prms('a=(b/c)')
        self.assertEqual(len(prms['a']), 2)

class TestSuite(unittest.TestSuite):
    """A unittest.TestSuite class which contains all of the aipy.coord unit tests."""

    def __init__(self):
        unittest.TestSuite.__init__(self)

        loader = unittest.TestLoader()
        self.addTests(loader.loadTestsFromTestCase(TestParseAnts))
        self.addTests(loader.loadTestsFromTestCase(TestParseSrcs))
        self.addTests(loader.loadTestsFromTestCase(TestParseChans))
        self.addTests(loader.loadTestsFromTestCase(TestParsePrms))

if __name__ == '__main__':
    unittest.main()

