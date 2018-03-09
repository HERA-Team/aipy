#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function, division
import unittest
import aipy.healpix as h
import numpy as np

hpm = h.HealpixMap(32, 'RING')

c1 = np.ones((20,), dtype=np.double)
c2 = np.ones((20,), dtype=np.double)
c3 = np.ones((20,), dtype=np.double)

print(hpm[0], hpm[0,0], hpm[0,0,0])
hpm[0] = 1
hpm[0,0] = 1
hpm[0,0,0] = 1
print(hpm[0], hpm[0,0], hpm[0,0,0])

c = np.arange(15)
print(hpm[c])
hpm[c] += c
print(hpm[c])
hpm[c] = c
print(hpm[c])

c1 = np.arange(-.5,.5,.1)
c2 = np.arange(-.5,.5,.1)
c3 = np.arange(-.5,.5,.1)

print(hpm[c1,c2,c3])
hpm[c1,c2,c3] += 1
print(hpm[c1,c2,c3])

hpm.set_map(np.arange(48, dtype=np.double))
print(hpm.get_map())
hpm.change_scheme('NEST')
print(hpm.get_map())
