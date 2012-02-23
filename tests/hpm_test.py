#! /usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import aipy.healpix as h, numpy as n

hpm = h.HealpixMap(32, 'RING')

c1 = n.ones((20,), dtype=n.double)
c2 = n.ones((20,), dtype=n.double)
c3 = n.ones((20,), dtype=n.double)

print hpm[0], hpm[0,0], hpm[0,0,0]
hpm[0] = 1
hpm[0,0] = 1
hpm[0,0,0] = 1
print hpm[0], hpm[0,0], hpm[0,0,0]

c = n.arange(15)
print hpm[c]
hpm[c] += c
print hpm[c]
hpm[c] = c
print hpm[c]

c1 = n.arange(-.5,.5,.1)
c2 = n.arange(-.5,.5,.1)
c3 = n.arange(-.5,.5,.1)

print hpm[c1,c2,c3]
hpm[c1,c2,c3] += 1
print hpm[c1,c2,c3]

hpm.set_map(n.arange(48, dtype=n.double))
print hpm.get_map()
hpm.change_scheme('NEST')
print hpm.get_map()
