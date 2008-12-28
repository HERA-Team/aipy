#! /usr/bin/env python

import aipy as a, numpy as n

crd1= n.array([0.,0.,1])
crd2= n.array([
    [1.,0.,0],
    [0.,1.,0],
    [0.,0.,1],
]).transpose()
freqs = n.array([.15, .16])

poly = n.array([
    [1.00, .00],
    [.10, .00],
])

b = a.fit.BeamCosSeries(freqs, poly)

print b.response(crd1)
print b.response(crd2)
