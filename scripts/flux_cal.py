#! /usr/bin/env python
"""
A script for applying the flux calibration modeled in "aipy.loc" to the data.

Author: Aaron Parsons
Date: 6/03/07
Revisions: None
"""

import sys, aipy, numpy, os
from optparse import OptionParser

p = OptionParser()
p.set_usage('flux_cal.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
opts, args = p.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
aa = aipy.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

opts, args = p.parse_args(sys.argv[1:])

def cal_func(uv, p, d):
    i, j = aa.bl2ij(p[-1])
    bp = aa.ants[j].gain * numpy.conj(aa.ants[i].gain)
    bp = numpy.where(bp < numpy.average(bp) / 2, 1, bp)
    return p, d / bp

# Process all files passed from the command line.
for filename in args:
    print filename
    if os.path.exists(filename+'f'):
        print 'File exists: skipping'
        continue
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(filename+'f', status='new')
    aipy.miriad.pipe_uv(uvi, uvo, mfunc=cal_func)
    del(uvi); del(uvo)
