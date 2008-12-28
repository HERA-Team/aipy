#! /usr/bin/env python
"""
A script for rotating UV data to a particular source.

Author: Aaron Parsons
Date: 6/03/07
Revisions: None
"""

import aipy, numpy, sys, os
from optparse import OptionParser

p = OptionParser()
p.set_usage('phs2src.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-s', '--source', dest='source',
    help='The source to phase to.')
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
opts, args = p.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
aa = aipy.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'],
    use_bp=False)
del(uv)

cat = aipy.src.get_catalog([opts.source])

# A pipe to use for phasing to a source
def phs(uv, p, d):
    bl = p[-1]
    i, j = aipy.miriad.bl2ij(bl)
    if i == j: return p, d
    aa.set_jultime(p[-2])
    cat.compute(aa)
    try:
        d = aa.phs2src(d, cat.values()[0], bl)
    except:
        print 'Warning: source below horizon'
        d *= 0
    return p, d

# Process data
for filename in args:
    print filename
    uvofile = filename + '.' + opts.source
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(uvofile, status='new')
    aipy.miriad.pipe_uv(uvi, uvo, mfunc=phs)
    del(uvi); del(uvo)

