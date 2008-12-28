#! /usr/bin/env python
"""
A script for rotating UV data to a particular source.

Author: Aaron Parsons
Date: 6/03/07
Revisions:
    12/11/07 arp    Ported to use new miriad file interface
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
    uvw, t, (i,j) = p
    if i == j: return p, d
    aa.set_jultime(t)
    cat.compute(aa)
    try: d = aa.phs2src(d, cat.values()[0], i, j)
    except(aipy.ant.PointingError): d *= 0
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
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=phs)

