#! /usr/bin/env python
"""
Rotating UV data to a particular source.  Optionally adjust amplitudes to
remove model passband and beam effects.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('phs2src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, loc=True, src=True)
o.add_option('--passband', dest='passband', action='store_true',
    help='Divide out model passband.')
o.add_option('--beam', dest='beam', action='store_true',
    help='Divide out model beam response toward source.')
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
src = a.scripting.parse_srcs(opts.src)
del(uv)

# A pipe to use for phasing to a source
curtime,cache,top = None, {}, None
def phs(uv, p, d, f):
    global curtime, cache, top
    uvw, t, (i,j) = p
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        src.compute(aa)
        if opts.beam:
            top = a.coord.azalt2top((src.az,src.alt))
            cache = {}
    if opts.passband:
        # Remove shape of model passband
        passband = aa.ants[i].gain * n.conjugate(aa.ants[j].gain)
        d = n.where(passband == 0, 0, d/passband)
    if i == j: return p, d, f
    try:
        d = aa.phs2src(d, src, i, j)
        if opts.beam:
            p1,p2 = a.miriad.pol2str[uv['pol']]
            # Cache beam response, since it is an expensive operation
            if not cache.has_key(i): cache[i] = {}
            if not cache[i].has_key(p1):
                r = aa.ants[i].bm_response(top, pol=p1)
                cache[i][p1] = r.flatten()
            if not cache.has_key(j): cache[j] = {}
            if not cache[j].has_key(p2):
                r = aa.ants[j].bm_response(top, pol=p2)
                cache[j][p2] = r.flatten()
            # Calculate beam strength
            w = cache[i][p1] * cache[j][p2]
            d = n.where(w > 0, d/w, 0)
    except(a.ant.PointingError): d *= 0
    return p, d, f

# Process data
for filename in args:
    print filename
    uvofile = filename + '.' + opts.src
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=phs, raw=True)

