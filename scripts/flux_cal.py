#! /usr/bin/env python
"""
A script for dividing out the passband, primary beam, and/or source spectrum 
scaling.  When dividing by a primary beam or source spectrum, it is recommended
a single source have been isolated in the data set.

Author: Aaron Parsons
"""

import aipy as a, numpy as np, os, sys, optparse, pickle
from matplotlib import pylab as p

o = optparse.OptionParser()
o.set_usage('flux_cal.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)
o.add_option('-b', '--beam', dest='beam', action='store_true',
    help='Normalize by the primary beam response in the direction of the specified source.')
o.add_option('-p', '--passband', dest='passband', action='store_true',
    help='Normalize by the passband response.')
o.add_option('-f', '--srcflux', dest='srcflux', action='store_true',
    help='Normalize by the spectrum of the specified source.')
opts, args = o.parse_args(sys.argv[1:])

assert(not (opts.src is None and (opts.beam or opts.srcflux)))
uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

if opts.srcflux:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
    s = cat.values()[0]
    print 'Calibrating for source with',
    print 'strength', s._jys,
    print 'measured at', s.mfreq, 'GHz',
    print 'with index', s.index
    src_spec = None
else: src_spec = 1

curtime = None
def mfunc(uv, p, d, f):
    global curtime, src_spec
    uvw,t,(i,j) = p
    pol = a.miriad.pol2str[uv['pol']]
    aa.set_active_pol(pol)
    if t != curtime:
        curtime = t
        aa.set_jultime(t)
        if opts.srcflux:
            s.compute(aa)
            src_spec = s.get_jys()
            s_eq = cat.get_crds('eq', ncrd=3)
            aa.sim_cache(s_eq)
    if opts.passband:
        passband = aa.passband(i,j)
    else: passband = 1
    if opts.beam: bm_resp = aa.bm_response(i,j).squeeze()
    else: bm_resp = 1
    gain = passband * bm_resp * src_spec
    d /= np.where(gain == 0, 1, gain)
    return p, d, f

ext = '.'
if opts.passband: ext += 'p'
if opts.beam: ext += 'b'
if opts.srcflux: ext += 'f'
for filename in args:
    uvofile = filename + ext
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, 
        append2hist='FLUXCAL: srcs=%s passband=%s beam=%s srcflux=%s\n' % \
        (opts.src, opts.passband, opts.beam, opts.srcflux))
