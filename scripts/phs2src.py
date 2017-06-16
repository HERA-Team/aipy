#! /usr/bin/env python
"""
Rotate zenith UV data to a particular source.  Can specify 'zen' to phase data
to zenith, or nothing at all to just remove delay/offset phase components.
"""

import aipy as a, numpy as np, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('phs2src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--setphs', dest='setphs', action='store_true',
    help='Instead of rotating phase, assign a phase corresponding to the specified source.')
o.add_option('--rot_uvw',action='store_true',
    help='Rotate the uvw coordinates to the source. Useful for exporting beyond AIPY')
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
if not opts.src is None:
    if not opts.src.startswith('zen'):
        srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
        src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]
    else: src = 'z'
else: src = None
del(uv)

# A pipe to use for phasing to a source
curtime = None
def phs(uv, p, d, f):
    global curtime
    uvw, t, (i,j) = p
    pol = a.miriad.pol2str[uv['pol']]
    aa.set_active_pol(pol)
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        if not src is None and not type(src) == str: src.compute(aa)
    if i == j: return p, d, f
    try:
        if opts.setphs: d = aa.unphs2src(np.abs(d), src, i, j)
        elif src is None: d *= np.exp(-1j*np.pi*aa.get_phs_offset(i,j))
        else: d = aa.phs2src(d, src, i, j)
        if opts.rot_uvw: 
            uvw = aa.gen_uvw(i,j,src=src)
            p = (uvw,t,(i,j))
    except(a.phs.PointingError): d *= 0
    return p, d, f

# Process data
for filename in args:
    if not opts.src is None: uvofile = filename + '.' + opts.src
    else: uvofile = filename + 'P'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=phs, raw=True)

