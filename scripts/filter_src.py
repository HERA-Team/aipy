#! /usr/bin/env python
"""
A script for removing a source using a fringe-rate/delay transform.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, os, sys, optparse

o = optparse.OptionParser()
o.set_usage('filter_src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, loc=True)
o.add_option('-f', '--fng_w', dest='fng_w', type=int, default=1,
    help='The number of fringe bins to null.')
o.add_option('-d', '--dly_w', dest='dly_w', type=float, default=5,
    help='The number of delay bins to null.')
o.add_option('-e', '--extract', dest='extract', action='store_true',
    help='Extract the source instead of removing it.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
src = a.scripting.parse_srcs(opts.src)
del(uv)

for uvfile in args:
    print 'Working on', uvfile
    phs_dat = {}
    if opts.extract: uvofile = uvfile+'.e_'+opts.src
    else: uvofile = uvfile+'.'+opts.src
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)

    curtime = None
    print '    Performing delay transform and cleaning...'
    for (uvw,t,(i,j)),d,f in uvi.all(raw=True):
        if i == j: continue
        if curtime != t:
            curtime = t
            aa.set_jultime(t)
            src.compute(aa)
        bl = aa.ij2bl(i,j)
        try:
            flags = n.logical_not(f).astype(n.float)
            gain = n.sqrt(n.average(flags**2))
            ker = n.fft.ifft(flags)
            d = n.where(f, 0, d)
            d = aa.phs2src(d, src, i, j)
            d = n.fft.ifft(d)
            if not n.all(d == 0):
                d, info = a.deconv.clean1d(d, ker, tol=opts.clean)
                d += info['res'] / gain
        except(a.ant.PointingError): d = n.zeros_like(d)
        try: phs_dat[bl].append(d)
        except(KeyError): phs_dat[bl] = [d]

    print '    Performing fringe rate transform and cleaning...'
    for bl in phs_dat:
        d = n.array(phs_dat[bl])
        flags = n.where(d[:,0] != 0, 1., 0.)
        gain = n.sqrt(n.average(flags**2))
        ker = n.fft.ifft(flags)
        d = n.fft.ifft(d, axis=0)
        for chan in range(d.shape[1]):
            d[:,chan], info = a.deconv.clean1d(d[:,chan], ker, tol=opts.clean)
            d[:,chan] += info['res'] / gain
        x1, x2 = opts.fng_w, -opts.fng_w+1
        if x2 == 0: x2 = d.shape[0]
        y1, y2 = opts.dly_w, -opts.dly_w
        if y2 == 0: y2 = d.shape[1]
        d[x1:x2,:] = 0
        d[:,y1:y2] = 0
        d = n.fft.fft(d, axis=0)
        d = n.fft.fft(d, axis=1)
        phs_dat[bl] = d

    cnt = {}
    for bl in phs_dat: cnt[bl] = 0

    # Generate a pipe for removing average phase bias from data
    def rm_mfunc(uv, p, d):
        uvw, t, (i,j) = p
        if i == j: return p, d
        bl = aa.ij2bl(i,j)
        aa.set_jultime(t)
        src.compute(aa)
        data = phs_dat[bl][cnt[bl],:]
        try: data = aa.unphs2src(data, src, i, j)
        except(a.ant.PointingError):
            if opts.extract: d *= 0
            data = 0
        cnt[bl] += 1
        if opts.extract: return p, n.ma.array(data, mask=d.mask)
        else: return p, d - data

    print '    Writing out filtered data...'
    # Apply the pipe to the data
    uvi.rewind()
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=rm_mfunc, append2hist='FILTER_SRC: src=%s fng_w=%d dly_w=%d extract=%s clean=%f\n' % (opts.src, opts.fng_w, opts.dly_w, opts.extract, opts.clean))
    del(uvi); del(uvo)
