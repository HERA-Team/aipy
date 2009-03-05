#! /usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, os, sys, optparse

def gen_skypass_delay(aa, sdf, nchan, max_bl_frac=1.5):
    bin_dly = 1. / (sdf * nchan)
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j <= i: continue
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)
        max_bl = max_bl_frac * n.sqrt(n.dot(max_bl, max_bl))
        dly,off = aa.get_phs_offset(i,j)[-2:]
        uthresh, lthresh = (dly + max_bl)/bin_dly + 1, (dly - max_bl)/bin_dly
        uthresh, lthresh = int(n.round(uthresh)), int(n.round(lthresh))
        #f = n.ones((nchan,), dtype=n.float)
        #f[uthresh:lthresh] = 0
        filters[bl] = (uthresh,lthresh)
    return filters

o = optparse.OptionParser()
o.set_usage('filter_src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, loc=True)
o.add_option('-r', '--drw', dest='drw', type=int, default=-1,
    help='The number of fringe bins to null.  Default is -1 = no fringe filtering.')
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-e', '--extract', dest='extract', action='store_true',
    help='Extract the source instead of removing it.')
o.add_option('-t', '--together', dest='together', action='store_true',
    help='Perform delay-rate transform over all files, not file-by-file.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
if not opts.src is None:
    aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
    srclist,cutoff = a.scripting.parse_srcs(opts.src)
    src = a.loc.get_catalog(opts.loc, srclist, cutoff).values()[0]
else:
    aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
    src = None
    filters = gen_skypass_delay(aa, uv['sdf'], uv['nchan'])
del(uv)

phs_dat = {}
for uvfile in args:
    if not opts.together: phs_dat = {}
    if src is None: uvofile = uvfile + 'f'
    elif opts.extract: uvofile = uvfile+'.e_'+opts.src
    else: uvofile = uvfile+'.'+opts.src
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        if opts.together:
            print uvofile, 'exists, aborting.'
            sys.exit(0)
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)

    curtime = None
    print '    Performing delay transform and cleaning...'
    for (uvw,t,(i,j)),d,f in uvi.all(raw=True):
        if i == j: continue
        if curtime != t:
            curtime = t
            if not src is None:
                aa.set_jultime(t)
                src.compute(aa)
        bl = a.miriad.ij2bl(i,j)
        try:
            flags = n.logical_not(f).astype(n.float)
            gain = n.sqrt(n.average(flags**2))
            ker = n.fft.ifft(flags)
            d = n.where(f, 0, d)
            if not src is None: d = aa.phs2src(d, src, i, j)
            d = n.fft.ifft(d)
            if not n.all(d == 0):
                d, info = a.deconv.clean(d, ker, tol=opts.clean)
                d += info['res'] / gain
        except(a.ant.PointingError): d = n.zeros_like(d)
        try: phs_dat[bl].append(d)
        except(KeyError): phs_dat[bl] = [d]

    if not opts.together:
        if opts.drw != -1:
            print '    Performing delay-rate transform and cleaning...'
        for bl in phs_dat:
            d = n.array(phs_dat[bl])
            if opts.drw != -1:
                flags = n.where(d[:,0] != 0, 1., 0.)
                gain = n.sqrt(n.average(flags**2))
                ker = n.fft.ifft(flags)
                d = n.fft.ifft(d, axis=0)
                for chan in range(d.shape[1]):
                    d[:,chan],info = a.deconv.clean(d[:,chan],ker,tol=opts.clean)
                    d[:,chan] += info['res'] / gain
                x1, x2 = opts.drw, -opts.drw+1
                if x2 == 0: x2 = d.shape[0]
                d[x1:x2,:] = 0
            if opts.dw != -1: 
                y1, y2 = opts.dw, -opts.dw
                if y2 == 0: y2 = d.shape[1]
            else: y1, y2 = filters[bl]
            d[:,y1:y2] = 0
            if opts.drw != -1: d = n.fft.fft(d, axis=0)
            d = n.fft.fft(d, axis=1)
            phs_dat[bl] = d

        cnt = {}
        for bl in phs_dat: cnt[bl] = 0

        # Generate a pipe for removing average phase bias from data
        def rm_mfunc(uv, p, d):
            uvw, t, (i,j) = p
            if i == j: return p, d
            bl = a.miriad.ij2bl(i,j)
            if not src is None:
                aa.set_jultime(t)
                src.compute(aa)
            data = phs_dat[bl][cnt[bl],:]
            if not src is None:
                try: data = aa.unphs2src(data, src, i, j)
                except(a.ant.PointingError): data *= 0
            cnt[bl] += 1
            if opts.extract: return p, n.ma.array(data, mask=d.mask)
            else: return p, d - data

        print '    Writing out filtered data...'
        # Apply the pipe to the data
        uvi.rewind()
        uvo = a.miriad.UV(uvofile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=rm_mfunc, append2hist='FILTER_SRC: src=%s drw=%d dw=%d extract=%s clean=%f\n' % (opts.src, opts.drw, opts.dw, opts.extract, opts.clean))

if opts.together:
    print '    Performing delay-rate transform and cleaning...'
    for bl in phs_dat:
        d = n.array(phs_dat[bl])
        if opts.drw != -1:
            flags = n.where(d[:,0] != 0, 1., 0.)
            gain = n.sqrt(n.average(flags**2))
            ker = n.fft.ifft(flags)
            d = n.fft.ifft(d, axis=0)
            for chan in range(d.shape[1]):
                d[:,chan],info = a.deconv.clean(d[:,chan],ker,tol=opts.clean)
                d[:,chan] += info['res'] / gain
            x1, x2 = opts.drw, -opts.drw+1
            if x2 == 0: x2 = d.shape[0]
            d[x1:x2,:] = 0
        if opts.dw != -1: 
            y1, y2 = opts.dw, -opts.dw
            if y2 == 0: y2 = d.shape[1]
        else: y1, y2 = filters[bl]
        d[:,y1:y2] = 0
        if opts.drw != -1: d = n.fft.fft(d, axis=0)
        d = n.fft.fft(d, axis=1)
        phs_dat[bl] = d

    cnt = {}
    for bl in phs_dat: cnt[bl] = 0

    # Generate a pipe for removing average phase bias from data
    def rm_mfunc(uv, p, d):
        uvw, t, (i,j) = p
        if i == j: return p, d
        bl = a.miriad.ij2bl(i,j)
        if not src is None:
            aa.set_jultime(t)
            src.compute(aa)
        data = phs_dat[bl][cnt[bl],:]
        if not src is None:
            try: data = aa.unphs2src(data, src, i, j)
            except(a.ant.PointingError):
                if opts.extract: d *= 0
                data = 0
        cnt[bl] += 1
        if opts.extract: return p, n.ma.array(data, mask=d.mask)
        else: return p, d - data

    for uvfile in args:
        if src is None: uvofile = uvfile + 'f'
        elif opts.extract: uvofile = uvfile+'.e_'+opts.src
        else: uvofile = uvfile+'.'+opts.src
        print uvfile,'->',uvofile
        if not opts.together: phs_dat = {}
        uvi = a.miriad.UV(uvfile)
        uvo = a.miriad.UV(uvofile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=rm_mfunc, append2hist='FILTER_SRC: src=%s drw=%d dw=%d extract=%s clean=%f\n' % (opts.src, opts.drw, opts.dw, opts.extract, opts.clean))
