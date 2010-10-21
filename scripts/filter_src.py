#! /usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, os, sys, optparse, math

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
        filters[bl] = (uthresh,lthresh)
    return filters

o = optparse.OptionParser()
o.set_usage('filter_src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)
o.add_option('-r', '--drw', dest='drw', type=int, default=-1,
    help='The number of delay-rate bins to null.  Default is -1 = no fringe filtering.')
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-p','--passband', dest='passband', action='store_true',
    help='Divide by the model passband before transforming.')
# This doesn't work yet b/c it is complicated to implement.  Is it worth it?
#o.add_option('-b','--beam', dest='beam', action='store_true',
#    help='Divide by the model beam response before transforming.')
o.add_option('-e', '--extract', dest='extract', action='store_true',
    help='Extract the source instead of removing it.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
if not opts.src is None:
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
    src = cat[opts.src]
else:
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    src = None
    filters = gen_skypass_delay(aa, uv['sdf'], uv['nchan'])
del(uv)

phs_dat = {}
for uvfile in args:
    phs_dat = {}
    if src is None: uvofile = uvfile + 'f'
    elif opts.extract: uvofile = uvfile+'.e_'+opts.src
    else: uvofile = uvfile+'.'+opts.src
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)

    curtime = None
    print '    Performing delay transform and cleaning...'
    src_up = False
    for (uvw,t,(i,j)),d,f in uvi.all(raw=True):
        if i == j: continue
        if curtime != t:
            curtime = t
            if not src is None:
                aa.set_jultime(t)
                cat.compute(aa)
                #if not opts.beam:
                #    s_eqs = cat.get_crds('eq', ncrd=3)
                #    aa.sim_cache(s_eqs)
        bl = a.miriad.ij2bl(i,j)
        try:
            flags = n.logical_not(f).astype(n.float)
            # Check first if source is up
            if not src is None:
                d = aa.phs2src(d, src, i, j)
                # Put beam into kernel, where it gets divided out of the data
                #if opts.beam:
                #    pol = a.miriad.pol2str[uvi['pol']]
                #    flags *= aa.bm_response(i,j, pol).squeeze()
            # Put passband into kernel, where it gets divided out of the data
            if opts.passband: flags *= aa.passband(i,j)
            gain = n.sqrt(n.average(flags**2))
            ker = n.fft.ifft(flags)
            d = n.where(f, 0, d)
            src_up = True
            d = n.fft.ifft(d)
            if not n.all(d == 0):
                d, info = a.deconv.clean(d, ker, tol=opts.clean)
                d += info['res'] / gain
        except(a.phs.PointingError): d = n.zeros_like(d)
        try: phs_dat[bl].append(d)
        except(KeyError): phs_dat[bl] = [d]

    if opts.drw != -1:
        print '    Performing delay-rate transform and cleaning...'
    print '    Applying filter'
    for bl in phs_dat:
        d = n.array(phs_dat[bl])
        if not src_up:
            phs_dat[bl] = d
            continue
        # create some padding data on either end to mitigate wrap-around
        # effects of a delay-rate filter
        padlen = math.ceil(d.shape[0] * .1)
        d = n.concatenate([n.flipud(d[:padlen]), d, n.flipud(d[-padlen:])])
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
        # unpad the data
        d = d[padlen:-padlen]
        phs_dat[bl] = d

    cnt = {}
    for bl in phs_dat: cnt[bl] = 0

    # Generate a pipe for removing average phase bias from data
    curtime = None
    def rm_mfunc(uv, p, d):
        global curtime
        uvw, t, (i,j) = p
        if i == j: return p, d
        bl = a.miriad.ij2bl(i,j)
        data = phs_dat[bl][cnt[bl],:]
        if not src is None:
            if curtime != t:
                curtime = t
                aa.set_jultime(t)
                cat.compute(aa)
                #if opts.beam:
                #    s_eqs = cat.get_crds('eq', ncrd=3)
                #    aa.sim_cache(s_eqs)
            try:
                data = aa.unphs2src(data, src, i, j)
                # Remultiply by beam b/c the kernel divided it out
                #if opts.beam:
                #    pol = a.miriad.pol2str[uv['pol']]
                #    data *= aa.bm_response(i,j, pol).squeeze()
            except(a.phs.PointingError): data *= 0
        # Remultiply by passband b/c the kernel divided it out
        if opts.passband: data *= aa.passband(i,j)
        cnt[bl] += 1
        if opts.extract: return p, n.ma.array(data, mask=d.mask)
        else: return p, d - data

    print '    Writing out filtered data...'
    uvi.rewind()
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    # Apply the pipe to the data
    uvo.pipe(uvi, mfunc=rm_mfunc, append2hist='FILTER_SRC: src=%s drw=%d dw=%d extract=%s clean=%f passband=%s beam=%s\n' % (opts.src, opts.drw, opts.dw, opts.extract, opts.clean, opts.passband, False))#opts.beam))

