#! /usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as np, os, sys, optparse, math

def gen_skypass_delay(aa, sdf, nchan, pol, max_bl_frac=1.5):
    aa.set_active_pol(pol)
    bin_dly = 1. / (sdf * nchan)
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j <= i: continue
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)
        max_bl = max_bl_frac * np.sqrt(np.dot(max_bl, max_bl))
        dly,off = aa.get_phs_offset(i,j)[-2:]
        uthresh, lthresh = (dly + max_bl)/bin_dly + 1, (dly - max_bl)/bin_dly
        uthresh, lthresh = int(np.round(uthresh)), int(np.round(lthresh))
        filters[bl] = (uthresh,lthresh)
    return filters

o = optparse.OptionParser()
o.set_usage('filter_src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)
o.add_option('-p','--pol',dest='pol',default='xx,xy,yx,yy',help='Polarizations to use.')
o.add_option('-r', '--drw', dest='drw', type=int, default=-1,
    help='The number of delay-rate bins to null.  Default is -1 = no fringe filtering.')
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-P','--passband', dest='passband', action='store_true',
    help='Divide by the model passband before transforming.')
o.add_option('-e', '--extract', dest='extract', action='store_true',
    help='Extract the source instead of removing it.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

pols = opts.pol.split(',')

uv = a.miriad.UV(args[0])
if not opts.src is None:
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
    src = cat[opts.src]
else:
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    src = None
    filters = {}
    for pol in pols:
        filters[pol] = gen_skypass_delay(aa, uv['sdf'], uv['nchan'],pol)
del(uv)

for uvfile in args:
    phs_dat = {}
    if src is None: uvofile = uvfile + 'f'
    elif opts.extract: uvofile = uvfile+'.e_'+opts.src
    else: uvofile = uvfile+'.f_'+opts.src
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)

    curtime = None
    print '    Performing delay transform and cleaning...'
    src_up = False
    for (uvw,t,(i,j)),d,f in uvi.all(raw=True):
        pol = a.miriad.pol2str[uvi['pol']]
        if not pol in pols: continue
        aa.set_active_pol(pol)
        if i == j: continue
        if not pol in phs_dat.keys(): phs_dat[pol] = {}
        if curtime != t:
            curtime = t
            if not src is None:
                aa.set_jultime(t)
                cat.compute(aa)
        bl = a.miriad.ij2bl(i,j)
        try:
            flags = np.logical_not(f).astype(np.float)
            # Check first if source is up
            if not src is None:
                d = aa.phs2src(d, src, i, j)
            # Put passband into kernel, where it gets divided out of the data
            if opts.passband: flags *= aa.passband(i,j)
            gain = np.sqrt(np.average(flags**2))
            ker = np.fft.ifft(flags)
            d = np.where(f, 0, d)
            src_up = True
            d = np.fft.ifft(d)
            if not np.all(d == 0):
                d, info = a.deconv.clean(d, ker, tol=opts.clean)
                d += info['res'] / gain
        except(a.phs.PointingError): d = np.zeros_like(d)
        try: phs_dat[pol][bl].append(d)
        except(KeyError): phs_dat[pol][bl] = [d]

    if opts.drw != -1:
        print '    Performing delay-rate transform and cleaning...'
    for pol in phs_dat:
        for bl in phs_dat[pol]:
            d = np.array(phs_dat[pol][bl])
            if not src_up:
                phs_dat[pol][bl] = d
                continue
            # create some padding data on either end to mitigate wrap-around
            # effects of a delay-rate filter
            padlen = math.ceil(d.shape[0] * .1)
            d = np.concatenate([np.flipud(d[:padlen]), d, np.flipud(d[-padlen:])])
            if opts.drw != -1:
                flags = np.where(d[:,0] != 0, 1., 0.)
                gain = np.sqrt(np.average(flags**2))
                ker = np.fft.ifft(flags)
                d = np.fft.ifft(d, axis=0)
                for chan in range(d.shape[1]):
                    d[:,chan],info = a.deconv.clean(d[:,chan],ker,tol=opts.clean)
                    d[:,chan] += info['res'] / gain
                x1, x2 = opts.drw, -opts.drw+1
                if x2 == 0: x2 = d.shape[0]
                d[x1:x2,:] = 0
            if opts.dw != -1: 
                y1, y2 = opts.dw, -opts.dw
                if y2 == 0: y2 = d.shape[1]
            else: y1, y2 = filters[pol][bl]
            d[:,y1:y2] = 0
            if opts.drw != -1: d = np.fft.fft(d, axis=0)
            d = np.fft.fft(d, axis=1)
            # unpad the data
            d = d[padlen:-padlen]
            phs_dat[pol][bl] = d

    cnt = {}
    for pol in phs_dat:
        cnt[pol] = {}
        for bl in phs_dat[pol]: cnt[pol][bl] = 0

    # Generate a pipe for removing average phase bias from data
    curtime = None
    def rm_mfunc(uv, p, d):
        global curtime
        uvw, t, (i,j) = p
        if i == j: return p, d
        pol = a.miriad.pol2str[uv['pol']]
        aa.set_active_pol(pol)
        bl = a.miriad.ij2bl(i,j)
        data = phs_dat[pol][bl][cnt[pol][bl],:]
        if not src is None:
            if curtime != t:
                curtime = t
                aa.set_jultime(t)
                cat.compute(aa)
            try: data = aa.unphs2src(data, src, i, j)
            except(a.phs.PointingError): data *= 0
        # Remultiply by passband b/c the kernel divided it out
        if opts.passband:
            data *= aa.passband(i,j)
        cnt[pol][bl] += 1
        if opts.extract: return p, np.ma.array(data, mask=d.mask)
        else: return p, d - data

    print '    Writing out filtered data...'
    uvi.rewind()
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    # Apply the pipe to the data
    uvo.pipe(uvi, mfunc=rm_mfunc, append2hist='FILTER_SRC: src=%s drw=%d dw=%d extract=%s clean=%f passband=%s beam=%s\n' % (opts.src, opts.drw, opts.dw, opts.extract, opts.clean, opts.passband, False))
