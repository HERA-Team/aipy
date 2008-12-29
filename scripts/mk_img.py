#! /usr/bin/env python
"""
This is a general-purpose script for making images from MIRIAD UV files.  Data
(optionally selected for baseline, channel) are read from the file, phased
to a provided position, normalized for passband/primary beam effects, gridded
to a UV matrix, imaged, and optionally deconvolved by a corresponding PSF to
produce a clean image.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, optparse, ephem, os

o = optparse.OptionParser()
o.set_usage('mk_img.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, loc2=True,
    src=True, dec=True)
o.add_option('-r', '--residual', dest='residual', action='store_true',
    help='Display only the residual in the 4th panel (otherwise display sum of clean image and residual).')
o.add_option('-d', '--deconv', dest='deconv', default='mem',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).')
o.add_option('-o', '--output', dest='output', default='bim',
    help='Comma delimited list of data to generate FITS files for.  Can be: dim (dirty image), dbm (dirty beam), cim (clean image), rim (residual image), bim (best image = clean + residuals), uvs (uv sampling), or bms (beam sampling).  Default is bim.')
o.add_option('--npts', dest='npts', type='int', default=200,
    help='If no src is provided, facet the sphere into this many pointings for making a map.  Default 200.')
o.add_option('--prefix', dest='prefix', default='im_',
    help='A string to prefix to all filenames.')
o.add_option('--fmt', dest='fmt', default='%04d',
    help='A format string for counting successive images written to files.  Default is %04d (i.e. 0001).')
o.add_option('--var', dest='var', type='float', default=.6,
    help='Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.  For annealing, interpreted as cooling speed.')
o.add_option('--maxiter', dest='maxiter', type='int', default=200,
    help='Number of allowable iterations per deconvolve attempt.')
o.add_option('-u', '--uniform', dest='uniform', type='float', default=0,
    help="Use uniform (rather than natural) weighting for uv bins that have a weight above the specified fraction of the maximum weighting.")
o.add_option('--skip_amp', dest='skip_amp', action='store_true',
    help='Do not use amplitude information to normalize visibilities.')
o.add_option('--skip_bm', dest='skip_bm', action='store_true',
    help='Do not weight visibilities by the strength of the primary beam.')
o.add_option('--size', dest='size', type='int', default=200,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--buf_thresh', dest='buf_thresh', default=1.8e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
locs = a.scripting.files_to_locs(opts.loc, args, sys.argv)
aas = {}
for L in locs:
    uv = a.miriad.UV(locs[L][0])
    (j,t,j),j = uv.read()
    chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    aa = a.loc.get_aa(L, uv['sdf'], uv['sfreq'], uv['nchan'])
    aa.select_chans(chans)
    aas[L] = aa
    afreqs = aa.ants[0].beam.afreqs
    cfreq = n.average(afreqs)
    aa.set_jultime(t)
outputs = opts.output.split(',')

if not opts.src is None: cat = a.scripting.parse_srcs(opts.src, force_cat=True)
else:
    ras,decs = a.map.facet_centers(opts.npts, ncrd=2)
    srcs = [a.ant.RadioFixedBody(ra,dec,name=str(i)) 
        for i,(ra,dec) in enumerate(zip(ras,decs))]
    cat = a.ant.SrcCatalog(srcs)

if opts.no_w: im = a.img.Img(opts.size, opts.res, mf_order=0)
else: im = a.img.ImgW(opts.size, opts.res, mf_order=0)
DIM = int(opts.size/opts.res)
del(uv)

for imgcnt, s in enumerate(cat.values()):
    s.compute(aa)
    print '%d / %d' % (imgcnt + 1, len(cat.values()))
    print 'Pointing (ra, dec):', s._ra, s._dec
    if abs(aa.lat - s._dec) > n.pi/2:
        print '    Source never rises: skipping'
        continue
    src = a.fit.SrcCatalog([s])

    def to_fits(t,i):
        filename = (opts.prefix+t+opts.fmt+'.fits') % imgcnt
        print 'Saving data to', filename
        while len(i.shape) < 4: i.shape = i.shape + (1,)
        cen = ephem.Equatorial(s.ra, s.dec, epoch=aa.epoch)
        cen = ephem.Equatorial(cen, epoch=ephem.J2000)
        L,M = im.get_LM()
        a.img.to_fits(filename, i, clobber=True,
            object=s.src_name, obs_date=str(aa.date),
            ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
            d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[-1,-1]*a.img.rad2deg,
            freq=n.average(aa.ants[0].beam.afreqs),
        )

    fdim = (opts.prefix+'dim'+opts.fmt+'.fits') % imgcnt
    fdbm = (opts.prefix+'dbm'+opts.fmt+'.fits') % imgcnt
    if os.path.exists(fdim) and os.path.exists(fdbm):
        print 'Found %s and %s, so skipping gridding/phasing.' % (fdim,fdbm)
        dim, kwds = a.img.from_fits(fdim)
        dbm, kwds = a.img.from_fits(fdbm)
        dim,dbm = dim.squeeze(), dbm.squeeze()
    else:
        # Gather data
        us,vs,ws = [],[],[]
        ds,wgts = [], []
        cnt, curtime = 0, None
        if opts.no_w: im = a.img.Img(opts.size, opts.res, mf_order=0)
        else: im = a.img.ImgW(opts.size, opts.res, mf_order=0)
        for L in locs:
          aa = aas[L]
          # Read each file
          for filename in locs[L]:
            sys.stdout.write('.'); sys.stdout.flush()
            uv = a.miriad.UV(filename)
            a.scripting.uv_selector(uv, opts.ant, opts.pol)
            # Read all data from each file
            for (crd,t,(i,j)),d,f in uv.all(raw=True):
                if curtime != t:
                    curtime = t
                    cnt = (cnt + 1) % opts.decimate
                    if cnt == 0:
                        aa.set_jultime(t)
                        src.compute(aa)
                        s_eq = src.get_crds('eq', ncrd=3)
                        aa.sim_cache(s_eq)
                if cnt != 0: continue
                d,f = d.take(chans), f.take(chans)
                if not opts.skip_amp: d /= aa.passband(i,j)
                try:
                    # Throws PointingError if not up:
                    d = aa.phs2src(d, s, i, j)
                    u,v,w = aa.gen_uvw(i,j,src=s)
                    if not opts.skip_bm:
                        # Calculate beam strength for weighting purposes
                        wgt = aa.bm_response(i,j,pol=opts.pol).squeeze()
                        # Optimal SNR: down-weight beam-attenuated data 
                        # by another factor of the beam response.
                        d *= wgt; wgt *= wgt
                    else: wgt = n.ones(d.shape, dtype=n.float)
                except(a.ant.PointingError): continue
                valid = n.logical_not(f)
                d = d.compress(valid)
                if len(d) == 0: continue
                ds.append(d)
                us.append(u.compress(valid))
                vs.append(v.compress(valid))
                ws.append(w.compress(valid))
                wgts.append(wgt.compress(valid))
                # If data buffer is full, grid data
                if len(ds) * len(chans) > opts.buf_thresh:
                    sys.stdout.write('|'); sys.stdout.flush()
                    ds = n.concatenate(ds)
                    us = n.concatenate(us)
                    vs = n.concatenate(vs)
                    ws = n.concatenate(ws)
                    wgts = n.concatenate(wgts).flatten()
                    # Grid data into UV matrix
                    (us,vs,ws),ds,wgts = im.append_hermitian((us,vs,ws),ds,wgts)
                    im.put((us,vs,ws), ds, wgts)
                    us,vs,ws = [],[],[]
                    ds,wgts = [],[]

        # Grid remaining data into UV matrix
        sys.stdout.write('|\n'); sys.stdout.flush()
        ds = n.concatenate(ds)
        us = n.concatenate(us)
        vs = n.concatenate(vs)
        ws = n.concatenate(ws)
        wgts = n.concatenate(wgts).flatten()
        if len(ds) == 0: raise ValueError('No data to use.')
        (us,vs,ws),ds,wgts = im.append_hermitian((us,vs,ws),ds,wgts)
        im.put((us,vs,ws), ds, wgts)
        if opts.uniform > 0: im.uniform_wgt(thresh=opts.uniform)

        # Form dirty images/beams
        uvs = a.img.recenter(n.abs(im.uv).astype(n.float), (DIM/2,DIM/2))
        bms = a.img.recenter(n.abs(im.bm[0]).astype(n.float), (DIM/2,DIM/2))
        dim = im.image((DIM/2, DIM/2))
        dbm = im.bm_image(term=0, center=(DIM/2,DIM/2))
        if 'uvs' in outputs: to_fits('uvs', uvs)
        if 'bms' in outputs: to_fits('bms', bms)
        if 'dim' in outputs: to_fits('dim', dim)
        if 'dbm' in outputs: to_fits('dbm', dbm)

    dbm = a.img.recenter(dbm, (DIM/2+1,DIM/2))
    bm_gain = a.img.beam_gain(dbm)
    print 'Gain of dirty beam:', bm_gain

    # Deconvolve
    if opts.deconv == 'mem':
        cim,info = a.deconv.maxent_findvar(dim, dbm, f_var0=opts.var,
            maxiter=opts.maxiter, verbose=True, tol=opts.tol, maxiterok=True)
    elif opts.deconv == 'lsq':
        cim,info = a.deconv.lsq(dim, dbm, 
            maxiter=opts.maxiter, verbose=True, tol=opts.tol)
    elif opts.deconv == 'cln':
        cim,info = a.deconv.clean(dim, dbm, 
            maxiter=opts.maxiter, verbose=True, tol=opts.tol)
    elif opts.deconv == 'ann':
        cim,info = a.deconv.anneal(dim, dbm, maxiter=opts.maxiter, 
            cooling=lambda i,x: opts.tol*(1-n.cos(i/50.))*(x**2), verbose=True)
    else:
        cim,info = n.zeros_like(dim), {'res':dim}

    rim = info['res'] / bm_gain
    bim = rim + cim
    if 'cim' in outputs: to_fits('cim', cim)
    if 'rim' in outputs: to_fits('rim', rim)
    if 'bim' in outputs: to_fits('bim', bim)

