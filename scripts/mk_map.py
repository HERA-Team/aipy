#! /usr/bin/env python
"""
This is a general-purpose script for making faceted sky maps from MIRIAD UV 
files.  Data (optionally selected for baseline, channel) are read from the 
file, phased to each of 264 phase centers on a 15 degree grid, normalized for 
passband/primary beam effects, gridded to a UV matrix, imaged, and optionally 
deconvolved by a corresponding PSF to produce a clean image.  The clean image
and residual are then combined and gridded onto a Healpix map with a weighting
dictated by the distance from phase center.  At each iteration, the resulting
map is saved to the output fits file.

Author: Aaron Parsons
"""

import sys, numpy as n, os, aipy as a, optparse

o = optparse.OptionParser()
o.set_usage('mk_map.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, 
    loc=True, dec=True)
o.add_option('--no_res', dest='no_res', action='store_true',
    help='Do not include the residual data (after deconvolution) in the map.')
o.add_option('-d', '--deconv', dest='deconv', default='mem',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).  Default "mem".')
o.add_option('--var', dest='var', type='float', default=.6,
    help='Starting variance guess for MEM deconvolution, as a fraction of the total image variance.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.')
o.add_option('--maxiter', dest='maxiter', type='int', default=200,
    help='Number of allowable iterations per deconvolve attempt.')
o.add_option('--skip_amp', dest='skip_amp', action='store_true',
    help='Do not use amplitude information to normalize visibilities.')
o.add_option('--skip_bm', dest='skip_bm', action='store_true',
    help='Do not weight visibilities by the strength of the primary beam.')
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data to healpix map.')
o.add_option('--size', dest='size', type='int', default=200,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--buf_thresh', dest='buf_thresh', default=1.8e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')
o.add_option('-m', '--map', dest='map',
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
a.scripting.uv_selector(uv, opts.ant, opts.pol)
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
DIM = int(opts.size/opts.res)
p1,p2 = opts.pol
del(uv)

# Open skymap
if os.path.exists(opts.map): skymap = a.map.Map(fromfits=opts.map)
else: skymap = a.map.Map(nside=256)
skymap.set_interpol(opts.interpolate)

# Define pointing centers of facets
NPTS = 100
ras,decs = a.map.facet_centers(NPTS, ncrd=2)

for i, (ra,dec) in enumerate(zip(ras,decs)):
    src = a.ant.RadioFixedBody(ra, dec)
    print '%d / %d' % (i + 1, NPTS)
    print 'Pointing (ra, dec):', src._ra, src._dec
    if abs(aa.lat - src._dec) > n.pi/2:
        print '    Source never rises: skipping'
        continue
    cnt, curtime = 0, None
    uvw, dat, wgt = [], [], []
    cache = {}
    if not opts.no_w: im = a.img.ImgW(opts.size, opts.res)
    else: im = a.img.Img(opts.size, opts.res)
    for filename in args:
        sys.stdout.write('.'); sys.stdout.flush()
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if curtime != t:
                curtime = t
                cnt = (cnt + 1) % opts.decimate
                if cnt == 0:
                    aa.set_jultime(t)
                    src.compute(aa)
                    top = a.coord.azalt2top((src.az,src.alt))
                    cache = {}
            if cnt != 0: continue
            d = d.take(chans)
            f = f.take(chans)
            if not opts.skip_amp:
                d /= aa.ants[i].gain * n.conjugate(aa.ants[j].gain)
            try:
                d = aa.phs2src(d, src, i, j)
                xyz = aa.gen_uvw(i,j,src=src)
                if not opts.skip_bm:
                    # Cache beam response, since it is an expensive operation
                    if not cache.has_key(i): cache[i] = {}
                    if not cache[i].has_key(p1):
                        r = aa.ants[i].bm_response(top, pol=p1)
                        cache[i][p1] = r.flatten()
                    if not cache.has_key(j): cache[j] = {}
                    if not cache[j].has_key(p2):
                        r = aa.ants[j].bm_response(top, pol=p2)
                        cache[j][p2] = r.flatten()
                    # Calculate beam strength for weighting purposes
                    w = cache[i][p1] * cache[j][p2]
                    # For optimal SNR, down-weight data that is already
                    # attenuated  by beam  by another factor of the beam 
                    # response (modifying  weight accordingly).
                    d *= w; w *= w
                else: w = n.ones(d.shape, dtype=n.float)
            except(a.ant.PointingError): continue
            valid = n.logical_not(f)
            d = d.compress(valid)
            if len(d) == 0: continue
            dat.append(d)
            uvw.append(xyz.compress(valid, axis=0))
            wgt.append(w.compress(valid))
            # If buffer gets big, grid data to UV matrix.
            if len(dat) * len(chans) > opts.buf_thresh:
                sys.stdout.write('|'); sys.stdout.flush()
                dat = n.concatenate(dat)
                uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
                wgt = n.concatenate(wgt).flatten()
                uvw,dat,wgt = im.append_hermitian(uvw,dat,wgt)
                im.put(uvw, dat, wgt)
                uvw, dat, wgt = [], [], []
    if len(uvw) == 0:
        print ' never up.'
        continue
    # Grid data into UV matrix
    sys.stdout.write('|\n'); sys.stdout.flush()
    dat = n.concatenate(dat)
    uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
    wgt = n.concatenate(wgt).flatten()
    uvw,dat,wgt = im.append_hermitian(uvw,dat,wgt)
    im.put(uvw, dat, wgt)

    im_img = im.image((DIM/2, DIM/2))
    bm_img = im.bm_image()
    bm_gain = n.sqrt((bm_img**2).sum())
    print 'Gain of dirty beam:', bm_gain

    if opts.deconv == 'none': img = im_img / bm_gain
    elif opts.deconv == 'mem':
        cl_img,info = a.deconv.maxent_findvar(im_img, bm_img, f_var0=opts.var,
            maxiter=opts.maxiter, verbose=True, tol=opts.tol, maxiterok=True)
    elif opts.deconv == 'lsq':
        cl_img,info = a.deconv.lsq(im_img, bm_img,
            maxiter=opts.maxiter, verbose=True, tol=opts.tol)
    elif opts.deconv == 'cln':
        cl_img,info = a.deconv.clean(im_img, bm_img,
            maxiter=opts.maxiter, verbose=True, tol=opts.tol)
    elif opts.deconv == 'ann':
        cl_img,info = a.deconv.anneal(im_img, bm_img, maxiter=opts.maxiter,
            cooling=lambda i,x: opts.tol*(1-n.cos(i/50.))*(x**2), verbose=True)
    if opts.deconv != 'none':
        rs_img = info['res'] / bm_gain
        if not opts.no_res: img = cl_img + rs_img
        else: img = cl_img
    # Get coordinates of image pixels in original (J2000) epoch
    ex,ey,ez = im.get_eq(src._ra, src._dec, center=(DIM/2,DIM/2))
    tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
    # Define a weighting for gridding data into the skymap
    map_wgts = n.exp(-(tx**2 + ty**2) / .1**2)
    map_wgts.shape = (map_wgts.size,)
    valid = n.logical_not(map_wgts.mask)
    ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
    map_wgts = map_wgts.compress(valid)
    img = img.flatten()
    img = img.compress(valid)
    # Put the data into the skymap
    skymap.add((ex,ey,ez), map_wgts, img)
    skymap.to_fits(opts.map, clobber=True)

