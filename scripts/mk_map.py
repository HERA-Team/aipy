#! /usr/bin/env python
import sys, numpy as n, os, aipy as a, random, optparse

o = optparse.OptionParser()
o.set_usage('mk_map.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-l', '--loc', dest='loc',
    help='Use location-specific info for this location.')
o.add_option('-d', '--deconv', dest='deconv', default='mem',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).')
o.add_option('-r', '--residual', dest='residual', action='store_true',
    help='Include the residual data (after deconvolution) in the map.')
o.add_option('--var', dest='var', type='float', default=.8,
    help='Starting variance guess for deconvolution, as a fraction of the total image variance.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.')
o.add_option('--maxiter', dest='maxiter', type='int', default=200,
    help='Number of allowable iterations per deconvolve attempt.')
o.add_option('-a', '--amp', dest='amp', action='store_true',
    help='Use amplitude information to normalize visibilities.')
o.add_option('-w', '--wgt_bm', dest='wgt_bm', action='store_true',
    help='Weight visibilities by the strength of the primary beam.')
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data to healpix map.')
o.add_option('-b', '--baselines', dest='baselines', default='cross',
    help='Select which antennas/baselines to include in plot.  Options are: "all", "auto", "cross", "<ant1 #>_<ant2 #>" (a specific baseline), or "<ant1 #>,..." (a list of active antennas).')
o.add_option('-p', '--pol', dest='pol', default=None,
    help='Choose which polarization (xx, yy, xy, yx) to plot.')
o.add_option('-c', '--chan', dest='chan', default='all',
    help='Select which channels (taken after any delay/fringe transforms) to plot.  Options are: "all", "<chan1 #>,..." (a list of active channels), or "<chan1 #>_<chan2 #>" (a range of channels).  If "all" or a range are selected, a 2-d image will be plotted.  If a list of channels is selected, an xy plot will be generated.')
o.add_option('-x', '--decimate', dest='decimate', default=1, type='int',
    help='Take every Nth time data point.')
o.add_option('--size', dest='size', type='int', default=200,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--buf', dest='buf_thresh', default=1.8e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')
o.add_option('-m', '--map', dest='map',
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')

def convert_arg_range(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    arg = arg.split(',')
    return [map(float, option.split('_')) for option in arg]

def data_selector(antopt, active_pol, uv):
    """Call uv.select with appropriate options based on command-line argument 
    for ants."""
    if antopt.startswith('all'): pass
    elif antopt.startswith('auto'): uv.select('auto', 0, 0)
    elif antopt.startswith('cross'): uv.select('auto', 0, 0, include=0)
    else:
        antopt = convert_arg_range(antopt)
        for opt in antopt:
            try: a1,a2 = opt
            except(ValueError): a1,a2 = opt + [-1]
            uv.select('antennae', a1, a2)
    uv.select('polarization', active_pol, 0)

def gen_chans(chanopt, uv):
    """Return an array of active channels based on command-line arguments."""
    if chanopt == 'all': chans = n.arange(uv['nchan'])
    else:
        chanopt = convert_arg_range(chanopt)
        if len(chanopt[0]) != 1:
            chanopt = [n.arange(x,y, dtype=n.int) for x,y in chanopt]
        chans = n.concatenate(chanopt)
    return chans.astype(n.int)

opts, args = o.parse_args(sys.argv[1:])

# Get antenna array information
uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
if opts.pol is None: active_pol = uv['pol']
else: active_pol = a.miriad.str2pol[opts.pol]
chans = gen_chans(opts.chan, uv)
del(uv)

aa.select_chans(chans)

# Open skymap
if os.path.exists(opts.map): skymap = a.map.Map(fromfits=opts.map)
else: skymap = a.map.Map(nside=256)
skymap.set_interpol(opts.interpolate)

# Some imaging constants
DIM = int(opts.size/opts.res)

# Loop through RA and DEC, imaging on 15 degree grid
ras_decs = []
for ra1 in range(0,24):
  ra2 = 0
  for dec in range(-75, 90, 15):
    ras_decs.append((ra1, ra2, dec))
random.shuffle(ras_decs)

for i, (ra1, ra2, dec) in enumerate(ras_decs):
    src = a.ant.RadioFixedBody('%2d:%02d:00' % (ra1, ra2),
        '%2d:00:00' % dec, name='%2d:%02d, %3d' % (ra1, ra2, dec))
    print '%d / %d' % (i + 1, len(ras_decs))
    print 'Pointing (ra, dec):', src.src_name
    cnt, curtime = 0, None
    uvw, dat, wgt = [], [], []
    im = a.img.ImgW(opts.size, opts.res)
    for filename in args:
        sys.stdout.write('.'); sys.stdout.flush()
        uv = a.miriad.UV(filename)
        data_selector(opts.baselines, active_pol, uv)
        for (crd,t,(i,j)),d in uv.all():
            if curtime != t:
                curtime = t
                cnt = (cnt + 1) % opts.decimate
                if cnt == 0:
                    aa.set_jultime(t)
                    src.compute(aa)
            if cnt != 0: continue
            d = d.take(chans)
            if opts.amp: d /= aa.ants[i].gain * aa.ants[j].gain
            try:
                d = aa.phs2src(d, src, i, j)
                xyz = aa.gen_uvw(i,j,src=src)
                if opts.wgt_bm:
                    src_top = a.coord.azalt2top((src.az,src.alt))
                    w = aa.ants[i].response(src_top, pol=opts.pol[0]) *\
                        aa.ants[j].response(src_top, pol=opts.pol[1])
                    w = w.flatten()
# For optimal SNR, down-weight data that is already attenuated by beam 
# by another factor of the beam response (modifying weight accordingly).
                    #d *= w; w *= w
                else: w = n.ones(d.shape, dtype=n.float)
            except(a.ant.PointingError): continue
            valid = n.logical_not(d.mask)
            d = d.compress(valid).data
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

    if opts.deconv == 'none': img = im_img #/ n.sqrt((bm_img**2).sum())
    elif opts.deconv == 'mem':
        cl_img,info = a.deconv.maxent_findvar(im_img, bm_img, f_var0=opts.var,
            maxiter=opts.maxiter, verbose=True, tol=opts.tol)
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
        rs_img = info['res']
        if opts.residual: img = n.abs(bm_gain*cl_img + rs_img)
        else: img = bm_gain * cl_img
    ex,ey,ez = im.get_eq(src.ra, src.dec, center=(DIM/2,DIM/2))
    #crds.shape = (3, crds.size/3)
    #crds = crds.transpose()
    tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
    map_wgts = n.exp(-(tx**2 + ty**2) / .1**2)
    map_wgts.shape = (map_wgts.size,)
    valid = n.logical_not(map_wgts.mask)
    #crds = crds.compress(valid, axis=0)
    ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
    map_wgts = map_wgts.compress(valid)
    img = img.flatten()
    img = img.compress(valid)
    #crds = n.asarray(crds)
    #print crds.shape, map_wgts.shape, img.shape
    #print skymap.wgt[crds].shape
    #skymap.add(crds, map_wgts, img)
    skymap.add((ex,ey,ez), map_wgts, img)
    skymap.to_fits(opts.map, clobber=True)

