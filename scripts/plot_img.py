#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('plot_img.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-s', '--src', dest='src',
    help='The name of a source to phase to, or ra,dec position.')
o.add_option('-l', '--loc', dest='loc', 
    help='Use location-specific info for this location.')
o.add_option('-r', '--residual', dest='residual', action='store_true',
    help='Display only the residual in the 4th panel (otherwise display sum of clean image and residual).')
o.add_option('-d', '--deconv', dest='deconv', default='mem',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).')
o.add_option('--var', dest='var', type='float', default=1.,
    help='Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.  For annealing, interpreted as cooling speed.')
o.add_option('--maxiter', dest='maxiter', type='int', default=200,
    help='Number of allowable iterations per deconvolve attempt.')
o.add_option('--skip_amp', dest='skip_amp', action='store_true',
    help='Do not use amplitude information to normalize visibilities.')
o.add_option('--skip_bm', dest='skip_bm', action='store_true',
    help='Do not weight visibilities by the strength of the primary beam.')
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
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--dyn_rng', dest='dyn_rng', type='float', default=3.,
    help="Dynamic range in color of image (log).")
o.add_option('--smooth', dest='smooth', type='float', default=0.,
    help="Before deconvolving, smooth img and beam by a gaussian with width specified in pixels.")
o.add_option('--phs', dest='phs', action='store_true',
    help="If plotting the UV matrix, show phases.")
o.add_option('--buf_thresh', dest='buf_thresh', default=1.5e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')

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
assert(opts.pol in ['xx', 'yy', 'xy', 'yx'])
p1,p2 = opts.pol

uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
if opts.pol is None: active_pol = uv['pol']
else: active_pol = a.miriad.str2pol[opts.pol]
chans = gen_chans(opts.chan, uv)
del(uv)

aa.select_chans(chans)
if opts.src.find(',') == -1: src = a.src.get_src(opts.src)
else:
    ra, dec = opts.src.split(',')
    src = a.ant.RadioFixedBody(ra, dec)
if opts.no_w: im = a.img.Img(opts.size, opts.res)
else: im = a.img.ImgW(opts.size, opts.res)
DIM = int(opts.size/opts.res)

cnt, curtime = 0, None
# Gather data
uvw, dat, wgt = [], [], []
cache = {}
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
                top = a.coord.azalt2top((src.az,src.alt))
                cache = {}
        if cnt != 0: continue
        d = d.take(chans)
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
                    cache[j][p1] = r.flatten()
                # Calculate beam strength for weighting purposes
                w = cache[i][p1] * cache[j][p2]
                # For optimal SNR, down-weight data that is already attenuated 
                # by beam  by another factor of the beam response (modifying 
                # weight accordingly).
                #d *= w; w *= w
            else: w = n.ones(d.shape, dtype=n.float)
        except(a.ant.PointingError): continue
        valid = n.logical_not(d.mask)
        d = d.compress(valid).data
        if len(d) == 0: continue
        dat.append(d)
        uvw.append(xyz.compress(valid, axis=0))
        wgt.append(w.compress(valid))
        if len(dat) * len(chans) > opts.buf_thresh:
            sys.stdout.write('|'); sys.stdout.flush()
            dat = n.concatenate(dat)
            uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
            wgt = n.concatenate(wgt).flatten()
            uvw,dat,wgt = im.append_hermitian(uvw,dat,wgt)
            im.put(uvw, dat, wgt)
            uvw, dat, wgt = [], [], []

if len(uvw) == 0: raise ValueError('No data to plot')
# Grid data into UV matrix
sys.stdout.write('|\n'); sys.stdout.flush()
dat = n.concatenate(dat)
uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
wgt = n.concatenate(wgt).flatten()
uvw,dat,wgt = im.append_hermitian(uvw,dat,wgt)
im.put(uvw, dat, wgt)

im_img = im.image((DIM/2, DIM/2))
bm_img = im.bm_image()

if opts.smooth != 0:
    sm = a.img.gaussian_beam(opts.smooth, shape=im_img.shape)
    im_img = n.abs(a.img.convolve2d(im_img, sm))
    bm_img = n.abs(a.img.convolve2d(bm_img, sm))

if opts.deconv == 'mem':
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
if opts.deconv != 'none': rs_img = info['res']

p.subplot(221)
im_img = n.log10(im_img.clip(1e-15,n.Inf))
mx = im_img.max()
p.imshow(im_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
p.title('Dirty Image')

p.subplot(222)
bm_img = im.bm_image((DIM/2,DIM/2))
bm_gain = n.sqrt((bm_img**2).sum())
print 'Gain of dirty beam:', bm_gain
bm_img = n.log10(bm_img.clip(1e-15,n.Inf))
mx = bm_img.max()
p.imshow(bm_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
p.title('Dirty Beam')

p.subplot(223)
if opts.deconv != 'none':
    cl_img *= bm_gain
    if not opts.residual: rs_img += cl_img
    cl_img = n.log10(cl_img.clip(1e-15,n.Inf))
else:
    if opts.phs: cl_img = n.angle(im.uv)
    else: cl_img = n.log10(n.abs(im.uv))
mx = cl_img.max()
p.imshow(cl_img.copy(), vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
if opts.deconv != 'none': p.title('%s Image' % opts.deconv.upper())
else: p.title('UV Sampling')
if opts.deconv != 'none':
    # Generate a little info about where the strongest src is
    eq = im.get_eq(ra=src.ra, dec=src.dec, center=(DIM/2,DIM/2))
    ra,dec = a.coord.eq2radec(eq)
    print 'Phase center:', (src.ra, src.dec)
    # Print top 10 srcs (destructively for cl_img)
    for i in range(10):
        src_loc = cl_img.argmax()
        src_ra,src_dec = ra.flat[src_loc], dec.flat[src_loc]
        eq = ephem.Equatorial(src_ra, src_dec, epoch=aa.epoch)
        eq = ephem.Equatorial(eq, epoch=ephem.J2000)
        x,y = n.indices(cl_img.shape)
        print 'Src %2d:' % i, eq.get(), 'J2000,',
        print 'at pixel', (y.flat[src_loc], x.flat[src_loc])
        cl_img.flat[src_loc] = 0

p.subplot(224)
if opts.deconv != 'none':
    rs_img = n.log10(n.abs(rs_img).clip(1e-15,n.Inf))
else: rs_img = n.log10(n.abs(im.bm))
mx = rs_img.max()
p.imshow(rs_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
if opts.deconv != 'none': p.title('Residual Image')
else: p.title('Beam Sampling')

p.show()
