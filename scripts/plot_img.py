#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('plot_img.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-s', '--src', dest='src',
    help='The name of a source to phase to, or ra,dec position.')
o.add_option('-l', '--loc', dest='loc', 
    help='Use location-specific info for this location.')
o.add_option('-d', '--deconv', dest='deconv', default='mem',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).')
o.add_option('--var', dest='var', type='float', default=1.,
    help='Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.  For annealing, interpreted as cooling speed.')
o.add_option('--maxiter', dest='maxiter', type='int', default=200,
    help='Number of allowable iterations per deconvolve attempt.')
o.add_option('-a', '--ants', dest='ants', default='cross',
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
o.add_option('--buf_thresh', dest='buf_thresh', default=1e6, type='float',
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

uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'],
    use_bp=False)
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
for filename in args:
    sys.stdout.write('.'); sys.stdout.flush()
    uv = a.miriad.UV(filename)
    data_selector(opts.ants, active_pol, uv)
    for (crd,t,(i,j)),d in uv.all():
        if curtime != t:
            curtime = t
            cnt = (cnt + 1) % opts.decimate
            if cnt == 0:
                aa.set_jultime(t)
                src.compute(aa)
        if cnt != 0: continue
        d = d.take(chans)
        try:
            d = aa.phs2src(d, src, i, j)
            xyz = aa.gen_uvw(i,j,src=src)
            #w = aa.ants[0].response((src.az, src.alt), pol=2)**2
            w = n.ones(d.shape, dtype=n.float)
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
            dat /= wgt; wgt = wgt * 0 + 1
            uvw,dat,wgt = im.append_hermitian(uvw,dat,wgt)
            im.put(uvw, dat, wgt)
            uvw, dat, wgt = [], [], []

if len(uvw) == 0: raise ValueError('No data to plot')
# Grid data into UV matrix
sys.stdout.write('|\n'); sys.stdout.flush()
dat = n.concatenate(dat)
uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
wgt = n.concatenate(wgt).flatten()
# For optimal SNR, down-weight data which is already attenuated by beam 
# by another factor of the beam response (modifying weight accordingly).
#data *= wgt; wgt *= wgt
dat /= wgt; wgt = wgt * 0 + 1
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
if opts.deconv != None: rs_img = info['res']

p.subplot(221)
im_img = n.log10(im_img + 1e-15)
mx = im_img.max()
p.imshow(im_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
p.title('Dirty Image')

p.subplot(222)
bm_img = im.bm_image((DIM/2,DIM/2))
bm_img = n.log10(bm_img + 1e-15)
mx = bm_img.max()
p.imshow(bm_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
p.title('Dirty Beam')

p.subplot(223)
if opts.deconv != None:
    cl_img = n.log10(cl_img + 1e-15)
    # Generate a little info about where the strongest src is
    src_loc = cl_img.argmax()
    eq = im.get_eq(ra=src.ra, dec=src.dec, center=(DIM/2,DIM/2))
    ra,dec = a.coord.eq2radec(eq)
    ra,dec = ra.flat[src_loc], dec.flat[src_loc]
    eq = ephem.Equatorial(ra, dec)
    print 'Phase center:', (src.ra, src.dec)
    x,y = n.indices(cl_img.shape)
    print 'Max src:', eq.get(), 'at pixel', (y.flat[src_loc], x.flat[src_loc])
    # ... and then finish plotting
else: cl_img = n.log10(n.abs(im.uv))
mx = cl_img.max()
p.imshow(cl_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
if opts.deconv != None: p.title('%s Image' % opts.deconv.upper())
else: p.title('UV Sampling')

p.subplot(224)
if opts.deconv != None: rs_img = n.log10(n.abs(rs_img) + 1e-15)
else: rs_img = n.log10(n.abs(im.bm))
mx = rs_img.max()
p.imshow(rs_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
if opts.deconv != None: p.title('Residual Image')
else: p.title('Beam Sampling')

p.show()
