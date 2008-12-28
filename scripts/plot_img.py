#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, sys, optparse

o = optparse.OptionParser()
o.set_usage('phs2src.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-s', '--src', dest='src',
    help='The name of a source to phase to, or ra,dec position.')
o.add_option('-l', '--loc', dest='loc', 
    help='Use location-specific info for this location.')
o.add_option('-d', '--deconv', dest='deconv', action='store_true',
    help='Attempt to deconvolve the dirty image by the dirty beam.')
o.add_option('--var', dest='var', type='float', default=0.,
    help='Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.')
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
        if cnt != 0: continue
        aa.set_jultime(t)
        src.compute(aa)
        d = d.take(chans)
        try:
            d, xyz = aa.phs2src(d, src, i, j, with_coord=True)
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

if opts.deconv:
    # Try various noise levels until one works
    cl_img, info = None, None
    loop_num = -1
    if opts.var != 0: im_var = opts.var
    else: im_var = n.var(im_img)
    while cl_img is None:
        print '________________________________________________________'
        if loop_num == -1: var0 = im_var
        else: var0 = im_var / (1.5**loop_num)
        while loop_num < 0 or var0 < im_var * (1.5**loop_num):
            print 'Trying var0=%f' % var0
            c, i = a.deconv.maxent(im_img, bm_img,
                var0=var0, maxiter=opts.maxiter, verbose=False, tol=opts.tol)
            print 'Success =', i['success'],
            print 'Term:', i['term'], 'Score:', i['score']
            # Check if fit converged
            if i['success'] and i['term'] == 'tol':
                cl_img, info = c, i
                break
            else:
                if not cl_img is None: break
                if loop_num == -1: break
                var0 *= 1.2 ** (1./(2*(loop_num+1)))
        loop_num += 1
    print 'Done with MEM.'
    rs_img = info['res']

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
if opts.deconv: cl_img = n.log10(cl_img + 1e-15)
else: cl_img = n.abs(im.uv)
mx = cl_img.max()
p.imshow(cl_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
if opts.deconv: p.title('Clean Image')
else: p.title('UV Sampling')

p.subplot(224)
if opts.deconv: rs_img = n.log10(n.abs(rs_img) + 1e-15)
else: rs_img = n.abs(im.bm)
mx = rs_img.max()
p.imshow(rs_img, vmin=mx-opts.dyn_rng, vmax=mx, aspect='auto')
p.colorbar(shrink=.5, fraction=.05)
if opts.deconv: p.title('Residual Image')
else: p.title('Beam Sampling')

p.show()
