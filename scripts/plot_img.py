#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, sys, optparse

o = optparse.OptionParser()
o.set_usage('phs2src.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-s', '--src', dest='src',
    help='The source to phase to.')
o.add_option('-l', '--loc', dest='loc', 
    help='Use location-specific info for this location.')
o.add_option('-a', '--ants', dest='ants', default='all',
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

uv = aipy.miriad.UV(args[0])
aa = aipy.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'],
    use_bp=False)
if opts.pol is None: active_pol = uv['pol']
else: active_pol = a.miriad.str2pol[opts.pol]
chans = gen_chans(opts.chan, uv)
del(uv)

aa.select_chans(chans)
src = aipy.src.get_src(opts.src)
im = aipy.img.ImgW(opts.size, opts.res)
DIM = int(opts.size/opts.res)

cnt, curtime = 0, None
# Gather data
uvw, dat, wgt = [], [], []
for filename in sys.argv[1:]:
    sys.stdout.write('.'); sys.stdout.flush()
    uv = aipy.miriad.UV(filename)
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
            w = aa.ants[0].response((src.az, src.alt), pol=2)**2
        except(aipy.ant.PointingError): break
        valid = n.logical_not(d.mask)
        d = d.compress(clean).data
        if len(d) == 0: continue
        dat.append(d)
        uvw.append(xyz.compress(valid, axis=0))
        wgt.append(w.compress(valid))
        if len(data) * len(active_chans) > opts.buf_thresh:
            sys.stdout.write('|'); sys.stdout.flush()
            data = n.concatenate(data)
            uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
            wgt = n.concatenate(wgt).flatten()
            data /= wgt; wgt = wgt * 0 + 1
            uvw,data,wgt = im.append_hermitian(uvw,data,wgt)
            im.put(uvw, data, wgt)
            uvw, data, wgt = [], [], []

if len(uvw) == 0: raise ValueError('No data to plot')
# Grid data into UV matrix
sys.stdout.write('|\n'); sys.stdout.flush()
data = n.concatenate(data)
uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
wgt = n.concatenate(wgt).flatten()
# For optimal SNR, down-weight data which is already attenuated by beam 
# by another factor of the beam response (modifying weight accordingly).
#data *= wgt; wgt *= wgt
data /= wgt; wgt = wgt * 0 + 1
uvw,data,wgt = im.append_hermitian(uvw,data,wgt)
im.put(uvw, data, wgt)

im_img = im.image((DIM/2, DIM/2))
bm_img = im.bm_image()
# Try various noise levels until one works
cl_img, info = None, None
loop_num = 0
im_var = n.var(im_img)
if im_var == 0: raise ValueError('No flux in image')
while cl_img is None:
    print loop_num, im_var
    var0 = im_var / 10 / (2**loop_num)
    while var0 < im_var * 10 * (2**loop_num):
        print 'Trying var0=%f' % var0
        c, i = a.deconv.maxent(im_img, bm_img,
            var0=var0, maxiter=200, verbose=False, tol=1e-6)
        print 'Success =', i['success'],
        print 'Term:', i['term'], 'Score:', i['score']
        # Check if fit converged
        if i['success'] and i['term'] == 'tol':
            cl_img, info = c, i
            break
        else:
            if not cl_img is None: break
            var0 *= 2. ** (1./(loop_num+1))
    loop_num += 1
print 'Done with MEM.'
rs_img = info['res'] / n.abs(bm_img).sum()

pylab.subplot(221)
im_img = numpy.log10(im_img + 1e-15)
mx = im_img.max()
pylab.imshow(im_img, vmin=mx-3, vmax=mx)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Dirty Image')

pylab.subplot(222)
bm_img = numpy.log10(bm_img + 1e-15)
mx = bm_img.max()
pylab.imshow(bm_img, vmin=mx-3, vmax=mx)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Dirty Beam')

pylab.subplot(223)
cl_img = numpy.log10(cl_img + 1e-15)
mx = cl_img.max()
pylab.imshow(cl_img, vmin=mx-3, vmax=mx)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Clean Image')

pylab.subplot(224)
rs_img = numpy.log10(rs_img + 1e-15)
mx = rs_img.max()
pylab.imshow(rs_img, vmin=mx-3, vmax=mx)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Residual Image')

pylab.show()
