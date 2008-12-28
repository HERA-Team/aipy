#! /usr/bin/env python
import sys, numpy, pylab, aipy

sfreq = 0.1126
sdf = 0.00234 / 2
nchan = 64
aa = aipy.loc.get_aa('pwa303', sdf, sfreq, nchan, use_bp=False)

chans = range(14,20) + range(24,64)
aa.select_chans(chans)
cat = aipy.src.get_catalog()
size = 100
res = .5
src = aipy.src.get_src('cyg')
#src = aipy.ant.RadioFixedBody('20:00','40')
#src = aipy.ant.RadioFixedBody('4:00','-15')

dim = int(size/res)

im = aipy.img.ImgW(size, res)
cnt = 0
curtime = -1
# Gather data
uvw, dat, wgt = [], [], []
for filename in sys.argv[1:]:
    print filename
    uv = aipy.miriad.UV(filename)
    uv.select('auto', 0, 0, include=False)
    for p,d in uv.all():
        uvw, t, (i,j) = p
        if curtime != t:
            curtime = t
            cnt = (cnt + 1) % 100
        if cnt != 0: continue
        d = d.take(chans)
        if numpy.any(d.mask): continue
        aa.set_jultime(t)
        src.compute(aa)
        try:
            d, xyz = aa.phs2src(d, src, i, j, with_coord=True)
            w = aa.ants[0].response((src.az, src.alt), pol=2)**2
        except(aipy.ant.PointingError): break
        dat.append(d)
        uvw.append(xyz)
        wgt.append(w)
uvw, dat, wgt = numpy.array(uvw), numpy.array(dat), numpy.array(wgt)
uvw.shape = (uvw.size / 3, 3)
dat = dat.flatten()
wgt = wgt.flatten()
uvw, dat, wgt = im.append_hermitian(uvw, dat, wgt)
im.put(uvw, dat, wgt)

bm_width = 1.2
eff_bm = aipy.img.gaussian_beam(bm_width, shape=im.uv.shape)

im_im = im.image((dim/2, dim/2))
im_bm = im.bm_image()
eff_bm /= eff_bm.sum() / im_bm.sum()
#im_cl, im_rs = aipy.img.clean(im_im, im_bm, gain=.1, verbose=True)
im_cl, im_rs = aipy.deconv.maxent(im_im, im_bm, var0=numpy.var(im_im),maxiter=200, verbose=True)
#im_cl, im_rs = aipy.deconv.lsq(im_im, im_bm, tol=1e-4, verbose=True)
im_sm = aipy.img.convolve2d(im_cl, eff_bm).real
im_rs = numpy.abs(im_rs)

pylab.subplot(221)
im_im = numpy.log10(im_im + 1e-15)
mx = im_im.max()
pylab.imshow(im_im, vmin=mx-3, vmax=mx)
#pylab.imshow(im_im)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Dirty')

pylab.subplot(222)
im_cl = numpy.log10(im_cl + 1e-15)
mx = im_cl.max()
pylab.imshow(im_cl, vmin=mx-3, vmax=mx)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Clean')

pylab.subplot(223)
im_rs = numpy.log10(im_rs + 1e-15)
mx = im_rs.max()
pylab.imshow(im_rs, vmin=mx-3, vmax=mx)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Residual')

pylab.subplot(224)
im_sm = numpy.log10(im_sm + 1e-15)
mx = im_im.max()
pylab.imshow(im_sm, vmin=mx-3, vmax=mx)
pylab.colorbar(shrink=.5, fraction=.05)
pylab.title('Smooth')

pylab.show()
