#! /usr/bin/env python
import sys, numpy, pylab, random, os, aipy

sfreq = 0.1126
sdf = 0.00234 / 2
nchan = 64
aa = aipy.loc.get_aa('pwa303', sdf, sfreq, nchan, use_bp=False)

chans = range(15,20) + range(23,64)
#chans = range(25,64)
aa.select_chans(chans)

size = float(100)
res = float(.5)
dim = int(size/res)

# Sky map in cylindrical coordinates
if os.path.exists('skymap'): skymap = aipy.img.SkyMap(fromfile='skymap')
else: skymap = aipy.img.SkyMap(res=1/size)

for ra1 in range(0,24):
#  for dec in range(-75, 90, 15):
  #for dec in range(-45, 75, 15):
    ra2 = 0
    srcs = []
    for dec in range(-75, 90, 15):
        srcs.append(aipy.ant.RadioFixedBody('%2d:%02d:00' % (ra1, ra2), 
            '%2d:00:00' % dec, name='%2d:%02d, %3d' % (ra1, ra2, dec)))
    #srcs.append(aipy.ant.RadioFixedBody('%2d:%02d:00' % (ra1, ra2), 
    #     '%2d:00:00' % dec, name='%2d:%02d, %3d' % (ra1, ra2, dec)))
    snames = [s.src_name for s in srcs]
    cat = aipy.ant.SrcCatalog(srcs)
    cnt = 0
    curtime = -1
    uvw, data, wgt = {}, {}, {}
    for s in cat:
        uvw[s] = []
        data[s] = []
        wgt[s] = []
    for filename in sys.argv[1:]:
        print filename
        uv = aipy.miriad.UV(filename)
        uv.select_data('auto', 0, 0, include_it=False)
        #uv.select_data('antennae', 0, 1, include_it=True)
        while True:
            p, d = uv.read_data()
            if d.size == 0: break
            t, bl = p[-2:]
            if curtime != t:
                curtime = t
                cnt = (cnt + 1) % 4
            if cnt != 0: continue
            aa.set_jultime(t)
            cat.compute(aa)
            d = d.take(chans)
            if numpy.any(d.mask): continue
            for s in cat:
                try:
                    ds, xyz = aa.phs2src(d, cat[s], bl, with_coord=True)
                    w = aa.ants[0].response((cat[s].az, cat[s].alt), pol=2)**2
                except(aipy.ant.PointingError): continue
                data[s].append(ds)
                uvw[s].append(xyz)
                wgt[s].append(w)
    for n, s in enumerate(snames):
        print cat[s].ra, cat[s].dec
        if len(uvw[s]) == 0: continue
        im = aipy.img.ImgW(size, res)
        #im = aipy.img.Img(size, res)
        im_uvw = numpy.array(uvw[s])
        im_data = numpy.array(data[s]).flatten()
        im_wgt = numpy.array(wgt[s]).flatten()
        im_uvw.shape = (im_uvw.size / 3, 3)
        im_uvw,im_data,im_wgt = im.append_hermitian(im_uvw,im_data,im_wgt)
        im.put(im_uvw, im_data, im_wgt)
        im.uv /= len(chans)      # Put back in Jys
        im_img = im.image((dim/2, dim/2))
        bm_img = im.bm_image()
        var0 = numpy.var(im_img)
        cl_img, res_img = aipy.deconv.maxent(im_img, bm_img, var0=var0, 
            maxiter=200, verbose=True)
        if res_img.max() > im_img.max(): continue
        cl_img = numpy.fliplr(cl_img)
        cl_img = numpy.flipud(cl_img)
        ras, decs = im.get_coords(cat[s].ra, cat[s].dec)
        wgts =aipy.img.gaussian_beam(dim/16, shape=im.uv.shape,
            center=(dim/2,dim/2))
        wgts *= bm_img.sum()
        skymap.add_data(ras, decs, cl_img, wgts)
        skymap.tofile('skymap')
skymap.tofile('skymap')
