#! /usr/bin/env python
import sys, numpy, os, aipy
from optparse import OptionParser

p = OptionParser()
p.set_usage('make_skymap.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
p.add_option('-x', '--decimate', dest='decimate', default=1, type='int',
    help='Only use every Nth integration.')
p.add_option('-s', '--size', dest='size', default=100., type='float',
    help='Size (in wavelengths) of UV matrix (actual matrix dimensions are size/resolution).  Default 100.')
p.add_option('-r', '--res', dest='res', default=.5, type='float',
    help='Resolution (in wavelengths) of UV matrix (actual matrix dimensions are size/resolution).  Default .5')
p.add_option('-p', '--pol', dest='pol', default='xx',
    help='The polarization to map (xx, yy, xy, yx).  Default xx')
p.add_option('-b', '--buf_thresh', dest='buf_thresh', default=1e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')
p.add_option('-c', '--chans', dest='chans', default='',
    help='Comma-delimited ranges (e.g. 1-3,5-10) of channels to manually flag before any other statistical flagging.')
p.add_option('-m', '--map', dest='map', 
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')
opts, args = p.parse_args(sys.argv[1:])

# Get antenna array information
uv = aipy.miriad.UV(args[0])
aa = aipy.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'],use_bp=False)
del(uv)

# Select channels to use in images
chans = opts.chans.split(',')
active_chans = []
for c in chans:
    if len(c) == 0: continue
    c = map(int, c.split('-'))
    active_chans.extend(range(c[0], c[1]+1))
aa.select_chans(active_chans)

# Open skymap (natively in cylindrical coordinates)
if os.path.exists(opts.map): skymap = aipy.img.SkyMap(fromfile=opts.map)
else: skymap = aipy.img.SkyMap(res=1/opts.size)

# Some imaging constants
dim = int(opts.size/opts.res)
pol = aipy.miriad.pol_code[opts.pol]

# Loop through RA and DEC, imaging on 15 degree grid
for ra1 in range(0,24):
  ra2 = 0
  for dec in range(-75, 90, 15):
    s = aipy.ant.RadioFixedBody('%2d:%02d:00' % (ra1, ra2), 
        '%2d:00:00' % dec, name='%2d:%02d, %3d' % (ra1, ra2, dec))
    print s.name
    cnt, curtime = 0, None
    uvw, data, wgt = [], [], []
    im = aipy.img.ImgW(opts.size, opts.res)
    for filename in args:
        print filename
        uv = aipy.miriad.UV(filename)
        uv.select_data('auto', 0, 0, include_it=False)  # Only use cross-corrs
        uv.select_data('polarization', pol, 0, include_it=True) # Only use 1 pol
        while True:
            p, d = uv.read_data()
            if d.size == 0: break
            t, bl = p[-2:]
            i, j = aipy.miriad.bl2ij(bl)
            if curtime != t:
                curtime = t
                cnt +=1
            if cnt % opts.decimate != 0: continue
            aa.set_jultime(t)
            s.compute(aa)
            d = d.take(active_chans)
            try:
                d, xyz = aa.phs2src(d, s, bl, with_coord=True)
                w = aa.ants[j].response((s.az, s.alt), pol=opts.pol[0]) * \
                    aa.ants[i].response((s.az, s.alt), pol=opts.pol[1])
            except(aipy.ant.PointingError): break # Finish file if below horiz
            clean = numpy.logical_not(d.mask)
            d = d.compress(clean).data
            if len(d) == 0: continue
            xyz = xyz.compress(clean, axis=0)
            w = w.compress(clean)
            data.append(d); uvw.append(xyz); wgt.append(w)
            # If buffer gets big, go ahead and grid data to UV matrix.
            if len(data) * len(active_chans) > opts.buf_thresh:
                print 'Gridding...'
                data = numpy.concatenate(data)
                uvw = numpy.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
                wgt = numpy.concatenate(wgt).flatten()
                uvw,data,wgt = im.append_hermitian(uvw,data,wgt)
                im.put(uvw, data, wgt)
                uvw, data, wgt = [], [], []
    if len(uvw) == 0: continue
    # Grid data into UV matrix
    print 'Gridding...'
    data = numpy.concatenate(data)
    uvw = numpy.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
    wgt = numpy.concatenate(wgt).flatten()
    # Fix weighting based on noise for optimal SNR
    # Downweight data which is already attenuated by beam with an additional
    # factor of the beam response, and modify weight accordingly.
    data *= wgt
    wgt *= wgt
    uvw,data,wgt = im.append_hermitian(uvw,data,wgt)
    im.put(uvw, data, wgt)
    im_img = im.image((dim/2, dim/2))
    bm_img = im.bm_image()
    var0 = numpy.var(im_img)
    cl_img, res_img = aipy.deconv.maxent(im_img, bm_img, var0=var0, 
        maxiter=200, verbose=True, tol=1e-5)
    # Check if fit converged
    if res_img.max() > im_img.max(): continue
    # Flip data around to go in increasing order of ra left->right
    # and increasing dec bottom->top
    cl_img = numpy.fliplr(cl_img); cl_img = numpy.flipud(cl_img)
    ras, decs = im.get_coords(s.ra, s.dec)
    # Weight data in skymap according to nearness to pointing center...
    wgts =aipy.img.gaussian_beam(dim/16, shape=im.uv.shape,
        center=(dim/2,dim/2))
    # and the strength of the beam in that direction.
    wgts *= numpy.abs(bm_img).sum()
    skymap.add_data(ras, decs, cl_img, wgts)
    skymap.tofile(opts.map)
skymap.tofile(opts.map)
