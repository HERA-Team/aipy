#! /usr/bin/env python
import sys, numpy as n, os, aipy as a
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
    help='Comma-delimited ranges (e.g. 1-3,5-10) of channels to manually use before any other statistical flagging.')
p.add_option('-m', '--map', dest='map', 
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')
opts, args = p.parse_args(sys.argv[1:])

# Get antenna array information
uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'], use_bp=False)
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
if os.path.exists(opts.map): skymap = a.img.SkyMap(fromfile=opts.map)
else: skymap = a.img.SkyMap(res=1/opts.size)

# Some imaging constants
dim = int(opts.size/opts.res)
pol = a.miriad.pol_code[opts.pol]

# Loop through RA and DEC, imaging on 15 degree grid
for ra1 in range(0,24):
  ra2 = 0
  for dec in range(-75, 90, 15):
    s = a.ant.RadioFixedBody('%2d:%02d:00' % (ra1, ra2), 
        '%2d:00:00' % dec, name='%2d:%02d, %3d' % (ra1, ra2, dec))
    print 'Pointing (ra, dec):', s.src_name
    cnt, curtime = 0, None
    uvw, data, wgt = [], [], []
    im = a.img.ImgW(opts.size, opts.res)
    for filename in args:
        sys.stdout.write('.'); sys.stdout.flush()
        uv = a.miriad.UV(filename)
        uv.select('auto', 0, 0, include=False)  # Only use cross-corrs
        uv.select('polarization', pol, 0, include=True) # Only use 1 pol
        for p,d in uv.all():
            uvw, t, (i,j) = p
            if curtime != t:
                curtime = t
                cnt +=1
            if cnt % opts.decimate != 0: continue
            aa.set_jultime(t)
            s.compute(aa)
            d = d.take(active_chans)
            try:
                d, xyz = aa.phs2src(d, s, i, j, with_coord=True)
                w = aa.ants[j].response((s.az, s.alt), pol=opts.pol[0]) * \
                    aa.ants[i].response((s.az, s.alt), pol=opts.pol[1])
            except(a.ant.PointingError): break # Finish file if below horiz
            clean = n.logical_not(d.mask)
            d = d.compress(clean).data
            if len(d) == 0: continue
            xyz = xyz.compress(clean, axis=0)
            w = w.compress(clean)
            data.append(d); uvw.append(xyz); wgt.append(w)
            # If buffer gets big, grid data to UV matrix.
            if len(data) * len(active_chans) > opts.buf_thresh:
                sys.stdout.write('|'); sys.stdout.flush()
                data = n.concatenate(data)
                uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
                wgt = n.concatenate(wgt).flatten()
                uvw,data,wgt = im.append_hermitian(uvw,data,wgt)
                im.put(uvw, data, wgt)
                uvw, data, wgt = [], [], []
    if len(uvw) == 0: continue
    # Grid data into UV matrix
    sys.stdout.write('|\n'); sys.stdout.flush()
    data = n.concatenate(data)
    uvw = n.concatenate(uvw); uvw.shape = (uvw.size / 3, 3)
    wgt = n.concatenate(wgt).flatten()
    # For optimal SNR, down-weight data which is already attenuated by beam 
    # by another factor of the beam response (modifying weight accordingly).
    data *= wgt; wgt *= wgt
    uvw,data,wgt = im.append_hermitian(uvw,data,wgt)
    im.put(uvw, data, wgt)
    im_img = im.image((dim/2, dim/2))
    bm_img = im.bm_image()
    var0 = n.var(im_img) * 10
    # Try various noise levels until one works
    cl_img, info = None, None
    while True:
        print 'Trying var0=%f' % var0
        c, i = a.deconv.maxent(im_img, bm_img, 
            var0=var0, maxiter=200, verbose=False, tol=1e-5)
        print 'Success =', i['success'], 
        print 'Term:', i['term'], 'Score:', i['score']
        # Check if fit converged
        if i['success'] and i['term'] == 'tol':
            cl_img, info = c, i
            var0 /= 4
        else:
            if not cl_img is None: break
            var0 *= 10
    print 'Done with MEM.'
    rs_img = info['res'] / n.abs(bm_img).sum()
    # Flip data around to go in increasing order of ra left->right
    # and increasing dec bottom->top
    cl_img = n.fliplr(cl_img); cl_img = n.flipud(cl_img)
    rs_img = n.fliplr(rs_img); rs_img = n.flipud(rs_img)
    img = cl_img + rs_img
    ras, decs = im.get_coords(s.ra, s.dec)
    # Weight data in skymap according to nearness to pointing center...
    wgts =a.img.gaussian_beam(dim/16, shape=im.uv.shape,
        center=(dim/2,dim/2))
    # and the strength of the beam in that direction.
    wgts *= n.abs(bm_img).sum()
    skymap.add_data(ras, decs, cl_img, wgts)
    skymap.tofile(opts.map)
skymap.tofile(opts.map)
