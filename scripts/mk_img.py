#! /usr/bin/env python
"""
This is a general-purpose script for making images from MIRIAD UV files.  Data
(optionally selected for baseline, channel) are read from the file, phased
to a provided position, normalized for passband/primary beam effects, gridded
to a UV matrix, imaged, and optionally deconvolved by a corresponding PSF to
produce a clean image.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, optparse, ephem, os

o = optparse.OptionParser()
o.set_usage('mk_img.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, loc2=True,
    src=True, dec=True)
o.add_option('-o', '--output', dest='output', default='dim,dbm',
    help='Comma delimited list of data to generate FITS files for.  Can be: dim (dirty image), dbm (dirty beam), uvs (uv sampling), or bms (beam sampling).  Default is dim,dbm.')
o.add_option('--npts', dest='npts', type='int', default=200,
    help='If no src is provided, facet the sphere into this many pointings for making a map.  Default 200.')
o.add_option('--snap', dest='snap', type='int', default=-1,
    help='Number of integrations to use in "snapshot" images.  Default is to not do snapshoting (i.e. all integrations go into one image).')
o.add_option('--cnt', dest='cnt', type='int', default=0,
    help='Start counting output images from this number.  Default 0.')
o.add_option('--fmt', dest='fmt', default='im%04d',
    help='A format string for counting successive images written to files.  Default is im%04d (i.e. im0001).')
o.add_option('-u', '--uniform', dest='uniform', type='float', default=0,
    help="Use uniform (rather than natural) weighting for uv bins that have a weight above the specified fraction of the maximum weighting.")
o.add_option('--skip_amp', dest='skip_amp', action='store_true',
    help='Do not use amplitude information to normalize visibilities.')
o.add_option('--skip_bm', dest='skip_bm', action='store_true',
    help='Do not weight visibilities by the strength of the primary beam.')
o.add_option('--size', dest='size', type='int', default=300,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--buf_thresh', dest='buf_thresh', default=2e6, type='float',
    help='Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
locs = a.scripting.files_to_locs(opts.loc, args, sys.argv)
aas = {}
for L in locs:
    uv = a.miriad.UV(locs[L][0])
    (j,t,j),j = uv.read()
    chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    aa = a.loc.get_aa(L, uv['sdf'], uv['sfreq'], uv['nchan'])
    aa.select_chans(chans)
    aas[L] = aa
    afreqs = aa.ants[0].beam.afreqs
    cfreq = n.average(afreqs)
    aa.set_jultime(t)
    del(uv)
outputs = opts.output.split(',')
if opts.no_w: Img = a.img.Img
else: Img = a.img.ImgW

# Get all sources that will be used as phase centers.  If no sources are
# specified, define phase centers for faceting a sphere.
if opts.src == 'zen':
    srcs = [a.ant.RadioFixedBody(aa.sidereal_time(), aa.lat, name='zen')]
    cat = a.ant.SrcCatalog(srcs)
elif not opts.src is None: 
    cat = a.scripting.parse_srcs(opts.src, force_cat=True)
else:
    ras,decs = a.map.facet_centers(opts.npts, ncrd=2)
    srcs = [a.ant.RadioFixedBody(ra,dec,name=str(i)) 
        for i,(ra,dec) in enumerate(zip(ras,decs))]
    cat = a.ant.SrcCatalog(srcs)

# Generate the image object that will be used.
us,vs,ws,ds,wgts = [],[],[],[],[]
im = Img(opts.size, opts.res, mf_order=0)
DIM = int(opts.size/opts.res)

# Define a quick function writing an image to a FITS file
def fname(ftag, cnt): return '%s.%s.fits' % (opts.fmt % cnt, ftag)
def to_fits(ftag,i,src,cnt):
    filename = fname(ftag,cnt)
    print 'Saving data to', filename
    while len(i.shape) < 4: i.shape = i.shape + (1,)
    cen = ephem.Equatorial(src.ra, src.dec, epoch=aa.epoch)
    cen = ephem.Equatorial(cen, epoch=ephem.J2000)
    L,M = im.get_LM()
    a.img.to_fits(filename, i, clobber=True,
        object=src.src_name, obs_date=str(aa.date),
        ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
        d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[1,1]*a.img.rad2deg,
        freq=n.average(aa.ants[0].beam.afreqs))

def grid_it(im,us,vs,ws,ds,wgts):
    sys.stdout.write('|'); sys.stdout.flush()
    if len(ds) == 0: raise ValueError('No data to use.')
    ds,wgts = n.concatenate(ds), n.concatenate(wgts).flatten()
    us,vs,ws = n.concatenate(us), n.concatenate(vs), n.concatenate(ws)
    # Grid data into UV matrix
    (us,vs,ws),ds,wgts = im.append_hermitian((us,vs,ws),ds,wgts)
    im.put((us,vs,ws), ds, wgts)

def img_it(im):
    if opts.uniform > 0: im.uniform_wgt(thresh=opts.uniform)
    # Form dirty images/beams
    uvs = a.img.recenter(n.abs(im.uv).astype(n.float), (DIM/2,DIM/2))
    bms = a.img.recenter(n.abs(im.bm[0]).astype(n.float), (DIM/2,DIM/2))
    dim = im.image((DIM/2, DIM/2))
    dbm = im.bm_image(term=0, center=(DIM/2,DIM/2))
    return uvs,bms, dim,dbm

# Loop through all specified sources, generating images
imgcnt = opts.cnt
for srccnt, s in enumerate(cat.values()):
    s.compute(aa)
    print '%d / %d' % (srccnt + 1, len(cat.values()))
    print 'Pointing (ra, dec):', s._ra, s._dec
    if abs(aa.lat - s._dec) > n.pi/2:
        print '    Source never rises: skipping'
        continue
    src = a.fit.SrcCatalog([s])
    # Gather data
    cnt,snapcnt,curtime = 0, 0, None
    for L in locs:
      aa = aas[L]
      # Read each file
      for filename in locs[L]:
        sys.stdout.write('.'); sys.stdout.flush()
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        # Read all data from each file
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if curtime != t:
                curtime = t
                cnt = (cnt + 1) % opts.decimate
                if cnt == 0:
                    aa.set_jultime(t)
                    # Make snapshot images (if specified)
                    if opts.snap > 0:
                        snapcnt = (snapcnt + 1) % opts.snap
                        if snapcnt == 0:
                            try:
                                grid_it(im,us,vs,ws,ds,wgts)
                                uvs,bms,dim,dbm = img_it(im)
                            except(ValueError):
                                uvs = n.abs(im.uv)
                                bms,dim,dbm = uvs,uvs,uvs
                            for k in ['uvs','bms','dim','dbm']:
                                if k in outputs: to_fits(k, eval(k), s,imgcnt)
                            imgcnt += 1
                            us,vs,ws,ds,wgts = [],[],[],[],[]
                            im = Img(opts.size, opts.res, mf_order=0)
                            if opts.src == 'zen':
                                s = a.ant.RadioFixedBody(aa.sidereal_time(), 
                                    aa.lat, name='zen')
                                src = a.fit.SrcCatalog([s])
                    src.compute(aa)
                    s_eq = src.get_crds('eq', ncrd=3)
                    aa.sim_cache(s_eq)
            if cnt != 0: continue
            d,f = d.take(chans), f.take(chans)
            if not opts.skip_amp: d /= aa.passband(i,j)
            try:
                # Throws PointingError if not up:
                d = aa.phs2src(d, s, i, j)
                u,v,w = aa.gen_uvw(i,j,src=s)
                if not opts.skip_bm:
                    # Calculate beam strength for weighting purposes
                    wgt = aa.bm_response(i,j,pol=opts.pol).squeeze()
                    # Optimal SNR: down-weight beam-attenuated data 
                    # by another factor of the beam response.
                    d *= wgt; wgt *= wgt
                else: wgt = n.ones(d.shape, dtype=n.float)
            except(a.ant.PointingError): continue
            valid = n.logical_not(f)
            d = d.compress(valid)
            if len(d) == 0: continue
            ds.append(d)
            us.append(u.compress(valid))
            vs.append(v.compress(valid))
            ws.append(w.compress(valid))
            wgts.append(wgt.compress(valid))
            # If data buffer is full, grid data
            if len(ds) * len(chans) > opts.buf_thresh:
                grid_it(im,us,vs,ws,ds,wgts)
                us,vs,ws,ds,wgts = [],[],[],[],[]

    # Grid remaining data into UV matrix
    grid_it(im,us,vs,ws,ds,wgts)
    uvs,bms,dim,dbm = img_it(im)
    for k in ['uvs','bms','dim','dbm']:
        if k in outputs: to_fits(k, eval(k), s, imgcnt)
    imgcnt += 1
    us,vs,ws,ds,wgts = [],[],[],[],[]
    im = Img(opts.size, opts.res, mf_order=0)
    
