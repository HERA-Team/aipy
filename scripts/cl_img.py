#! /usr/bin/env python
"""
This is a general-purpose script for deconvolving dirty images by a 
corresponding PSF to produce a clean image.
"""

import aipy as a, numpy as n, sys, optparse, ephem, os

o = optparse.OptionParser()
o.set_usage('cl_img.py [options] *.dim.fits *.dbm.fits')
o.set_description(__doc__)
o.add_option('-d', '--deconv', dest='deconv', default='cln',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).')
o.add_option('-o', '--output', dest='output', default='bim',
    help='Comma delimited list of data to generate FITS files for.  Can be: cim (clean image), rim (residual image), or bim (best image = clean + residuals). Default is bim.')
o.add_option('--var', dest='var', type='float', default=.6,
    help='Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.  For annealing, interpreted as cooling speed.')
o.add_option('--div', dest='div', action='store_true',
    help='Allow clean to diverge (i.e. allow residual score to increase)')
o.add_option('--minuv', dest='minuv', type='float', default=0.,
    help='Minimum uv length to include (post-gridding)')
o.add_option('--maxuv', dest='maxuv', type='float', default=-1,
    help='Maximum uv length to include (post-gridding)')
o.add_option('-r', '--rewgt', dest='rewgt',  default='natural',
    help='Reweighting to apply to dim/dbm data before cleaning.  Options are: natural, uniform(LEVEL), radial, or midrange, where LEVEL is the fractional cutoff for using uniform weighting (recommended range .01 to .1).  Default is natural')
o.add_option('--maxiter', dest='maxiter', type='int', default=200,
    help='Number of allowable iterations per deconvolve attempt.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
outputs = opts.output.split(',')
keys = {}
for arg in args: keys[arg[:-len('.dim.fits')]] = None
keys = keys.keys()
keys.sort()

# Define a quick function writing an image to a FITS file
def to_fits(prefix,ftag,i,kwds):
    filename = '%s.%s.fits' % (prefix, ftag)
    print 'Saving data to', filename
    a.img.to_fits(filename, i.astype(n.float32), clobber=True, **kwds)

# Loop through all specified files
for cnt, k in enumerate(keys):
    fdim,fdbm = k+'.dim.fits', k+'.dbm.fits'
    print '-----------------------------------------------------------'
    print '%d / %d' % (cnt + 1, len(keys))
    print 'Deconvolving %s by %s' % (fdim, fdbm)
    if n.all([os.path.exists('%s.%s.fits' % (k,ftag)) for ftag in outputs]):
        print 'All output files exist already... skipping.'
        continue
    dim, kwds = a.img.from_fits(fdim)
    dbm, kwds = a.img.from_fits(fdbm)
    dim = dim.astype(n.float32)
    dbm = dbm.astype(n.float32)
    if n.all(dbm == 0):
        print 'No data in image, so skipping.'
        for ftag in ['cim','rim','bim']:
            if ftag in outputs: to_fits(k, ftag, dbm, kwds)
        continue
    print kwds
    size = 1/(kwds['d_ra']*a.img.deg2rad)
    res = size/dim.shape[0]
    im = a.img.Img(size=size, res=res)
    u,v = im.get_uv()
    r = n.sqrt(u**2+v**2)
    dim,dbm = dim.squeeze(), dbm.squeeze()
    DIM = dim.shape[0]
    uvs,bms = n.fft.fft2(dim), n.fft.fft2(dbm)
    if opts.rewgt.startswith('natural'): pass
    else:
        if opts.rewgt.startswith('uniform'): 
            level = float(opts.rewgt.split('(')[-1][:-1])
            abms = n.abs(bms)
            thresh = abms.max() * level
            divisor = abms.clip(thresh, n.Inf)
            uvs /= divisor; bms /= divisor
        elif opts.rewgt.startswith('radial'):
            #x,y = n.indices(dim.shape)
            #x = a.img.recenter(x - DIM/2, (DIM/2,DIM/2))
            #y = a.img.recenter(y - DIM/2, (DIM/2,DIM/2))
            #r = n.sqrt(x**2 + y**2)
            uvs *= r; bms *= r
        elif opts.rewgt.startswith('midrange'):
            x,y = n.indices(dim.shape)
            x = a.img.recenter(x - DIM/2, (DIM/2,DIM/2))
            y = a.img.recenter(y - DIM/2, (DIM/2,DIM/2))
            wgt = x**2 + y**2
            wgt = wgt.astype(n.float32) / (DIM/2)**2
            wgt = wgt * (1-wgt)
            wgt = n.where(wgt < 0, 0, wgt)
            uvs *= wgt; bms *= wgt
        else: raise ValueError('Unrecognized rewgt: %s' % opts.rewgt)
    #mask = n.where(r < opts.minuv, 0, 1)
    #mask = n.where(r > opts.maxuv, 0, mask)
    mask = 1
    if opts.minuv > 0: mask *= 1 - n.exp(-r**2/opts.minuv**2)
    if opts.maxuv > 0: mask *= n.exp(-r**2/opts.maxuv**2)
    dim = n.fft.ifft2(uvs * mask).real
    dbm = n.fft.ifft2(bms * mask).real
    
    dbm = a.img.recenter(dbm, (DIM/2,DIM/2))
    bm_gain = a.img.beam_gain(dbm)
    print 'Gain of dirty beam:', bm_gain
    if opts.deconv == 'mem':
        cim,info = a.deconv.maxent_findvar(dim, dbm, f_var0=opts.var,
            maxiter=opts.maxiter, verbose=True, tol=opts.tol, maxiterok=True)
    elif opts.deconv == 'lsq':
        cim,info = a.deconv.lsq(dim, dbm, 
            maxiter=opts.maxiter, verbose=True, tol=opts.tol)
    elif opts.deconv == 'cln':
        cim,info = a.deconv.clean(dim, dbm, 
            maxiter=opts.maxiter, stop_if_div=not opts.div, 
            verbose=True, tol=opts.tol)
    elif opts.deconv == 'ann':
        cim,info = a.deconv.anneal(dim, dbm, maxiter=opts.maxiter, 
            cooling=lambda i,x: opts.tol*(1-n.cos(i/50.))*(x**2), verbose=True)
    else:
        cim,info = n.zeros_like(dim), {'res':dim}
    rim = info['res']
    bim = rim / bm_gain + cim
    for ftag in ['cim','rim','bim']:
        if ftag in outputs: to_fits(k, ftag, eval(ftag), kwds)
    

