#! /usr/bin/env python
"""
This is a general-purpose script for deconvolving dirty images by a 
corresponding PSF to produce a clean image.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, optparse, ephem, os

o = optparse.OptionParser()
o.set_usage('cl_img.py [options] *.dim.fits *.dbm.fits')
o.set_description(__doc__)
o.add_option('-d', '--deconv', dest='deconv', default='mem',
    help='Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).')
o.add_option('-o', '--output', dest='output', default='bim',
    help='Comma delimited list of data to generate FITS files for.  Can be: cim (clean image), rim (residual image), or bim (best image = clean + residuals). Default is bim.')
o.add_option('--var', dest='var', type='float', default=.6,
    help='Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.')
o.add_option('--tol', dest='tol', type='float', default=1e-6,
    help='Tolerance for successful deconvolution.  For annealing, interpreted as cooling speed.')
o.add_option('--taper', dest='taper', type='float', default=n.Inf,
    help='Width of tapering window (in pixels) to apply to both dirty beam and dirty image before deconvolving.')
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
    dim,dbm = dim.squeeze(), dbm.squeeze()
    DIM = dim.shape[0]
    if opts.taper < n.Inf:
        taper = a.img.gaussian_beam(opts.taper,
            shape=dim.shape,center=(DIM/2,DIM/2))
        taper = 2*taper.clip(0,.5)
        dbm *= taper
        #dim *= taper
    #dbm = a.img.recenter(dbm, (DIM/2+1,DIM/2))
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
            maxiter=opts.maxiter, verbose=True, tol=opts.tol)
    elif opts.deconv == 'ann':
        cim,info = a.deconv.anneal(dim, dbm, maxiter=opts.maxiter, 
            cooling=lambda i,x: opts.tol*(1-n.cos(i/50.))*(x**2), verbose=True)
    else:
        cim,info = n.zeros_like(dim), {'res':dim}
    rim = info['res']
    bim = rim / bm_gain + cim
    for ftag in ['cim','rim','bim']:
        if ftag in outputs: to_fits(k, ftag, eval(ftag), kwds)
    

