#! /usr/bin/env python
"""
This is a general-purpose script for making faceted, spherical maps (stored in
Healpix FITS format) from individual "flat" maps stored in FITS files.

Author: Aaron Parsons
"""

import sys, numpy as n, os, aipy as a, optparse

o = optparse.OptionParser()
o.set_usage('mk_map.py [options] *.fits')
o.set_description(__doc__)
a.scripting.add_standard_options(o, pol=True, chan=True)
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data to healpix map.')
o.add_option('-m', '--map', dest='map',
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')
o.add_option('--nside', dest='nside', type='int', default=256,
    help='NSIDE parameter for map, if creating a new file.')
o.add_option('--fwidth', dest='fwidth', type='float', default=10,
    help='Width of gaussian, in degrees, used to downweight data away from pointing center for each facet.  Default 10.')
opts, args = o.parse_args(sys.argv[1:])

# Open skymap
if os.path.exists(opts.map): skymap = a.map.Map(fromfits=opts.map)
else: skymap = a.map.Map(nside=opts.nside)
skymap.set_interpol(opts.interpolate)

prev_dra = None
for i, filename in enumerate(args):
    img, kwds = a.img.from_fits(filename)
    img = img.squeeze()
    s = a.phs.RadioFixedBody(kwds['ra']*a.img.deg2rad,
                             kwds['dec']*a.img.deg2rad)
    print '-----------------------------------------------------------'
    print 'Reading file %s (%d / %d)' % (filename, i + 1, len(args))
    print kwds
    print 'Pointing (ra, dec):', s._ra, s._dec
    print 'Image Power:', n.abs(img).sum()
    if prev_dra != kwds['d_ra']:
        prev_dra = kwds['d_ra']
        DIM = img.shape[0]
        RES = 1. / (kwds['d_ra'] * a.img.deg2rad * DIM)
        im = a.img.Img(DIM*RES, RES, mf_order=0)
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        # Define a weighting for gridding data into the skymap
        map_wgts = n.exp(-(tx**2+ty**2) / n.sin(opts.fwidth*a.img.deg2rad)**2)
        map_wgts.shape = (map_wgts.size,)
        valid = n.logical_not(map_wgts.mask)
        map_wgts = map_wgts.compress(valid)
    # Get coordinates of image pixels in original (J2000) epoch
    ex,ey,ez = im.get_eq(s._ra, s._dec, center=(DIM/2,DIM/2))
    ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
    img = img.flatten(); img = img.compress(valid)
    # Put the data into the skymap
    skymap.add((ex,ey,ez), map_wgts, img)
skymap.to_fits(opts.map, clobber=True)

