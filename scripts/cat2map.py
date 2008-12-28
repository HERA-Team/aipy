#! /usr/bin/env python
import aipy as a, ephem as e, numpy as n, sys, optparse

o = optparse.OptionParser()
o.set_usage('cat2map.py [options]')
o.set_description(__doc__)
o.add_option('-s', '--srcs', dest='srcs',
    help='The sources to add to the map.')
o.add_option('-n', '--nside', dest='nside', type='int', default=16,
    help='The NSIDE parameter for the Healpix Map.  Default is 16.')
o.add_option('-m', '--map', dest='map', default='out.fits',
    help='The location to save the map.')
o.add_option('-i', '--interp', dest='interp', action='store_true',
    help='Interpolate between pixels when placing sources.')
opts,args = o.parse_args(sys.argv[1:])

h = a.healpix.HealpixMap(nside=opts.nside)
h.map = n.zeros(h.Npix(), dtype=n.float64)
h.use_interp = opts.interp

for src in opts.srcs.split(','):
    print 'Adding', src
    ra,dec,flux = a.src.src_data[src][:3]
    crd = e.Equatorial(ra, dec, epoch=e.J2000)
    crd = a.coord.radec2eq(crd.get())
    crd.shape = (1,3)
    h[crd] = n.array([flux])

print 'Saving to', opts.map
h.to_fits(opts.map)

