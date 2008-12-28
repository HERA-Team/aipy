#! /usr/bin/env python
import aipy as a, ephem as e, numpy as n, sys, optparse, os

o = optparse.OptionParser()
o.set_usage('modmap.py [options]')
o.set_description(__doc__)
o.add_option('-s', '--srcs', dest='srcs', default='',
    help='The sources to add to the map.')
o.add_option('-n', '--nside', dest='nside', type='int', default=32,
    help='The NSIDE parameter for the Healpix Map.  Default is 32.')
o.add_option('-m', '--map', dest='map', default='out.fits',
    help='The location to save the map.')
o.add_option('-i', '--in_map', dest='in_map', 
    help='Existing map to add sources to.  Treated as temperature, i.e. if resolution is degraded, pixels are averaged.  Use mscale to re-scale if you want Jy per pixel.')
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Frequency to use for source strengths.  Default .150 GHz.')
o.add_option('--nindices', dest='nindices', type='int',default=1,
    help='Number of terms in the spectral index polynomial.')
o.add_option('--mscale', dest='mscale', type='float',default=1,
    help='Scaling to apply to data in map.  Useful for changing T into Jy per pixel, and reverse.')
o.add_option('--sscale', dest='sscale', type='float',default=1,
    help='Scaling to apply to sources.  Useful for changing Jy to T.  It may be useful to know that the size of a map pixel is 4*pi / (12*nside**2).')
o.add_option('--iepoch', dest='iepoch', default=e.J2000,
    help='The epoch of coordinates in the input map.  Default J2000.')
o.add_option('--oepoch', dest='oepoch', default=e.J2000, 
    help='The epoch of coordinates in the output map.  Default J2000.')
o.add_option('--isys', dest='isys', default='eq',
    help='Coordinates (eq,ec,ga) of input map.  Default is eq')
o.add_option('--osys', dest='osys', default='eq',
    help='Coordinates (eq,ec,ga) of output map.  Default is eq')

opts,args = o.parse_args(sys.argv[1:])

print 'Starting a new map with NSIDE=%d' % opts.nside
h = a.map.Map(nside=opts.nside, nindices=opts.nindices)
if not opts.in_map is None:
    print 'Importing data from %s' % opts.in_map
    h.from_map(a.map.Map(fromfits=opts.in_map))
    h.map.map *= opts.mscale
    px = n.arange(h.npix())
    m = a.coord.convert_m(opts.osys, opts.isys,
        iepoch=opts.oepoch, oepoch=opts.iepoch)
    x,y,z = n.dot(m,h.px2crd(px, ncrd=3))
    h.set_interpol(True)
    h.map[px] = h.map[x,y,z]
    h.wgt[px] = h.wgt[x,y,z]
    for i in h.ind: i[px] = i[x,y,z]
    h.set_interpol(False)
h.reset_wgt()

m = a.coord.convert_m('eq', opts.osys, oepoch=opts.oepoch)
for src in opts.srcs.split(','):
    if src == '': continue
    ra,dec,flux,mfreq,index,size = a.src.src_data[src]
    eq = e.Equatorial(ra, dec, epoch=e.J2000)
    eq = a.coord.radec2eq(eq.get())
    x,y,z = n.dot(m, eq)
    strength = flux * (opts.freq / mfreq)**index * opts.sscale
    print 'Adding', src, 'with strength %f Jy' % strength,
    print 'and index %f' % index
    if opts.nindices == 0:
        curflux = h[x,y,z]
        h.put((x,y,z),1,strength+curflux)
    else:
        curflux, i = h[x,y,z]
        h.put((x,y,z),1,strength+curflux,[index])

print 'Saving to', opts.map
h.to_fits(opts.map)

