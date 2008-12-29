#! /usr/bin/env python
"""
Script for displaying a projection of a spherical (Healpix) data set stored
in a *.fits file.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse
from matplotlib.toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_map.py [options] mapfile')
o.set_description(__doc__)
o.add_option('-p', '--projection', dest='projection', default='moll',
    help='Map projection to use: moll (default), mill, cyl, robin.')
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plotting mode (can be log,lin).')
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Interpolate between pixels.')
o.add_option('-c', '--cen', dest='cen', type='float', default=180.,
    help="Center longitude (in degrees) of map.")
o.add_option('-j', '--juldate', dest='juldate', type='float', default=2454489,
    help='Julian date used for locating moving sources (default 2454489).')
o.add_option('--srcs', dest='srcs', type='float',
    help="Cutoff flux for labeling known radio sources in plot.")
o.add_option('--isys', dest='isys', default='eq',
    help='Input coordinate system (in map).')
o.add_option('--osys', dest='osys', default='eq',
    help='Output coordinate system (plotted).')
o.add_option('--iepoch', dest='iepoch', type='float', default=ephem.J2000,
    help='Epoch of input coordinates (in map).')
o.add_option('--oepoch', dest='oepoch', type='float', default=ephem.J2000,
    help='Epoch of output coordinates (plotted).')
o.add_option( '--max', dest='max', type='float', default=None,
    help='Manually set the maximum color level (log10).')
o.add_option('--dyn_rng', dest='dyn_rng', type='float', default=None,
    help="Dynamic range in color of image (log10).")
o.add_option('--levels', dest='levels', type='int', default=15,
    help="Number of color levels to plot.")
o.add_option('--nobar', dest='nobar', action='store_true',
    help="Do not show colorbar.")
o.add_option('--res', dest='res', type='float', default=.005,
    help="Resolution of plot (in radians).")
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
o.add_option('--scaling', dest='scaling', type='float',default=1,
    help='Scaling of existing input map.')
opts,args = o.parse_args(sys.argv[1:])

CEN = (180. - opts.cen) * a.img.deg2rad
SZ = (int(n.pi/opts.res), int(2*n.pi/opts.res))
print 'Reading %s' % args[0]
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()
if not opts.nside is None:
    nh = a.healpix.HealpixMap(nside=opts.nside)
    crd = h.px2crd(n.arange(h.npix()))
    nh[crd] += h[crd]
    h = nh
h.set_interpol(opts.interpolate)

lats, lons = n.indices(SZ)
lats = n.pi/2 - lats.astype(n.float) * opts.res
lons = lons.astype(n.float) * opts.res
# Convert lat/long to xyz (it says radec, but this is genereic lat/long)
get_lons = lons - CEN
crd = a.coord.radec2eq(n.array([get_lons.flatten(), lats.flatten()]))
m = a.coord.convert_m(opts.osys, opts.isys, 
    iepoch=opts.oepoch, oepoch=opts.iepoch)
x,y,z = n.dot(m, crd)
try: data, indices = h[x,y,z]
except(ValueError): data = h[x,y,z]
data *= opts.scaling
# Convert to degrees for the basemap module
lats *= a.img.rad2deg
lons = lons * a.img.rad2deg - 180

o = ephem.Observer()
o.date = a.ant.juldate2ephem(opts.juldate)
if not opts.srcs is None:
    cat = a.src.get_catalog(cutoff=opts.srcs)
    cat.compute(o)
    # lat/lon coordinates of sources
    scrds = [ephem.Equatorial(s.ra,s.dec) for s in cat.values()]
    sflxs = cat.get_fluxes()
    snams = cat.keys()
if opts.osys == 'ga':
    scrds = [ephem.Galactic(s, epoch=opts.oepoch) for s in scrds]
elif opts.osys == 'ec':
    scrds = [ephem.Ecliptic(s, epoch=opts.oepoch) for s in scrds]
slats = n.array([float(s.get()[1]) for s in scrds]) * a.img.rad2deg
slons = n.array([float(s.get()[0] + CEN) for s in scrds]) * a.img.rad2deg
slons -= 180
slons = n.where(slons < -180, slons + 360, slons)
if opts.osys == 'eq': slons = 360 - slons

map = Basemap(projection=opts.projection,lat_0=0,lon_0=0)
map.drawmapboundary()
map.drawmeridians(n.arange(0,360,30))
map.drawparallels(n.arange(-90,90,30))
x, y = map(lons, lats)
if opts.mode.startswith('log'): data = n.log10(n.abs(data))
if opts.max is None: max = data.max()
else: max = opts.max
if opts.dyn_rng is None:
    min = data.min()
    if min < (max - 10): min = max-10
else: min = max - opts.dyn_rng
data = data.clip(min, max)
step = (max - min) / opts.levels
levels = n.arange(min-step, max+step, step)
data.shape = SZ
x.shape = SZ
y.shape = SZ
if opts.osys == 'eq': data = n.fliplr(data)
CS = map.contourf(x, y, data, levels, linewidths=0)
sx, sy = map(slons,slats)
if not opts.srcs is None:
    for name, xpt, ypt, flx in zip(snams, sx, sy, sflxs):
        if xpt >= 1e30 or ypt >= 1e30: continue
        #map.plot(sx, sy, 'ko', markerfacecolor=None)
        p.text(xpt+50000, ypt+50000, name, size=5+2*int(n.round(n.log10(flx))))
if not opts.nobar: p.colorbar(shrink=.5, format='%.2f')
p.show()
