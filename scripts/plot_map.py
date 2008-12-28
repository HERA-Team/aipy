#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse
from matplotlib.toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_map.py [options] mapfile')
o.set_description(__doc__)
o.add_option('-p', '--projection', dest='projection', default='moll',
    help='Map projection to use: moll (default), mill, cyl, robin.')
o.add_option('-j', '--juldate', dest='juldate', type='float', default=2454489,
    help='Julian date used for locating moving sources (default 2454489).')
o.add_option('-i', '--isys', dest='isys', default='eq',
    help='Input coordinate system (in map).')
o.add_option('-o', '--osys', dest='osys', default='eq',
    help='Output coordinate system (plotted).')
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plotting mode (can be log,lin).')
o.add_option('--iepoch', dest='iepoch', type='float', default=ephem.J2000,
    help='Epoch of input coordinates (in map).')
o.add_option('--oepoch', dest='oepoch', type='float', default=ephem.J2000,
    help='Epoch of output coordinates (plotted).')
o.add_option( '--max', dest='max', type='float', default=None,
    help='Manually set the maximum color level (log10).')
o.add_option('--dyn_rng', dest='dyn_rng', type='float', default=5.,
    help="Dynamic range in color of image (log10).")
o.add_option('--levels', dest='levels', type='int', default=15,
    help="Number of color levels to plot.")
o.add_option('--res', dest='res', type='float', default=.005,
    help="Resolution of plot (in radians).")
o.add_option('-c', '--cen', dest='cen', type='float', default=180.,
    help="Center longitude (in degrees) of map.")
o.add_option('--raw', dest='raw', action='store_true',
    help="Plot a raw map (no weighting, etc.)")
o.add_option('--no_srcs', dest='no_srcs', action='store_true',
    help="Omit labels of known radio sources in plot.")
opts,args = o.parse_args(sys.argv[1:])

CEN = (180. - opts.cen) * a.img.deg2rad
SZ = (int(n.pi/opts.res), int(2*n.pi/opts.res))
if not opts.raw: skymap = a.img.SkyHMap(32, fromfits=args[0])
else:
    skymap = a.healpix.HealpixMap()
    skymap.from_fits(args[0])
    print 'ORDERING:', skymap.Order()
    print 'NSIDE:', skymap.Nside()
lats, lons = n.indices(SZ)
lats = n.pi/2 - lats.astype(n.float) * opts.res
lons = lons.astype(n.float) * opts.res
# Convert lat/long to xyz (it says radec, but this is genereic lat/long)
get_lons = lons - CEN
crd = a.coord.radec2eq(n.array([get_lons.flatten(), lats.flatten()]))
m = a.coord.convert_m(opts.osys, opts.isys, 
    iepoch=opts.oepoch, oepoch=opts.iepoch)
crd = n.dot(m, crd)
data = skymap[crd.transpose()]
# Convert to degrees for the basemap module
lats *= a.img.rad2deg
lons = lons * a.img.rad2deg - 180

cat = a.src.get_catalog(type='sim')
o = ephem.Observer()
o.date = a.ant.juldate2ephem(opts.juldate)
cat.compute(o)
# lat/lon coordinates of sources
scrds = [ephem.Equatorial(s.ra,s.dec) for s in cat.values()]
if opts.osys.startswith('ga'): scrds = [ephem.Galactic(s) for s in scrds]
elif opts.osys.startswith('ec'): scrds = [ephem.Ecliptic(s) for s in scrds]
slats = n.array([float(s.get()[1]) for s in scrds]) * a.img.rad2deg
slons = n.array([float(s.get()[0] + CEN) for s in scrds]) * a.img.rad2deg
slons -= 180
slons = n.where(slons < -180, slons + 360, slons)
snams = cat.keys()

map = Basemap(projection=opts.projection,lat_0=0,lon_0=0)
map.drawmapboundary()
map.drawmeridians(n.arange(0,360,30))
map.drawparallels(n.arange(-90,90,30))
x, y = map(lons, lats)
if opts.mode.startswith('log'): data = n.log10(n.abs(data))
max = opts.max
if max is None: max = data.max()
data = data.clip(max-opts.dyn_rng, max)
levels = n.arange(data.min()-.1, data.max()+.1, 
    (data.max()-data.min())/opts.levels)
data.shape = SZ
x.shape = SZ
y.shape = SZ
CS = map.contourf(x, y, data, levels, linewidths=0)
sx, sy = map(slons,slats)
if not opts.no_srcs:
    for name, xpt, ypt in zip(snams, sx, sy):
        if xpt >= 1e30 or ypt >= 1e30: continue
        #map.plot(sx, sy, 'ko', markerfacecolor=None)
        p.text(xpt+50000, ypt+50000, name)
p.colorbar(shrink=.5, format='%.2f')
p.show()
