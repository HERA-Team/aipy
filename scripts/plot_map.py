#! /usr/bin/env python
"""
Script for displaying a projection of a spherical (Healpix) data set stored
in a *.fits file.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse
try: from mpl_toolkits.basemap import Basemap
except(ImportError): from matplotlib.toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_map.py [options] mapfile')
o.set_description(__doc__)
o.add_option('-p', '--projection', dest='projection', default='moll',
    help='Map projection to use: moll (default), mill, cyl, robin.')
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plotting mode (can be log,lin).')
o.add_option('-c', '--cen', dest='cen', type='float', 
    help="Center longitude (in degrees) of map.")
o.add_option('-j', '--juldate', dest='juldate', type='float', 
    help='Julian date used for locating moving sources.')
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
o.add_option('--res', dest='res', type='float', default=.25,
    help="Resolution of plot (in degrees).")
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
o.add_option('--scaling', dest='scaling', type='float',default=1,
    help='Scaling of existing input map.')
opts,args = o.parse_args(sys.argv[1:])

if opts.cen is None:
    if opts.osys == 'eq': opts.cen = 180
    else: opts.cen = 0
map = Basemap(projection=opts.projection,lat_0=0,lon_0=opts.cen, rsphere=1.)
lons,lats,x,y = map.makegrid(360/opts.res,180/opts.res, returnxy=True)
lt = lats[:,0]
ln1 = n.ones_like(lt) * (lons[lons.shape[0]/2,0])
ln2 = n.ones_like(lt) * (lons[lons.shape[0]/2,-1])
x1,y1 = map(ln1,lt); x2,y2 = map(ln2,lt)
x = n.ma.array(x)
for c,(i,j) in enumerate(zip(x1,x2)): x[c] = n.ma.masked_outside(x[c], i, j)
mask = x.mask
if opts.osys == 'eq': lons = 360 - lons
lats *= a.img.deg2rad; lons *= a.img.deg2rad
print 'Reading %s' % args[0]
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()
if not opts.nside is None:
    nh = a.healpix.HealpixMap(nside=opts.nside)
    nh.from_hpm(h)
h.set_interpol(True)

crd = a.coord.radec2eq(n.array([lons.flatten(), lats.flatten()]))
m = a.coord.convert_m(opts.osys, opts.isys, 
    iepoch=opts.oepoch, oepoch=opts.iepoch)
x,y,z = n.dot(m, crd)
try: data, indices = h[x,y,z]
except(ValueError): data = h[x,y,z]
data *= opts.scaling
data.shape = lats.shape

if not opts.srcs is None:
    cat = a.src.get_catalog(cutoff=opts.srcs)
    o = ephem.Observer()
    if opts.juldate is None: del(cat['Sun'])
    else: o.date = a.ant.juldate2ephem(opts.juldate)
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
    slons = n.array([float(s.get()[0]) for s in scrds]) * a.img.rad2deg
    if opts.osys == 'eq': slons = 360 - slons
    slons = n.where(slons < -180, slons + 360, slons)
    slons = n.where(slons > 180, slons - 360, slons)

map.drawmapboundary()
map.drawmeridians(n.arange(0,360,30))
map.drawparallels(n.arange(-90,90,30)[1:], labels=[0,1,0,0], labelstyle='+/-')
if opts.mode.startswith('log'): data = n.log10(n.abs(data))
if opts.max is None: max = data.max()
else: max = opts.max
if opts.dyn_rng is None:
    min = data.min()
    if min < (max - 10): min = max-10
else: min = max - opts.dyn_rng
data = data.clip(min, max)
data = n.ma.array(data, mask=mask)
map.imshow(data, vmax=max, vmin=min)

if not opts.srcs is None:
    sx, sy = map(slons,slats)
    for name, xpt, ypt, flx in zip(snams, sx, sy, sflxs):
        if xpt >= 1e30 or ypt >= 1e30: continue
        #map.plot(sx, sy, 'ko', markerfacecolor=None)
        p.text(xpt+.001, ypt+.001, name, size=5+2*int(n.round(n.log10(flx))))
if not opts.nobar: p.colorbar(shrink=.5, format='%.2f')
else: p.subplots_adjust(.05,.05,.95,.95)

cnt = 1
def click(event):
    global cnt
    if event.button != 3: return
    lon,lat = map(event.xdata, event.ydata, inverse=True)
    if opts.osys == 'eq': lon = (360 - lon) % 360
    lon *= a.img.deg2rad; lat *= a.img.deg2rad
    ra,dec = ephem.hours(lon), ephem.degrees(lat)
    x,y,z = a.coord.radec2eq((ra,dec))
    flx = h[(x,y,z)]
    print '#%d (RA,DEC): (%s, %s), Jy: %f' % (cnt, ra, dec, flx)
    cnt += 1

#register this function with the event handler
p.connect('button_press_event', click)

p.show()
