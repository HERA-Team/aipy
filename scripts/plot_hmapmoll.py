#! /usr/bin/env python
import pylab as p, numpy, aipy, sys, os, ephem, optparse
from matplotlib.toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_hmapmoll.py [options] mapfile')
o.set_description(__doc__)
o.add_option('-j', '--juldate', dest='juldate', type='float', default=2454489,
    help='Julian date used for locating moving sources (default 2454489).')
opts,args = o.parse_args(sys.argv[1:])

RES = .005
SZ = (int(numpy.pi/RES), int(2*numpy.pi/RES))
skymap = aipy.img.SkyHMap(128, fromfits=args[0])
lats, lons = numpy.indices(SZ)
lats = lats.astype(numpy.float) * RES
lons = lons.astype(numpy.float) * RES
crd = numpy.array([lats.flatten(), lons.flatten()]).transpose()
data = skymap[crd]
lats = 90 - crd[:,0] * aipy.img.rad2deg
lons = crd[:,1] * aipy.img.rad2deg - 180

cat = aipy.src.get_catalog(type='ant')
o = ephem.Observer()
o.date = aipy.ant.juldate2ephem(opts.juldate)
cat.compute(o)
# lat/lon coordinates of sources
slats = numpy.array(map(lambda s: float(s.dec), cat.values()))
slons = numpy.array(map(lambda s: float(s.ra), cat.values()))
slats = aipy.img.rad2deg*(slats)
slons = aipy.img.rad2deg*(slons) - 180
snams = cat.keys()

map = Basemap(projection='moll',lat_0=0,lon_0=0)
#map = Basemap(projection='moll',lat_0=0,lon_0=180)
map.drawmapboundary()
map.drawmeridians(numpy.arange(0,360,30))
map.drawparallels(numpy.arange(-90,90,30))
x, y = map(lons, lats)
data = numpy.log10(data)
data = data.clip(data.max()-5, data.max())
levels = numpy.arange(data.min()-.1, data.max()+.5, (data.max()-data.min())/15.)
data.shape = SZ
x.shape = SZ
y.shape = SZ
CS = map.contourf(x, y, data, levels, linewidths=0)
sx, sy = map(slons,slats)
for name, xpt, ypt in zip(snams, sx, sy):
    if xpt >= 1e30 or ypt >= 1e30: continue
    #map.plot(sx, sy, 'ko', markerfacecolor=None)
    p.text(xpt+50000, ypt+50000, name)
p.show()
