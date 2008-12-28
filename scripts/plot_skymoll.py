#! /usr/bin/env python
"""Taking a skymap file, makes a movie of the Boolardy sky in 1hr steps."""
import pylab as p
import numpy, aipy, ephem, sys, os
from matplotlib.toolkits.basemap import Basemap

skymap = aipy.img.SkyMap(fromfile=sys.argv[-1])
a = skymap.get_map().transpose()
b = numpy.log10(a + 1e-15) - 1
lats, lons = numpy.indices(a.shape)
lats = aipy.img.degrees(lats *.01) - 90
lons = aipy.img.degrees(lons *.01) - 180

cat = aipy.src.get_catalog(type='ant')
o = ephem.Observer()
o.date = aipy.ant.juldate2ephem(2454303)
cat.compute(o)
# lat/lon coordinates of sources
slats = numpy.array(map(lambda s: float(s.dec), cat.values()))
slons = numpy.array(map(lambda s: float(s.ra), cat.values()))
slats = aipy.img.degrees(slats)
slons = aipy.img.degrees(slons) - 180
snams = cat.keys()

map = Basemap(projection='moll',lat_0=0,lon_0=0)
map.drawmapboundary()
map.drawmeridians(numpy.arange(0,360,30))
map.drawparallels(numpy.arange(-90,90,30))
x, y = map(lons, lats)
b = b.clip(b.max()-5, b.max())
levels = numpy.arange(b.min()-.1, b.max()+.1, (b.max()-b.min())/15.)
CS = map.contourf(x, y, b, levels, linewidths=0)
sx, sy = map(slons,slats)
for name, xpt, ypt in zip(snams, sx, sy):
    if xpt >= 1e30 or ypt >= 1e30: continue
    #map.plot(sx, sy, 'ko', markerfacecolor=None)
    p.text(xpt+50000, ypt+50000, name)
p.show()
