#! /usr/bin/env python

import numpy as n, aipy as a, pylab as p, random, sys
try: from mpl_toolkits.basemap import Basemap
except(ImportError): from matplotlib.toolkits.basemap import Basemap

NPTS = 100
random.seed(1)

ra,dec = a.map.facet_centers(NPTS, ncrd=2)

slats = dec * a.img.rad2deg
slons = ra * a.img.rad2deg
slons -= 180
slons = n.where(slons < -180, slons +360, slons)

map = Basemap(projection='moll', lat_0=0, lon_0=0)
map.drawmapboundary()
map.drawmeridians(n.arange(0,360,30))
map.drawparallels(n.arange(-90,90,30))
sx, sy = map(slons, slats)
map.plot(sx, sy, 'ko')
for c, (xpt,ypt) in enumerate(zip(sx, sy)):
    p.text(xpt+50000, ypt+50000, str(c))
p.show()
