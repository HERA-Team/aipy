#! /usr/bin/env python

import numpy as n, aipy as a, pylab as p, random, sys
from matplotlib.toolkits.basemap import Basemap

random.seed(1)

def pack_sphere(N):
    dz = 2. / N
    z = n.arange(-1+dz/2,1, dz)
    r = n.sqrt(1-z**2)
    dL = n.pi * (3 - n.sqrt(5))
    long = n.arange(0, dL * N, dL)
    return n.array([r*n.cos(long), r*n.sin(long), z])

def bit_reverse(N, nbits=None):
    if nbits is None: nbits = int(n.floor(n.log2(N))) + 1
    ans = 0
    for bit in range(nbits):
        ans += n.bitwise_and(N, 2**bit) * 2**(nbits-2*bit-1)
    return ans

def bit_reverse_order(N):
    nbits = int(n.floor(n.log2(N))) + 1
    indices = bit_reverse(n.arange(2**nbits), nbits=nbits)
    return indices.compress(indices < N)

def local_shuffle(L, width=2):
    for i in range(int(n.ceil(len(L) / float(width)))):
        chunk = L[width*i:width*(i+1)]
        random.shuffle(chunk)
        L[width*i:width*(i+1)] = chunk
        
pnts = pack_sphere(int(sys.argv[-1]))
ra,dec = a.coord.eq2radec(pnts)
ind1 = n.arange(len(ra))
local_shuffle(ind1)
ind2 = bit_reverse_order(len(ra))
indices = ind1.take(ind2)
ra = ra.take(indices)
dec = dec.take(indices)

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
