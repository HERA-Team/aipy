#! /usr/bin/env python
import pylab, numpy, sys, aipy.img

skymap = aipy.img.SkyMap(fromfile=sys.argv[-1])
m = skymap.get_map().transpose()
b = numpy.log10(m + 1e-15) - 1
mx = b.max()
#mx = 4.
ax = pylab.subplot(111)
pylab.xlabel('RA')
pylab.ylabel('DEC')
pylab.title('Multifrequency Sky Map, Boolardy')
ra1 = pylab.datestr2num('00:00:00')
ra2 = pylab.datestr2num('23:59:59')
dec1 = -90
dec2 = 90
im = pylab.imshow(b, extent=(ra1,ra2,dec1,dec2), vmin=mx-7, vmax=mx,
    aspect=.5*(ra2-ra1)/(dec2-dec1), origin='lower')

fmt = pylab.DateFormatter('%H:%M')
ax.xaxis.set_major_locator(pylab.HourLocator())
ax.xaxis.set_major_formatter(fmt)
pylab.setp(ax.get_xticklabels(), 
    'rotation', 45, 'horizontalalignment', 'right')
ax.yaxis.set_ticks(range(-90,90,10))
pylab.colorbar(shrink=.5, fraction=.05)
pylab.show()
