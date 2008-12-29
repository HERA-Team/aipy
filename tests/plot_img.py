#! /usr/bin/env python
"""
This is a general-purpose script for plotting simple FITS images.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, pylab as p, sys, optparse, ephem, math
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_img.py [options] *.fits')
o.set_description(__doc__)
a.scripting.add_standard_options(o, chan=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('-o', '--out_file', dest='out_file', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
o.add_option('-p', '--pol', dest='pol', type='int', default=0, 
    help='Polarization index if FITS file has multiple polarizations.  Default 0.')
o.add_option('--max', dest='max', default=None, type='float',
    help='Upper clip value on 2D plots.')
o.add_option('--dyn_rng', dest='dyn_rng', default=None, type='float',
    help='Dynamic range in scale of 2D plots.')
opts, args = o.parse_args(sys.argv[1:])

m2 = int(math.sqrt(len(args)))
m1 = int(math.ceil(float(len(args)) / m2))

for cnt, filename in enumerate(args):
    # Gather data
    d, kwds = a.img.from_fits(filename)
    print kwds

    # Parse command-line options
    compress_axes = []
    ra_ax,dec_ax = (0,1)
    try:
        for i,ax in enumerate(kwds['axes']):
            if ax.startswith('ra'): ra_ax = i
            elif ax.startswith('dec'): dec_ax = i
            elif ax.startswith('freq'):
                chans = a.scripting.parse_chans(opts.chan, d.shape[i])
                d = d.take(chans, axis=i)
                compress_axes.append(i)
            elif ax.startswith('stokes'):
                d = d.take([opts.pol], axis=i)
                compress_axes.append(i)
    except(KeyError): pass

    compress_axes.reverse()
    for ax in compress_axes: d = n.average(d, axis=ax)

    # Put array in (ra,dec) order for plotting
    d = d.transpose((ra_ax,dec_ax))

    # Generate plots
    if opts.mode.startswith('phs'): d = n.angle(d.filled(0))
    elif opts.mode.startswith('lin'): d = n.ma.absolute(d)
    elif opts.mode.startswith('real'): d = d.real
    elif opts.mode.startswith('imag'): d = d.imag
    elif opts.mode.startswith('log'):
        d = n.ma.absolute(d)
        d = n.ma.masked_less_equal(d, 0)
        d = n.ma.log10(d)

    if not opts.max is None: max = opts.max
    else: max = d.max()
    if not opts.dyn_rng is None: min = max - opts.dyn_rng
    else: min = d.min()

    p.subplot(m2, m1, cnt+1)
    xpx,ypx = d.shape
    dx1 = -(xpx/2 + .5) * kwds['d_ra'] * a.img.deg2rad
    dx2 = (xpx/2 - .5) * kwds['d_ra'] * a.img.deg2rad
    dy1 = -(ypx/2 + .5) * kwds['d_dec'] * a.img.deg2rad
    dy2 = (ypx/2 - .5) * kwds['d_dec'] * a.img.deg2rad
    map = Basemap(projection='ortho', lon_0=180, lat_0=kwds['dec'],
        rsphere=1, llcrnrx=dx1, llcrnry=dy1, urcrnrx=dx2,urcrnry=dy2)
    map.drawmeridians(n.arange(kwds['ra']-180,kwds['ra']+180,30))
    map.drawparallels(n.arange(-90,120,30))
    map.drawmapboundary()
    map.imshow(d, vmin=min, vmax=max)
    p.colorbar(shrink=.5, fraction=.05)
    p.title(filename)

# Add right-click functionality for finding locations/strengths in map.
cnt = 1
def click(event):
    global cnt
    if event.button != 3: return
    lon,lat = map(event.xdata, event.ydata, inverse=True)
    lon = (180 + kwds['ra'] - lon) % 360
    lon *= a.img.deg2rad; lat *= a.img.deg2rad
    ra,dec = ephem.hours(lon), ephem.degrees(lat)
    ypx = n.around(event.xdata / (kwds['d_ra'] * a.img.deg2rad))
    xpx = n.around(event.ydata / (kwds['d_dec'] * a.img.deg2rad))
    flx = d[xpx,ypx]
    if opts.mode.startswith('log'): flx = 10**flx
    print '#%d (RA,DEC): (%s, %s), PX: (%d,%d) Jy: %f' % \
        (cnt, ra, dec, xpx, ypx, flx)
    cnt += 1

#register this function with the event handler
p.connect('button_press_event', click)

p.show()
