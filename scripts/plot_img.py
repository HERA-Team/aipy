#! /usr/bin/env python
"""
This is a general-purpose script for plotting simple FITS images.
"""

import aipy as a, sys, optparse, os
import numpy as np, ephem, math
from matplotlib import pylab as p

o = optparse.OptionParser()
o.set_usage('plot_img.py [options] *.fits')
o.set_description(__doc__)
a.scripting.add_standard_options(o, chan=True, cmap=True, max=True, drng=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('-o', '--outfile', dest='outfile', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
o.add_option('-p', '--pol', dest='pol', type='int', default=0, 
    help='Polarization index if FITS file has multiple polarizations.  Default 0.')
o.add_option('--batch', dest='batch', action='store_true',
    help='Process files in batch mode (one plot each) and output to a <input file>.png file')
o.add_option('--nogrid', dest='nogrid', action='store_true',
    help='Do not display RA/DEC grid.')
o.add_option('-f', '--fft', dest='fft', action='store_true',
    help='Perform 2D FFT of image.')
opts, args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap)
if opts.batch: m1,m2 = 1,1
else:
    m2 = int(math.sqrt(len(args)))
    m1 = int(math.ceil(float(len(args)) / m2))

for cnt, filename in enumerate(args):
    print filename
    if opts.batch:
        cnt = 0
        outfile = filename+'.png'
        if os.path.exists(outfile):
            print 'Output file exists... skipping.'
            continue
    # Gather data
    d, kwds = a.img.from_fits(filename)
    print d.shape
    print kwds
    print '-----------------------------------------------------------'

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
    for ax in compress_axes: d = np.average(d, axis=ax)

    # Put array in (ra,dec) order for plotting
    d = d.transpose((ra_ax,dec_ax))

    # Generate plots
    if opts.fft:
        d = np.fft.fft2(d)
        d = a.img.recenter(d, (d.shape[0]/2, d.shape[1]/2))
    if opts.mode.startswith('phs'): d = np.angle(d)
    elif opts.mode.startswith('lin'): d = np.absolute(d)
    elif opts.mode.startswith('real'): d = d.real
    elif opts.mode.startswith('imag'): d = d.imag
    elif opts.mode.startswith('log'):
        d = np.ma.absolute(d)
        d = np.ma.masked_less_equal(d, 0)
        d = np.ma.log10(d)

    if not opts.max is None: max = opts.max
    else: max = d.max()
    if not opts.drng is None: min = max - opts.drng
    else: min = d.min()

    p.subplot(m2, m1, cnt+1)
    if not opts.nogrid:
        from mpl_toolkits.basemap import Basemap
        xpx,ypx = d.shape
        dx1 = -(xpx/2 + .5) * kwds['d_ra'] * a.img.deg2rad
        dx2 = (xpx/2 - .5) * kwds['d_ra'] * a.img.deg2rad
        dy1 = -(ypx/2 + .5) * kwds['d_dec'] * a.img.deg2rad
        dy2 = (ypx/2 - .5) * kwds['d_dec'] * a.img.deg2rad
        map = Basemap(projection='ortho', lon_0=180, lat_0=kwds['dec'],
            rsphere=1, llcrnrx=dx1, llcrnry=dy1, urcrnrx=dx2,urcrnry=dy2)
        map.drawmeridians(np.arange(kwds['ra']-180,kwds['ra']+180,30))
        map.drawparallels(np.arange(-90,120,30))
        map.drawmapboundary()
        map.imshow(d, vmin=min, vmax=max, cmap=cmap, interpolation='nearest')
    else: p.imshow(d, vmin=min, vmax=max, origin='lower', cmap=cmap, interpolation='nearest')
    p.colorbar(shrink=.5, fraction=.05)
    p.title(filename)

    if opts.batch:
        print 'Saving to', outfile
        p.savefig(outfile)
        p.clf()
        

# Add right-click functionality for finding locations/strengths in map.
cnt = 1
def click(event):
    global cnt
    if not event.button in [2,3]: return
    if not opts.nogrid:
        lon,lat = map(event.xdata, event.ydata, inverse=True)
        lon = (180 + kwds['ra'] - lon) % 360
        lon *= a.img.deg2rad; lat *= a.img.deg2rad
        ra,dec = ephem.hours(lon), ephem.degrees(lat)
        xpx = np.around((event.xdata-1-dx1) / (dx2 - dx1) * d.shape[0] - .5)
        ypx = np.around((event.ydata-1-dy1) / (dy2 - dy1) * d.shape[1] - .5)
        flx = d[ypx,xpx]
        if opts.mode.startswith('log'): flx = 10**flx
        print '#%d (RA,DEC): (%s, %s), PX: (%d,%d) Jy: %f' % (cnt, ra, dec, xpx, ypx, flx)
    else:
        xpx = np.around(event.xdata)
        ypx = np.around(event.ydata)
        flx = d[ypx,xpx]
        if opts.mode.startswith('log'): flx = 10**flx
        print '#%d PX: (%d,%d) Jy: %f' % (cnt, xpx, ypx, flx)
    cnt += 1

#register this function with the event handler
p.connect('button_press_event', click)

if not opts.batch:
    if opts.outfile != '':
        print 'Saving to', opts.outfile
        p.savefig(opts.outfile)
    else: p.show()
