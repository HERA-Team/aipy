#! /usr/bin/env python
"""
This is a general-purpose script for plotting simple FITS images.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, pylab as p, sys, optparse, ephem, math

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
    try:
        ra_rng = (kwds['ra'] - d.shape[0] * kwds['d_ra']/2, 
            kwds['ra'] + d.shape[0] * kwds['d_ra']/2)
        assert(ra_rng[0] != ra_rng[1])
    except(KeyError,AssertionError): ra_rng = (-d.shape[0]/2,d.shape[0]/2)
    try:
        dec_rng = (kwds['dec'] - d.shape[1] * kwds['d_dec']/2, 
            kwds['dec'] + d.shape[1] * kwds['d_dec']/2)
        assert(dec_rng[0] != dec_rng[1])
    except(KeyError,AssertionError): dec_rng = (-d.shape[1]/2,d.shape[1]/2)

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
    p.imshow(d, vmin=min, vmax=max, 
        extent=(ra_rng[0],ra_rng[1],dec_rng[0],dec_rng[1]), aspect='auto')
    p.colorbar(shrink=.5, fraction=.05)
    p.title(filename)
p.show()
