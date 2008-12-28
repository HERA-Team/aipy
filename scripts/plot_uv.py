#! /usr/bin/env python
"""
Creates waterfall plots from Miriad UV files.  Can tile multiple plots
on one window, or plot just a single baseline.

Author: Aaron Parsons
Date: 07/05/07
Revisions: None
"""

import aipy.miriad, sys, pylab, numpy, math
from optparse import OptionParser

p = OptionParser()
p.set_usage('plot_uv.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-p', '--phase', dest='phase', action='store_true',
    help='Plot phase.')
p.add_option('-l', '--linear', dest='linear', action='store_true',
    help='Plot on linear scale.')
p.add_option('-a', '--autos', dest='autos', action='store_true',
    help='Plot only auto-correlations.')
p.add_option('-c', '--cross', dest='cross', action='store_true',
    help='Plot only cross-correlations.')
p.add_option('-i', '--ant_i', dest='ant_i', default=-1, type='int',
    help='Only plot baselines including this antenna.')
p.add_option('-j', '--ant_j', dest='ant_j', default=-1, type='int',
    help='Only plot baselines which include both ant_i and this antenna.')
p.add_option('-s', '--save_file', dest='save_file', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
p.add_option('-k', '--stokes', dest='stokes', default=1, type='int',
    help='Choose which stokes parameter (1, 2, 3, 4 = xx, yy, xy, yx).')
p.add_option('-d', '--delay', dest='delay', action='store_true',
    help='Take FFT of frequency axis to go to delay (t) space.')
p.add_option('-f', '--fringe', dest='fringe', action='store_true',
    help='Take FFT of time axis to go to fringe (Hz) space.')

opts, args = p.parse_args(sys.argv[1:])
active_pol = -4 - opts.stokes

# Loop through UV files collecting relevant data
data = {}
for uvfile in args:
    print 'Reading', uvfile
    uv = aipy.miriad.UV(uvfile)
    # Read data from a single UV file
    while True:
        p, d = uv.read_data()
        if d.size == 0: break
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        # Only save data that is needed to plot
        if opts.ant_i != -1:
            if opts.ant_j != -1:
                if opts.ant_i != i or opts.ant_j != j: continue
            if opts.ant_i != i and opts.ant_i != j: continue
        if opts.autos and i != j: continue
        if opts.cross and i == j: continue
        if uv['pol'] != active_pol: continue
        if not data.has_key(bl): data[bl] = []
        d.shape = (1,) + d.shape
        data[bl].append(d)
    del(uv)

ks = data.keys()
ks.sort()
m2 = int(math.sqrt(len(ks)))
m1 = int(math.ceil(float(len(ks)) / m2))
# Generate all the plots
for n, k in enumerate(ks):
    d = numpy.ma.concatenate(data[k], axis=0)
    if opts.fringe:
        d = numpy.ma.array(numpy.fft.ifft(d.filled(0), axis=0))
        #d[0] = (d[1] + .5*d[2] + .5*d[-2] + d[-1]) / 3
        d = numpy.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], axis=0)
    if opts.delay:
        d = numpy.ma.array(numpy.fft.ifft(d.filled(0), axis=1))
        d = numpy.ma.concatenate([d[:,d.shape[1]/2:], d[:,:d.shape[1]/2]],
            axis=1)
    d = d[:,1:-1] - (d[:,:-2] + d[:,2:]) / 2
    if opts.phase: d = numpy.angle(d.filled(0))
    elif opts.linear: d = numpy.abs(d)
    else: d = numpy.ma.log10(numpy.abs(d)+1e-6)
    pylab.subplot(m1, m2, n+1)
    pylab.imshow(d, aspect='auto')
    pylab.title(str(aipy.miriad.bl2ij(k)))
# Save to a file or pop up a window
if opts.save_file != '': pylab.savefig(opts.save_file)
else: pylab.show()
