#! /usr/bin/env python
"""
Creates waterfall plots from Miriad UV files.  Can tile multiple plots
on one window, or plot just a single baseline.

Author: Aaron Parsons
Date: 07/05/07
Revisions:
    08/15/07    arp Ported to use select_data, added decimate option, and
                    unmask option.
"""

import aipy, sys, pylab, numpy, math
from optparse import OptionParser

p = OptionParser()
p.set_usage('plot_uv.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
p.add_option('-a', '--ants', dest='ants', default='all',
    help='Select which antennas/baselines to include in plot.  Options are: "all", "auto", "cross", "<ant1 #>,<ant2 #>" (a specific baseline), or "+<ant1 #>,..." (a list of active antennas).')
p.add_option('-c', '--chan', dest='chan', default='all',
    help='Select which channels (taken after any delay/fringe transforms) to plot.  Options are: "all", "+<chan1 #>,..." (a list of active channels), or "<chan1 #>,<chan2 #>" (a range of channels).  If "all" or a range are selected, a 2-d image will be plotted.  If a list of channels is selected, an xy plot will be generated.')
p.add_option('-p', '--pol', dest='pol', default=None,
    help='Choose which polarization (xx, yy, xy, yx) to plot.')
p.add_option('-x', '--decimate', dest='decimate', default=1, type='int',
    help='Take every Nth time data point.')
p.add_option('-u', '--unmask', dest='unmask', action='store_true',
    help='Plot masked data, too.')
p.add_option('-s', '--sum_chan', dest='sum_chan', default=None,
    help='Sum all active channels together.  Options are "coh" (coherent summing) or "inc" (incoherent).')
p.add_option('-o', '--out_file', dest='out_file', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
p.add_option('-d', '--delay', dest='delay', action='store_true',
    help='Take FFT of frequency axis to go to delay (t) space.')
p.add_option('-f', '--fringe', dest='fringe', action='store_true',
    help='Take FFT of time axis to go to fringe (Hz) space.')
p.add_option('-t', '--lst', dest='lst', action='store_true',
    help='Choose time axis to be lst.')
p.add_option('', '--dt', dest='dt', action='store_true',
    help='Remove a linear extrapolation from adjacent times.')
p.add_option('', '--df', dest='df', action='store_true',
    help='Remove a linear extrapolation from adjacent frequency channels.')
p.add_option('', '--plot_min', dest='plot_min', default=None, type='int', 
    help='Lower clip value on log plots.')

def data_selector(antopt, uv):
    """Convert the command-line argument for ants into a function which, when
    passed i, j as arguments, will return whether this baseline is to be
    plotted."""
    if antopt == 'all': pass
    elif antopt == 'auto': uv.select('auto', 0, 0)
    elif antopt == 'cross': uv.select('auto', 0, 0, include=0)
    elif antopt.startswith('+'):
        antopt = eval('['+antopt[1:]+']')
        for a in antopt: uv.select('antennae', a, -1)
    else:
        x, y = eval('('+antopt+')')
        uv.select('antennae', x, y)

def gen_chan_extractor(chanopt, get_y=None):
    """Convert the command-line argument for chan into a function which takes
    data for all channels as an argument, and returns a subarray containing
    only the relevant channels.  If data is provided to 'get_y', this function
    returns the channel indices extracted (useful for plotting)."""
    if not get_y is None: get_y = numpy.arange(get_y.size)
    if chanopt == 'all':
        if not get_y is None: return get_y
        return lambda chan: chan
    if chanopt.startswith('+'):
        chans = eval('['+chanopt[1:]+']')
        if not get_y is None: return get_y.take(chans)
        return lambda chan: chan.take(chans)
    else:
        x, y = eval('('+chanopt+')')
        if x < 0 and y > 0:
            def ex(a): return numpy.concatenate([a[x:], a[:y]])
            if not get_y is None: return ex(get_y)
            return lambda chan: ex(chan)
        elif y == 0:
            if not get_y is None: return get_y[x:]
            return lambda chan: chan[x:]
        else:
            if not get_y is None: return get_y[x:y]
            return lambda chan: chan[x:y]

opts, args = p.parse_args(sys.argv[1:])
uv = aipy.miriad.UV(args[0])
if opts.pol is None: active_pol = uv['pol']
else: active_pol = aipy.miriad.str2pol[opts.pol]
chan_extractor = gen_chan_extractor(opts.chan)
p, d = uv.read()
chans = gen_chan_extractor(opts.chan, get_y=d)
if opts.chan.startswith('+'): plots = 'plot'
else: plots = 'imshow'
del(uv)

# Loop through UV files collecting relevant data
plot_x = {}
times = []
plot_times = []
cnt = 0
for uvfile in args:
    print 'Reading', uvfile
    uv = aipy.miriad.UV(uvfile)
    # Only select data that is needed to plot
    data_selector(opts.ants, uv)
    uv.select('polarization', active_pol, 0)
    # Read data from a single UV file
    for p, d in uv.all():
        uvw,t,(i,j) = p
        bl = '%d,%d' % (i,j)
        # Implement Decimation
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            cnt = (cnt + 1) % opts.decimate
            if cnt == 0:
                if opts.lst: plot_times.append(uv['lst'])
                else: plot_times.append(t)
        if cnt != 0: continue
        # Do delay transform if required
        if opts.delay:
            if opts.unmask: d = d.data
            else: d = d.filled(0)
            d = numpy.ma.array(numpy.fft.ifft(d))
            d = numpy.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], 
                axis=0)
        # Get data structure ready to accept this data
        if not plot_x.has_key(bl): plot_x[bl] = []
        # Extract specific channels for plotting
        d = chan_extractor(d)
        # Implement unmasking data
        if not opts.delay and opts.unmask: d = d.data
        d.shape = (1,) + d.shape
        plot_x[bl].append(d)
    del(uv)

plot_times = numpy.array(plot_times)
bls = plot_x.keys()
bls.sort()
if len(bls) == 0:
    print 'No data to plot.'
    sys.exit(0)
m2 = int(math.sqrt(len(bls)))
m1 = int(math.ceil(float(len(bls)) / m2))
# Generate all the plots
for n, bl in enumerate(bls):
    d = numpy.ma.concatenate(plot_x[bl], axis=0)
    if opts.df: d = d[:,:-2]/2 + d[:,2:]/2 - d[:,1:-1]
    if opts.dt: d = d[:-2]/2 + d[2:]/2 - d[1:-1]
    if opts.fringe:
        d = numpy.ma.array(numpy.fft.ifft(d.filled(0), axis=0))
        d = numpy.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], axis=0)
    if opts.sum_chan == 'coh': d = d.sum(axis=1)
    elif opts.sum_chan == 'inc': d = numpy.abs(d.filled(0)).sum(axis=1)
    if opts.mode.startswith('phs'): d = numpy.angle(d.filled(0))
    elif opts.mode.startswith('lin'): d = numpy.ma.absolute(d)
    elif opts.mode.startswith('real'): d = d.real
    elif opts.mode.startswith('imag'): d = d.imag
    else: d = numpy.ma.log10(numpy.ma.absolute(d)+1e-50)
    pylab.subplot(m1, m2, n+1)
    if plots == 'imshow' and opts.sum_chan is None:
        if opts.mode.startswith('log'):
            if not opts.plot_min is None:
                pylab.imshow(d, aspect='auto', vmin=opts.plot_min)
            else: pylab.imshow(d, aspect='auto')
        else:
            pylab.imshow(d, aspect='auto')
        pylab.colorbar()
    else:
        if not opts.sum_chan is None:
            pylab.plot(plot_times, d, '.', label='chan(+)')
        else:
            for c, chan in enumerate(chans):
                pylab.plot(plot_times, d[:,c], '.', label='chan%d' % chan)
    pylab.title(bl)
if plots == 'plot': pylab.legend(loc='best')
# Save to a file or pop up a window
if opts.out_file != '': pylab.savefig(opts.out_file)
else: pylab.show()
