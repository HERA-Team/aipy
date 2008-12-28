#! /usr/bin/env python
"""
Creates waterfall plots from Miriad UV files.  Can tile multiple plots
on one window, or plot just a single baseline.

Author: Aaron Parsons
Date: 07/05/07
Revisions:
    08/15/07    arp Ported to use select_data, added decimate option, and
                    unmask option.
    01/20/08    arp Added physical coordinates to plots
    02/19/08    gf  (Griffin Foster) Added the -t option to look at time slices
"""

import aipy as a, numpy as n, pylab as p, math, sys
from optparse import OptionParser

o = OptionParser()
o.set_usage('plot_uv.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('-a', '--ants', dest='ants', default='all',
    help='Select which antennas/baselines to include in plot.  Options are: "all", "auto", "cross", "<ant1 #>_<ant2 #>" (a specific baseline), or "<ant1 #>,..." (a list of active antennas).')
o.add_option('-p', '--pol', dest='pol', default=None,
    help='Choose which polarization (xx, yy, xy, yx) to plot.')
o.add_option('-c', '--chan', dest='chan', default='all',
    help='Select which channels (taken after any delay/fringe transforms) to plot.  Options are: "all", "<chan1 #>,..." (a list of active channels), or "<chan1 #>_<chan2 #>" (a range of channels).  If "all" or a range are selected, a 2-d image will be plotted.  If a list of channels is selected, an xy plot will be generated.')
o.add_option('--sum_chan', dest='sum_chan', action='store_true',
    help='Sum active channels together.')
o.add_option('-x', '--decimate', dest='decimate', default=1, type='int',
    help='Take every Nth time data point.')
o.add_option('-t', '--time', dest='time', default='all', help='Select which time sample to plot. Options are: "all" (default), "<time1 #>_<time2 #>" (a range of times to plot), or "<time1 #>,<time2 #>" (a list of times to plot). If "all" or a range are selected, a 2-d image will be plotted. If a list of times is selected an xy plot will be generated.')
o.add_option('-u', '--unmask', dest='unmask', action='store_true',
    help='Plot masked data, too.')
o.add_option('-d', '--delay', dest='delay', action='store_true',
    help='Take FFT of frequency axis to go to delay (t) space.')
o.add_option('-f', '--fringe', dest='fringe', action='store_true',
    help='Take FFT of time axis to go to fringe (Hz) space.')
o.add_option('--dt', dest='dt', action='store_true',
    help='Remove a linear extrapolation from adjacent times.')
o.add_option('--df', dest='df', action='store_true',
    help='Remove a linear extrapolation from adjacent frequency channels.')
o.add_option('-o', '--out_file', dest='out_file', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
o.add_option('--plot_min', dest='plot_min', default=None, type='int', 
    help='Lower clip value on log plots.')
o.add_option('--time_axis', dest='time_axis', default='index',
    help='Choose time axis to be integration/fringe index (index), or physical coordinates (physical), or if doing xy plot in time-mode, (lst) is also available.  Default is index.')
o.add_option('--chan_axis', dest='chan_axis', default='index',
    help='Choose channel axis to be channel/delay index (index), or physical coordinates (physical).  Default is index.')
o.add_option('--clean', dest='clean', action='store_true',
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.')

def convert_arg_range(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    arg = arg.split(',')
    return [map(float, option.split('_')) for option in arg]

def data_selector(antopt, active_pol, uv):
    """Call uv.select with appropriate options based on command-line argument 
    for ants."""
    if antopt.startswith('all'): pass
    elif antopt.startswith('auto'): uv.select('auto', 0, 0)
    elif antopt.startswith('cross'): uv.select('auto', 0, 0, include=0)
    else:
        antopt = convert_arg_range(antopt)
        for opt in antopt:
            try: a1,a2 = opt
            except(ValueError): a1,a2 = opt + [-1]
            uv.select('antennae', a1, a2)
    uv.select('polarization', active_pol, 0)

def gen_chans(chanopt, uv, coords, is_delay):
    """Return an array of active channels and whether or not a range of
    channels is selected (as opposed to one or more individual channels)
    based on command-line arguments."""
    is_chan_range = True
    if chanopt == 'all': chans = n.arange(uv['nchan'])
    else:
        chanopt = convert_arg_range(chanopt)
        if coords != 'index':
            if is_delay:
                def conv(c):
                    return int(n.round(c * uv['sdf'] * uv['nchan'])) \
                        + uv['nchan']/2
            else:
                def conv(c): return int(n.round((c - uv['sfreq']) / uv['sdf']))
        else:
            if is_delay:
                def conv(c): return int(c) + uv['nchan']/2
            else:
                def conv(c): return c
        chanopt = [map(conv, c) for c in chanopt]
        if len(chanopt[0]) != 1: 
            chanopt = [n.arange(x,y, dtype=n.int) for x,y in chanopt]
        else: is_chan_range = False
        chans = n.concatenate(chanopt)
    return chans.astype(n.int), is_chan_range

def gen_times(timeopt, uv, coords, decimate, is_fringe):
    is_time_range = True
    if timeopt == 'all' or is_fringe:
        def time_selector(t, cnt): return True
    else:
        timeopt = convert_arg_range(timeopt)
        if len(timeopt[0]) != 1:
            def time_selector(t, cnt):
                if coords == 'index': t = cnt
                for opt in timeopt:
                    if (t >= opt[0]) and (t < opt[1]): return True
                return False
        else:
            is_time_range = False
            timeopt = [opt[0] for opt in timeopt]
            inttime = uv['inttime'] / a.const.s_per_day * decimate
            def time_selector(t, cnt):
                if coords == 'index': return cnt in timeopt
                for opt in timeopt:
                    if (t >= opt) and (t < opt + inttime): return True
                return False
    return time_selector, is_time_range

# Initialized based on command line arguments
opts, args = o.parse_args(sys.argv[1:])
uv = a.miriad.UV(args[0])
if opts.pol is None: active_pol = uv['pol']
else: active_pol = a.miriad.str2pol[opts.pol]
chans, is_chan_range = gen_chans(opts.chan, uv, opts.chan_axis, opts.delay)
freqs = n.arange(uv['sfreq'], uv['sfreq']+uv['nchan']*uv['sdf'], uv['sdf'])
freqs = freqs.take(chans)
delays = n.arange(-.5/uv['sdf'], .5/uv['sdf'], 1/(uv['sdf']*uv['nchan']))
delays = delays.take(chans)
time_sel, is_time_range = gen_times(opts.time, uv, opts.time_axis, 
    opts.decimate, opts.fringe)
inttime = uv['inttime'] * opts.decimate
del(uv)

# Loop through UV files collecting relevant data
plot_x = {}
plot_t = {'jd':[], 'lst':[], 'cnt':[]}
times = []

for uvfile in args:
    print 'Reading', uvfile
    uv = a.miriad.UV(uvfile)
    # Only select data that is needed to plot
    data_selector(opts.ants, active_pol, uv)
    # Read data from a single UV file
    for (uvw,t,(i,j)),d in uv.all():
        bl = '%d,%d' % (i,j)
        # Implement Decimation
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            use_this_time = ((len(times) - 1) % opts.decimate) == 0
            use_this_time &= time_sel(t, (len(times)-1) / opts.decimate)
            if use_this_time:
                plot_t['lst'].append(uv['lst'])
                plot_t['jd'].append(t)
                plot_t['cnt'].append((len(times)-1) / opts.decimate)
        if not use_this_time: continue
        # Do delay transform if required
        if opts.delay:
            if opts.unmask:
                d = d.data
                ker = n.zeros_like(d)
                ker[0] = 1.
                gain = 1.
            else:
                flags = n.logical_not(d.mask).astype(n.float)
                ker = n.fft.ifft(flags)
                d = d.filled(0)
                gain = flags.sum()
            flags = n.where
            d = n.fft.ifft(d)
            if opts.clean and not n.all(d == 0):
                d, info = a.deconv.clean1d(d, ker, tol=1e-10)
                d += info['res'] / gain
            d = n.ma.array(d)
            d = n.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], 
                axis=0)
        elif opts.unmask: d = d.data
        # Extract specific channels for plotting
        d = d.take(chans)
        d.shape = (1,) + d.shape
        if not plot_x.has_key(bl): plot_x[bl] = []
        plot_x[bl].append(d)
    del(uv)

bls = plot_x.keys()
bls.sort()
if len(bls) == 0:
    print 'No data to plot.'
    sys.exit(0)
m2 = int(math.sqrt(len(bls)))
m1 = int(math.ceil(float(len(bls)) / m2))

# Generate all the plots
for cnt, bl in enumerate(bls):
    d = n.ma.concatenate(plot_x[bl], axis=0)
    if opts.df: d = d[:,:-2]/2 + d[:,2:]/2 - d[:,1:-1]
    if opts.dt: d = d[:-2]/2 + d[2:]/2 - d[1:-1]
    if opts.fringe:
        d = n.ma.array(n.fft.ifft(d.filled(0), axis=0))
        d = n.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], axis=0)
    if opts.sum_chan:
        d = d.sum(axis=1)
        is_chan_range = False
    if opts.mode.startswith('phs'): d = n.angle(d.filled(0))
    elif opts.mode.startswith('lin'): d = n.ma.absolute(d)
    elif opts.mode.startswith('real'): d = d.real
    elif opts.mode.startswith('imag'): d = d.imag
    elif opts.mode.startswith('log'):
        d = n.ma.absolute(d.filled(0))
        d = n.ma.masked_less_equal(d, 0)
        d = n.ma.log10(d)
    else: raise ValueError('Unrecognized plot mode.')
    p.subplot(m1, m2, cnt+1)
    if is_chan_range and is_time_range:
        if opts.fringe:
            if opts.time_axis == 'index':
                t1 = len(plot_t['jd'])/2 - len(plot_t['jd'])
                t2 = len(plot_t['jd'])/2
                ylabel = 'Fringe Rate (bins)'
            else:
                t1 = -500/inttime
                t2 =  500/inttime - 1000 / (inttime * len(plot_t['jd']))
                ylabel = 'Fringe Rate (milliHz)'
        else:
            if opts.time_axis == 'index':
                t1,t2 = plot_t['cnt'][0], plot_t['cnt'][-1]
                ylabel = 'Time (integrations)'
            else:
                t1,t2 = plot_t['jd'][0], plot_t['jd'][-1]
                ylabel = 'Time (Julian Date)'
        if opts.delay:
            if opts.chan_axis == 'index':
                c1,c2 = len(chans)/2 - len(chans), len(chans)/2
                xlabel = 'Delay (bins)'
            else:
                c1,c2 = delays[0], delays[-1]
                xlabel = 'Delay (ns)'
        else:
            if opts.chan_axis == 'index':
                c1,c2 = 0, len(chans) - 1
                xlabel = 'Frequency (chan)'
            else:
                c1,c2 = freqs[0], freqs[-1]
                xlabel = 'Frequency (GHz)'
        if opts.mode.startswith('log'):
            if not opts.plot_min is None:
                p.imshow(d, extent=(c1,c2,t2,t1), aspect='auto', 
                    vmin=opts.plot_min)
            else: p.imshow(d, extent=(c1,c2,t2,t1), aspect='auto')
        else: p.imshow(d, extent=(c1,c2,t2,t1), aspect='auto')
        p.colorbar()
        p.xlabel(xlabel)
        p.ylabel(ylabel)
    elif is_chan_range and not is_time_range:
        if opts.delay:
            if opts.chan_axis == 'index':
                plot_chans = range(len(chans)/2 - len(chans), len(chans)/2)
                xlabel = 'Delay (bins)'
            else:
                plot_chans = delays
                xlabel = 'Delay (ns)'
        else:
            if opts.chan_axis == 'index':
                plot_chans = chans
                xlabel = 'Frequency (chan)'
            else:
                plot_chans = freqs
                xlabel = 'Frequency (GHz)'
        if opts.time_axis == 'index':
            plot_t = plot_t['cnt']
            label = '#%d'
        else:
            plot_t = plot_t['jd']
            label = 'jd%f'
        for i,t in enumerate(plot_t):
            p.plot(plot_chans, d[i,:], '-', label=label % t)
        p.xlabel(xlabel)
    elif not is_chan_range and is_time_range:
        if opts.time_axis == 'index': plot_times = range(len(plot_t['jd']))
        elif opts.time_axis == 'physical': plot_times = plot_t['jd']
        elif opts.time_axis == 'lst': plot_times = plot_t['lst']
        else: raise ValueError('Unrecognized time axis type.')
        if opts.sum_chan: p.plot(plot_times, d, '.', label='(+)')
        else:
            if opts.chan_axis == 'index': label = '#%d'
            else:
                chans = freqs
                label = '%f GHz'
            for c, chan in enumerate(chans):
                p.plot(plot_times, d[:,c], '.', label=label % chan)
    else: raise ValueError('Either time or chan needs to be a range.')
    p.title(bl)
if not is_time_range or not is_chan_range: p.legend(loc='best')

# Save to a file or pop up a window
if opts.out_file != '': p.savefig(opts.out_file)
else: p.show()

