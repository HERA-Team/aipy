#! /usr/bin/env python
"""
Creates waterfall plots from Miriad UV files.  Can tile multiple plots
on one window, or plot just a single baseline.  When taking the delay
transform (-d), channels are interepreted as selecting delays after the
the transform operation.  Similarly, the time select (-t) will be interpreted
as selecting fringe rates if the fringe-rate transform (-f) has been selected.
In both cases, the ranges specified are intepreted to be the units of the
output plot (i.e. as specified by --chan_axis and --time_axis).

Author: Aaron Parsons, Griffin Foster
"""

import aipy as a, numpy as n, pylab as p, math, sys, optparse

o = optparse.OptionParser()
o.set_usage('plot_uv.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True, chan=True, dec=True,
    cmap=True, max=True, drng=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('--sum_chan', dest='sum_chan', action='store_true',
    help='Sum active channels together.')
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
o.add_option('--time_axis', dest='time_axis', default='index',
    help='Choose time axis to be integration/fringe index (index), or physical coordinates (physical), or if doing xy plot in time-mode, (lst) is also available.  Default is index.')
o.add_option('--chan_axis', dest='chan_axis', default='index',
    help='Choose channel axis to be channel/delay index (index), or physical coordinates (physical).  Default is index.')
o.add_option('--clean', dest='clean', type='float',
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('--nolegend', dest='nolegend', action='store_true',
    help='Omit legend in last plot.')
o.add_option('--share', dest='share', action='store_true',
    help='Share plots in a single frame.')
o.add_option('--window', dest='window', default='kaiser3',
    help='Windowing function to use in delay transform.  Default is kaiser3.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))

def convert_arg_range(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    arg = arg.split(',')
    return [map(float, option.split('_')) for option in arg]

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

def data_mode(data, mode='abs'):
    if mode.startswith('phs'): data = n.angle(data.filled(0))
    elif mode.startswith('lin'):
        data = n.ma.absolute(data.filled(0))
        data = n.ma.masked_less_equal(data, 0)
    elif mode.startswith('real'): data = data.real
    elif mode.startswith('imag'): data = data.imag
    elif mode.startswith('log'):
        data = n.ma.absolute(data.filled(0))
        data = n.ma.masked_less_equal(data, 0)
        data = n.ma.log10(data)
    else: raise ValueError('Unrecognized plot mode.')
    return data

opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
cmap = p.get_cmap(opts.cmap)
uv = a.miriad.UV(args[0])
a.scripting.uv_selector(uv, opts.ant, opts.pol)
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

# Hold plotting handles
plots = {}
plt_data = {}

if opts.cal != None:
    uv = a.miriad.UV(args[0])
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    del(uv)
else: aa = None

for uvfile in args:
    print 'Reading', uvfile
    uv = a.miriad.UV(uvfile)
    # Only select data that is needed to plot
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    # Read data from a single UV file
    for (uvw,t,(i,j)),d in uv.all():
        bl = '%d,%d' % (i,j)
        # Implement Decimation
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            use_this_time = ((len(times) - 1) % opts.decimate) == 0
            use_this_time &= time_sel(t, (len(times)-1) / opts.decimate)
            if use_this_time:
                if aa == None: lst = uv['lst']
                else:
                    aa.set_jultime(t)
                    lst = aa.sidereal_time()
                plot_t['lst'].append(lst)
                plot_t['jd'].append(t)
                plot_t['cnt'].append((len(times)-1) / opts.decimate)
        if not use_this_time: continue
        # Do delay transform if required
        if opts.delay:
            w = a.dsp.gen_window(d.shape[-1], window=opts.window)
            if opts.unmask:
                d = d.data
                ker = n.zeros_like(d)
                ker[0] = 1.
                gain = 1.
            else:
                flags = n.logical_not(d.mask).astype(n.float)
                gain = n.sqrt(n.average(flags**2))
                ker = n.fft.ifft(flags*w)
                d = d.filled(0)
            d = n.fft.ifft(d*w)
            if not opts.clean is None and not n.all(d == 0):
                d, info = a.deconv.clean(d, ker, tol=opts.clean)
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
def sort_func(a, b):
    ai,aj = map(int, a.split(','))
    bi,bj = map(int, b.split(','))
    if bi > ai or (bi == ai and bj > aj): return -1
    return 1
bls.sort(cmp=sort_func)
if len(bls) == 0:
    print 'No data to plot.'
    sys.exit(0)
m2 = int(math.sqrt(len(bls)))
m1 = int(math.ceil(float(len(bls)) / m2))

# Generate all the plots
dmin,dmax = None, None
for cnt, bl in enumerate(bls):
    d = n.ma.concatenate(plot_x[bl], axis=0)
    if opts.df: d = d[:,:-2]/2 + d[:,2:]/2 - d[:,1:-1]
    if opts.dt: d = d[:-2]/2 + d[2:]/2 - d[1:-1]
    if opts.fringe:
        d = d.filled(0)
        flags = n.where(d[:,0] != 0, 1., 0.)
        gain = n.sqrt(n.average(flags**2))
        ker = n.fft.ifft(flags)
        d = n.fft.ifft(d, axis=0)
        if not opts.clean is None:
            for chan in range(d.shape[1]):
                d[:,chan],info = a.deconv.clean(d[:,chan],ker,tol=opts.clean)
                d[:,chan] += info['res'] / gain
        d = n.ma.concatenate([d[d.shape[0]/2:], d[:d.shape[0]/2]], axis=0)
    if opts.sum_chan:
        d = d.sum(axis=1)
        is_chan_range = False
    plt_data[cnt+1] = d
    d = data_mode(d, opts.mode)
    if not opts.share:
        p.subplot(m2, m1, cnt+1)
        dmin,dmax = None,None
        label = ''
    else: label = bl + ' '
    if is_chan_range and is_time_range:
        if opts.fringe:
            if opts.time_axis == 'index':
                if opts.time != 'all':
                    t1, t2 = map(float, opts.time.split('_'))
                    d = d[t1+d.shape[0]/2:t2+d.shape[0]/2]
                else:
                    t1 = len(plot_t['jd'])/2 - len(plot_t['jd'])
                    t2 = len(plot_t['jd'])/2
                ylabel = 'Delay-Rate (bins)'
            else:
                t1 = -500/inttime
                t2 =  500/inttime - 1000 / (inttime * len(plot_t['jd']))
                ylabel = 'Delay-Rate (milliHz)'
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
        if not opts.max is None: dmax = opts.max
        elif dmax is None: dmax = d.max()
        else: dmax = max(dmax,d.max())
        if not opts.drng is None: dmin = dmax - opts.drng
        elif dmin is None: dmin = d.min()
        else: dmin = min(dmin,d.min())
        plots[cnt+1] = p.imshow(d, extent=(c1,c2,t2,t1), 
            aspect='auto', interpolation='nearest', 
            vmax=dmax, vmin=dmin, cmap=cmap)
        p.colorbar(shrink=0.5)
        p.xlabel(xlabel); p.ylabel(ylabel)
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
            if cnt == 0: plot_t = plot_t['cnt']
            label += '#%d'
        else:
            if cnt == 0: plot_t = plot_t['jd']
            label += 'jd%f'
        for i,t in enumerate(plot_t):
            p.plot(plot_chans, d[i,:], '-', label=label % t)
        p.xlabel(xlabel)
        if not opts.max is None: dmax = opts.max
        elif dmax is None: dmax = d.max()
        else: dmax = max(dmax,d.max())
        if not opts.drng is None: dmin = dmax - opts.drng
        elif dmin is None: dmin = d.min()
        else: dmin = min(dmin,d.min())
        p.ylim(dmin,dmax)
    elif not is_chan_range and is_time_range:
        if opts.time_axis == 'index': plot_times = range(len(plot_t['jd']))
        elif opts.time_axis == 'physical': plot_times = plot_t['jd']
        elif opts.time_axis == 'lst': plot_times = plot_t['lst']
        else: raise ValueError('Unrecognized time axis type.')
        if opts.sum_chan: p.plot(plot_times, d, '-', label=label+'(+)')
        else:
            if opts.chan_axis == 'index': label += '#%d'
            else:
                chans = freqs
                label += '%f GHz'
            for c, chan in enumerate(chans):
                p.plot(plot_times, d[:,c], '-', label=label % chan)
        if not opts.max is None: dmax = opts.max
        elif dmax is None: dmax = d.max()
        else: dmax = max(dmax,d.max())
        if not opts.drng is None: dmin = dmax - opts.drng
        elif dmin is None: dmin = d.min()
        else: dmin = min(dmin,d.min())
        p.ylim(dmin,dmax)
    else: raise ValueError('Either time or chan needs to be a range.')
    if not opts.share: p.title(bl)
if not opts.nolegend and (not is_time_range or not is_chan_range): 
    p.legend(loc='best')

# Save to a file or pop up a window
if opts.out_file != '': p.savefig(opts.out_file)
else:
    def click(event):
        print [event.key]
        if event.key == 'm':
            mode = raw_input('Enter new mode: ')
            for k in plots:
                try:
                    d = data_mode(plt_data[k], mode)
                    plots[k].set_data(d)
                except(ValueError):
                    print 'Unrecognized plot mode'
            p.draw()
        elif event.key == 'd':
            max = raw_input('Enter new max: ')
            try: max = float(max)
            except(ValueError): max = None
            drng = raw_input('Enter new drng: ')
            try: drng = float(drng)
            except(ValueError): drng = None
            for k in plots:
                _max,_drng = max, drng
                if _max is None or _drng is None:
                    d = plots[k].get_array()
                    if _max is None: _max = d.max()
                    if _drng is None: _drng = _max - d.min()
                plots[k].set_clim(vmin=_max-_drng, vmax=_max)
            print 'Replotting...'
            p.draw()
    p.connect('key_press_event', click)
    p.show()
