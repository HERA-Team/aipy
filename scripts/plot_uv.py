#!/usr/bin/env python

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

from __future__ import print_function, division, absolute_import

import aipy as a, numpy as np, sys, optparse
from matplotlib import pylab as p

o = optparse.OptionParser()
o.set_usage('plot_uv.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True,ant=True, pol=True, chan=True, dec=True,
    cmap=True, max=True, drng=True, cal=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
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
o.add_option('--xlim', dest='xlim',
    help='Limits on the x axis (channel/delay) for plotting.')
o.add_option('--ylim', dest='ylim',
    help='Limits on the x axis (time/delay-rate) for plotting.')
o.add_option('--plot_each', dest='plot_each',
    help='Instead of a waterfall plot, plot each of the specified axis (chan,time)')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))

def convert_arg_range(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    arg = arg.split(',')
    return [map(float, option.split('_')) for option in arg]

def gen_times(timeopt, uv, coords, decimate):
    if timeopt == 'all':
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
            timeopt = [opt[0] for opt in timeopt]
            inttime = uv['inttime'] / a.const.s_per_day * decimate
            def time_selector(t, cnt):
                if coords == 'index': return cnt in timeopt
                for opt in timeopt:
                    if (t >= opt) and (t < opt + inttime): return True
                return False
    return time_selector

def data_mode(data, mode='abs'):
    if mode.startswith('phs'): data = np.angle(data.filled(0))
    elif mode.startswith('lin'):
        data = np.ma.absolute(data.filled(0))
        data = np.ma.masked_less_equal(data, 0)
    elif mode.startswith('real'): data = data.real
    elif mode.startswith('imag'): data = data.imag
    elif mode.startswith('log'):
        data = np.ma.absolute(data.filled(0))
        data = np.ma.masked_less_equal(data, 0)
        data = np.ma.log10(data)
    else: raise ValueError('Unrecognized plot mode.')
    return data

opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
cmap = p.get_cmap(opts.cmap)
if not opts.xlim == None: opts.xlim = map(float, opts.xlim.split('_'))
if not opts.ylim == None: opts.ylim = map(float, opts.ylim.split('_'))
uv = a.miriad.UV(args[0])
a.scripting.uv_selector(uv, opts.ant, opts.pol)
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
is_chan_range, is_time_range = True, True
if opts.plot_each == 'chan': is_chan_range = False
elif opts.plot_each == 'time': is_time_range = False
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
freqs = freqs.take(chans)
if opts.delay:
    if freqs.size == freqs[-1] - freqs[0] + 1:
        # XXX someday could allow for equal spaced chans
        raise ValueError('Channels must be contiguous to do delay transform (chan=%s)' % (opts.chan))
    delays = np.fft.fftfreq(freqs.size, freqs[1]-freqs[0])
    delays = np.fft.fftshift(delays)
time_sel = gen_times(opts.time, uv, opts.time_axis, opts.decimate)
inttime = uv['inttime'] * opts.decimate
if not opts.src is None:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]
del(uv)

# Loop through UV files collecting relevant data
plot_x = {}
plot_t = {'jd':[], 'lst':[], 'cnt':[]}
times = []

# Hold plotting handles
plots = {}
plt_data = {}

for uvfile in args:
    print('Reading', uvfile)
    uv = a.miriad.UV(uvfile)
    if not opts.cal is None:
        aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
        aa.set_active_pol(opts.pol)
        aa.select_chans(chans)
    else: aa = None
    # Only select data that is needed to plot
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    # Read data from a single UV file
    for (uvw,t,(i,j)),d in uv.all():
        bl = '%d,%d,%d' % (i,j,uv['pol'])
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            # Implement time selection
            use_this_time = time_sel(t, (len(times)-1) / opts.decimate)
            if use_this_time:
                if aa == None: lst = uv['lst']
                else:
                    aa.set_jultime(t)
                    lst = aa.sidereal_time()
                plot_t['lst'].append(lst)
                plot_t['jd'].append(t)
                plot_t['cnt'].append(len(times)-1)
        if not use_this_time: continue
        d = d.take(chans)
        if opts.ant.find('%d_%d' % (j,i)) != -1: d = d.conj() # obey antenna ordering, if possible
        #apply cal phases
        if not opts.cal is None:
            aa.set_jultime(t)
            if not opts.src is None:
                src.compute(aa)
                d = aa.phs2src(d, src, i, j)
            #else: took out this mode because it's not used, and prefer not to phase.
            #    d *= np.exp(-1j*np.pi*aa.get_phs_offset(i,j))
        # Do delay transform if required
        if opts.delay:
            w = a.dsp.gen_window(d.shape[-1], window=opts.window)
            if opts.unmask:
                flags = np.ones(d.shape, dtype=np.float)
                d = d.data
            else:
                flags = np.logical_not(d.mask).astype(np.float)
                d = d.filled(0)
            d = np.fft.ifft(d*w)
            ker = np.fft.ifft(flags*w)
            gain = a.img.beam_gain(ker)
            if not opts.clean is None and not np.all(d == 0):
                d, info = a.deconv.clean(d, ker, tol=opts.clean)
                d += info['res'] / gain
            d = np.ma.array(d)
            d = np.fft.fftshift(d, axes=0)
        elif opts.unmask: d = d.data
        d.shape = (1,) + d.shape
        if not plot_x.has_key(bl): plot_x[bl] = []
        plot_x[bl].append(d)
    del(uv)

bls = plot_x.keys()
def sort_func(a, b):
    ai,aj,pa = map(int, a.split(','))
    bi,bj,pb = map(int, b.split(','))
    if bi > ai or (bi == ai and bj > aj) or (bi == ai and bj == aj and pb < pa): return -1
    return 1
bls.sort(cmp=sort_func)
if len(bls) == 0:
    print('No data to plot.')
    sys.exit(0)
m2 = int(np.sqrt(len(bls)))
m1 = int(np.ceil(float(len(bls)) / m2))

# Generate all the plots
dmin,dmax = None, None
fig = p.figure()
if not opts.src is None:fig.suptitle(opts.src)
for cnt, bl in enumerate(bls):
    d = np.ma.concatenate(plot_x[bl], axis=0)
    i,j,pol = map(int,bl.split(','))
    if opts.df: d = d[:,:-2]/2 + d[:,2:]/2 - d[:,1:-1]
    if opts.dt: d = d[:-2]/2 + d[2:]/2 - d[1:-1]
    if opts.fringe:
        d = d.filled(0)
        w = a.dsp.gen_window(d.shape[0], window=opts.window); w.shape += (1,)
        wgts = np.where(d != 0, 1., 0.) * w
        gain = np.sqrt(np.average(wgts**2, axis=0))
        ker = np.fft.ifft(wgts, axis=0) # w already put in 2 lines above
        d = np.fft.ifft(d*w, axis=0)
        if not opts.clean is None:
            for chan in range(d.shape[1]):
                if gain[chan] == 0: continue
                d[:,chan],info = a.deconv.clean(d[:,chan],ker[:,chan],tol=opts.clean)
                d[:,chan] += info['res'] / gain[chan]
        d = np.fft.fftshift(d, axes=0)
        d = np.ma.array(d)
    plt_data[cnt+1] = d
    d = data_mode(d, opts.mode)
    if not opts.share:
        p.subplot(m2, m1, cnt+1)
        dmin,dmax = None,None
        label = ''
    else:
        pol = a.miriad.pol2str[pol]
        label = '%d%s,%d%s ' % (i,pol[0],j,pol[-1])
    if is_chan_range and is_time_range:
        if opts.fringe:
            if opts.time_axis == 'index':
                drates = np.fft.fftfreq(len(plot_t['cnt']), 1./len(plot_t['cnt']))
                step = drates[1] - drates[0]
                ylabel = 'Delay-Rate (bins)'
            else:
                drates = np.fft.fftfreq(len(plot_t['cnt']), inttime) * 1e3 # mHz
                step = drates[1] - drates[0]
                ylabel = 'Delay-Rate (milliHz)'
            drates = np.fft.fftshift(drates)
            t1,t2 = drates[0]-0.5*step,drates[-1]+0.5*step
        else:
            if opts.time_axis == 'index':
                step = 1.
                t1,t2 = plot_t['cnt'][0]-0.5*step, plot_t['cnt'][-1]+0.5*step
                ylabel = 'Time (integrations)'
            elif opts.time_axis=='lst':
                step = plot_t['lst'][1] - plot_t['lst'][0]
                t1,t2 = (plot_t['lst'][0]-0.5*step)*12/np.pi, (plot_t['lst'][-1]+0.5*step)*12/np.pi
                ylabel = 'Local Sideral time (hrs)'
            else:
                step = plot_t['jd'][1] - plot_t['jd'][0]
                t1,t2 = plot_t['jd'][0]-0.5*step, plot_t['jd'][-1]+0.5*step
                ylabel = 'Time (Julian Date)'
        if opts.delay:
            if opts.chan_axis == 'index':
                step = 1
                c1,c2 = len(chans)/2 - len(chans)-0.5*step, len(chans)/2-0.5*step
                xlabel = 'Delay (bins)'
            else:
                step = delays[1] - delays[0]
                c1,c2 = delays[0]-0.5*step, delays[-1]+0.5*step
                xlabel = 'Delay (ns)'
        else:
            if opts.chan_axis == 'index':
                step = 1
                c1,c2 = 0-0.5*step, len(chans)-0.5*step
                xlabel = 'Frequency (chan)'
            else:
                step = freqs[1] - freqs[0]
                c1,c2 = freqs[0]-0.5*step, freqs[-1]+0.5*step
                xlabel = 'Frequency (GHz)'
        if not opts.max is None: dmax = opts.max
        elif dmax is None: dmax = d.max()
        else: dmax = max(dmax,d.max())
        if not opts.drng is None: dmin = dmax - opts.drng
        elif dmin is None: dmin = d.min()
        else: dmin = min(dmin,d.min())
        plots[cnt+1] = p.imshow(d, extent=(c1,c2,t2,t1), origin='upper',
            aspect='auto', interpolation='nearest',
            vmax=dmax, vmin=dmin, cmap=cmap)
        p.colorbar(shrink=0.5)
        p.xlabel(xlabel); p.ylabel(ylabel)
        if not opts.xlim == None: p.xlim(*opts.xlim)
        if not opts.ylim == None: p.ylim(opts.ylim[1],opts.ylim[0]) # Reverse b/c otherwise ylim flips origin for unknown reasons
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
        for ti,t in enumerate(plot_t):
            p.plot(plot_chans, d[ti,:], '-', label=label % t)
        p.xlabel(xlabel)
        if not opts.xlim == None: p.xlim(*opts.xlim)
        if not opts.max is None: dmax = opts.max
        elif dmax is None: dmax = d.max()
        else: dmax = max(dmax,d.max())
        if not opts.drng is None: dmin = dmax - opts.drng
        elif dmin is None: dmin = d.min()
        else: dmin = min(dmin,d.min())
        if not opts.ylim == None: p.ylim(*opts.ylim)
        else: p.ylim(dmin,dmax)
    elif not is_chan_range and is_time_range:
        if opts.fringe:
            if opts.time_axis == 'index':
                drates = np.fft.fftfreq(len(plot_t['cnt']), 1./len(plot_t['cnt']))
                xlabel = 'Delay-Rate (bins)'
            else:
                print(inttime, len(plot_t['cnt']))
                drates = np.fft.fftfreq(len(plot_t['cnt']), inttime) * 1e3 # mHz
                xlabel = 'Delay-Rate (milliHz)'
            plot_times = np.fft.fftshift(drates)
        else:
            if opts.time_axis == 'index':
                plot_times = range(len(plot_t['jd']))
                xlabel = 'Time (integrations)'
            elif opts.time_axis == 'physical':
                plot_times = plot_t['jd']
                xlabel = 'Time (JD)'
            elif opts.time_axis == 'lst':
                plot_times = plot_t['lst']
                xlabel = 'Time (Sidereal Radians)'
            else: raise ValueError('Unrecognized time axis type.')
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
        p.xlabel(xlabel)
        if not opts.xlim == None: p.xlim(*opts.xlim)
        if not opts.ylim == None: p.ylim(*opts.ylim)
        else: p.ylim(dmin,dmax)
    else: raise ValueError('Either time or chan needs to be a range.')
    if not opts.share:
        pol = a.miriad.pol2str[pol]
        title = '%d%s,%d%s ' % (i,pol[0],j,pol[-1])
        p.title(title)
if not opts.nolegend and (not is_time_range or not is_chan_range):
    p.legend(loc='best')

# Save to a file or pop up a window
if opts.out_file != '': p.savefig(opts.out_file)
else:
    def click(event):
        print([event.key])
        if event.key == 'm':
            mode = raw_input('Enter new mode: ')
            for k in plots:
                try:
                    d = data_mode(plt_data[k], mode)
                    plots[k].set_data(d)
                except(ValueError):
                    print('Unrecognized plot mode')
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
            print('Replotting...')
            p.draw()
    p.connect('key_press_event', click)
    p.show()
