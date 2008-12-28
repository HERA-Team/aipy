#!/usr/bin/python
"""
A script for detecting and flagging RFI related effects in UV files.

Author: Aaron Parsons
Date: 5/31/07
Revisions:
    10/10/07    arp Lowered thresh for flagging an integration to coincidence
                    on 2 antennas.
    12/11/07 arp    Ported to use new miriad file interface
"""

__version__ = '0.0.3'

import numpy, aipy

def gaussian(amp, sig, off, x):
    """Generate gaussian value at x given amplitude, sigma, and x offset."""
    return amp * numpy.exp(-((x - off)/sig)**2)

def fit_gaussian(xs, ys):
    """Fit a gaussian (returning amplitude, sigma, offset) to a set of x
    values and their corresponding y values."""
    def fitfunc(args):
        amp, sig, off = args
        sum = ((ys - gaussian(amp, sig, off, xs))**2).sum()
        return sum
    return aipy.optimize.fmin(fitfunc, (ys.sum(), 50, ys.argmax()), disp=0)

def gen_rfi_thresh(data, nsig=2, per_bin=1000):
    """Generate a threshold at nsig times the standard deviation above the mean
    above which data is considered rfi."""
    try:
        histdata = numpy.log10(numpy.abs(data.compressed()))
        h, bvals = numpy.histogram(histdata, bins=histdata.size/per_bin)
    except(ValueError): return 0, 0
    if h.size == 0: return 0, 0
    # Fit a gaussian to histogram (better than just std-dev of data)
    amp, sig, off = fit_gaussian(numpy.arange(h.size), h)
    sig = abs(sig)
    hi_thresh = numpy.clip(numpy.round(off + nsig*sig), 0, len(bvals)-1)
    lo_thresh = numpy.clip(numpy.round(off - nsig*sig), 0, len(bvals)-1)
    hi_thresh = bvals[hi_thresh]
    lo_thresh = bvals[lo_thresh]
    return 10**hi_thresh, 10**lo_thresh

def flag_by_int(preflagged_auto, nsig=1):
    """Flag rfi for an autocorrelation.  Iteratively removes outliers and
    fits a smooth passband, then uses smooth passband to remove outliers
    (both positive and negative)."""
    pwr_vs_t = numpy.ma.average(abs(preflagged_auto), axis=1)
    spikey_pwr_vs_t = abs(pwr_vs_t - remove_spikes(pwr_vs_t))
    hi_thr, lo_thr = gen_rfi_thresh(spikey_pwr_vs_t, per_bin=20, nsig=nsig)
    mask = spikey_pwr_vs_t > hi_thr
    return mask.astype(numpy.int)

def remove_spikes(data, order=6, iter=3, return_poly=False):
    """Iteratively fits a smooth function by removing outliers and fitting
    a polynomial, then using the polynomial to remove other outliers."""
    xs = numpy.arange(data.size)
    if numpy.ma.isMA(data) and data.mask.size == data.size: mask = data.mask
    else: mask = numpy.zeros_like(data.data)
    im = numpy.logical_not(mask)
    nxs = numpy.compress(im, xs)
    ndata = numpy.compress(im, data)
    if len(nxs) != 0: p = numpy.polyfit(nxs, ndata, deg=order)
    else: p = numpy.polyfit(xs, data, deg=order)
    if iter != 0:
        residual = abs(data - numpy.polyval(p, xs))
        sig = numpy.sqrt(numpy.ma.average(residual**2))
        mask.put(numpy.where(residual > sig/2), 1)
        data = numpy.ma.array(data, mask=mask)
        p = remove_spikes(data, order=order, iter=iter-1, return_poly=True)
    if return_poly: return p
    else: return numpy.polyval(p, xs)

if __name__ == '__main__':
    import os, sys, pickle
    from optparse import OptionParser

    o = OptionParser()
    o.set_usage('xrfi.py [options] *.uv')
    o.set_description(__doc__)
    o.add_option('-n', '--nsig', dest='nsig', default=2., type='float',
        help='Number of standard deviations above mean to flag.  Default 2.')
    o.add_option('-m', '--flagmode', dest='flagmode', default='both',
        help='Can be val,int,both for flagging by value only, integration only, or both.  Default both.')
    o.add_option('-c', '--chans', dest='chans', default='',
        help='Comma-delimited ranges (e.g. 1-3,5-10) of channels to manually flag before any other statistical flagging.')
    o.add_option('-t', '--ch_thresh', dest='ch_thresh',type='float',default=.33,
        help='Fraction of the data in a channel which, if flagged, will result in the entire channel begin flagged.  Default .33')
    o.add_option('-i', '--infile', dest='infile', action='store_true',
        help='Apply xrfi flags generated with the -o option.')
    o.add_option('-o', '--outfile', dest='outfile', action='store_true',
        help='Rather than apply the flagging to the data, store them in a file (named by JD) to apply to a different file with the same JD.')


    opts, args = o.parse_args(sys.argv[1:])
    chans = opts.chans.split(',')
    flag_chans = []
    for c in chans:
        if len(c) == 0: continue
        c = map(int, c.split('-'))
        flag_chans.extend(range(c[0], c[1]+1))

    for uvfile in args:
        print 'Working on', uvfile
        uvofile = uvfile+'r'
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
        uvi = aipy.miriad.UV(uvfile)
        (uvw,jd,(i,j)),d,f = uvi.read(raw=True)
        uvi.rewind()
        if opts.infile:
            if not os.path.exists('%f.xrfi' % jd):
                print '%f.xrfi' % jd, 'does not exist.  Skipping...'
                continue
            f = open('%f.xrfi' % jd)
            mask = pickle.load(f)
            f.close()
            for m in mask: mask[m] = numpy.array(mask[m])
        else:
            # Gather all data and each time step
            window = None
            data = {}
            mask = {}
            times = []
            for (uvw,t,(i,j)), d, f in uvi.all(raw=True):
                if len(times) == 0 or times[-1] != t: times.append(t)
                mask[t] = mask.get(t, 0) | f
                bl = aipy.miriad.ij2bl(i,j)
                if i != j:
                    if window is None: 
                        window = d.size/2 - abs(numpy.arange(d.size) - d.size/2)
                    d = numpy.fft.fft(numpy.fft.ifft(d) * window)
                pol = uvi['pol']
                if not pol in data: data[pol] = {}
                if not bl in data[pol]: data[pol][bl] = {}
                data[pol][bl][t] = d

            # Generate a single mask for all baselines which masks if any
            # baseline has an outlier at that freq/time.  Data flagged
            # strictly on basis of nsigma above mean.
            if not opts.flagmode.startswith('int'):
                new_mask = {}
                for k in mask:
                    mask[k][flag_chans] = 1
                    new_mask[k] = mask[k].copy()
                for p in data:
                  for k in data[p]:
                    i, j = aipy.miriad.bl2ij(k)
                    if i == j: continue
                    data_times = data[p][k].keys()
                    d = numpy.ma.array([data[p][k][t] for t in data_times],
                        mask=[mask[t] for t in data_times])
                    hi_thr, lo_thr = gen_rfi_thresh(d, nsig=opts.nsig)
                    m = numpy.where(numpy.abs(d) > hi_thr,1,0)
                    for i, t in enumerate(data_times): new_mask[t] |= m[i]
                mask = new_mask
                # If more than half the data in a channel is flagged, flag the
                # whole thing
                msk_cnt = numpy.array([mask[t] for t in data_times]).sum(axis=0)
                ch_msk = numpy.where(msk_cnt > msk_cnt.max()*opts.ch_thresh,1,0)
                for k in mask: mask[k] |= ch_msk

            # Use autocorrelations to flag entire integrations which have
            # anomalous powers.  All antennas must agree for a integration
            # to get flagged.
            if not opts.flagmode.startswith('val'):
                new_mask = {}
                for p in data:
                  for k in data[p]:
                    i, j = aipy.miriad.bl2ij(k)
                    if i != j: continue
                    data_times = data[p][k].keys()
                    d = numpy.ma.array([data[p][k][t] for t in data_times],
                        mask=[mask[t] for t in data_times])
                    for i in numpy.where(flag_by_int(d))[0]:
                        t = data_times[i]
                        new_mask[t] = new_mask.get(t, 0) + 1
                for t in new_mask:
                    if new_mask[t] > 1: mask[t] |= 1
        if opts.outfile:
            for t in mask: mask[t] = list(mask[t])
            print 'Writing %f.xrfi' % jd
            f = open('%f.xrfi' % jd, 'w')
            pickle.dump(mask, f)
            f.close()
        else:
            # Generate a pipe for applying both the total_mask and the int_mask
            # to the data as it comes it.
            def rfi_mfunc(uv, preamble, data, flags):
                uvw, t, (i,j) = preamble
                return preamble, data, mask[t]

            uvi.rewind()
            uvo = aipy.miriad.UV(uvofile, status='new')
            uvo.init_from_uv(uvi)
            uvo.pipe(uvi, mfunc=rfi_mfunc, raw=True, append2hist='XRFI: ver %s, nsig %f, chans %s, mode %s, ch_thresh %f\n' %  (__version__, opts.nsig, opts.chans, opts.flagmode, opts.ch_thresh))
