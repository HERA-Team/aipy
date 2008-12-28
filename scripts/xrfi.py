#!/usr/bin/python
"""
A script for detecting and flagging RFI related effects in UV files.

Author: Aaron Parsons
Date: 5/31/07
Revisions: None
"""

import numpy, aipy.miriad, scipy.optimize

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
    return scipy.optimize.fmin(fitfunc, (ys.sum(), 50, ys.argmax()), disp=0)

def gen_rfi_thresh(data, nsig=2, per_bin=1000):
    """Generate a threshold at nsig times the standard deviation above the mean
    above which data is considered rfi."""
    if numpy.ma.isMA(data): histdata = numpy.log10(abs(data.filled(0))+1e-6)
    else: histdata = numpy.log10(abs(data)+1e-6)
    h, bvals = numpy.histogram(histdata, bins=histdata.size/per_bin)
    amp, sig, off = fit_gaussian(numpy.arange(h.size), h)
    hi_thresh = numpy.clip(numpy.round(off + nsig*sig), 0, len(bvals)-1)
    lo_thresh = numpy.clip(numpy.round(off - nsig*sig), 0, len(bvals)-1)
    hi_thresh = bvals[hi_thresh]
    lo_thresh = bvals[lo_thresh]
    return 10**hi_thresh, 10**lo_thresh

def flag_by_int(preflagged_auto):
    """Flag rfi for an autocorrelation.  Iteratively removes outliers and
    fits a smooth passband, then uses smooth passband to remove outliers
    (both positive and negative)."""
    pwr_vs_t = numpy.ma.average(abs(preflagged_auto), axis=1)
    spikey_pwr_vs_t = abs(pwr_vs_t - remove_spikes(pwr_vs_t))
    hi_thr, lo_thr = gen_rfi_thresh(spikey_pwr_vs_t, per_bin=20, nsig=1)
    mask = spikey_pwr_vs_t > hi_thr
    return mask

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
    import os, sys
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('xrfi.py [options] *.uv')
    p.add_option('-s', '--sky_removal', dest='sky_removal', action='store_true',
        help='Window to filter out sky-based sources before flagging rfi.  Applies Parzen window to lag-domain data.')
    p.set_description(__doc__)

    opts, args = p.parse_args(sys.argv[1:])

    for uvfile in args:
        print 'Working on', uvfile
        uvofile = uvfile+'r'
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
        uvi = aipy.miriad.UV(uvfile)

        # Gather all data and each time step
        window = None
        data = {}
        times = []
        while True:
            p, d = uvi.read_data()
            if d.size == 0: break
            i, j = aipy.miriad.bl2ij(p[-1])
            if i != j and opts.sky_removal:
                if window is None: 
                    window = d.size/2 - abs(numpy.arange(d.size) - d.size/2)
                img = numpy.fft.ifft(d)
                img *= window
                d = numpy.fft.fft(img)
            d.shape = (1,) + d.shape
            try: data[p[-1]].append(d)
            except(KeyError): data[p[-1]] = [d]
            if len(times) == 0 or times[-1] != p[-2]: times.append(p[-2])

        #import pylab
        #ks = data.keys()
        #ks.sort()
        #for n, k in enumerate(ks):
        #    d = numpy.ma.concatenate(data[k], axis=0)
        #    pylab.subplot(4, 3, n+1)
        #    pylab.imshow(numpy.log10(numpy.abs(d)+1e-6))
        #    pylab.title(str(aipy.miriad.bl2ij(k)))
        #pylab.show()

        # Generate a single mask for all baselines which masks if any
        # baseline has an outlier at that freq/time.  Data flagged
        # strictly on basis of nsigma above mean.
        total_mask = 0
        for k in data:
            data[k] = numpy.ma.concatenate(data[k], axis=0)
            i, j = aipy.miriad.bl2ij(k)
            if i == j: continue
            hi_thresh, lo_thresh = gen_rfi_thresh(data[k])
            total_mask |= numpy.where(abs(data[k]) > hi_thresh, 1, 0)

        # Use autocorrelations to flag entire integrations which have
        # anomalous powers.  All antennas must agree for a integration
        # to get flagged.
        int_mask = 1
        for k in data:
            data[k] = numpy.ma.array(data[k], mask=total_mask)
            i, j = aipy.miriad.bl2ij(k)
            if i == j: int_mask &= flag_by_int(data[k])
        for n in numpy.where(int_mask): total_mask[n] = 1

        # Generate a pipe for applying both the total_mask and the int_mask
        # to the data as it comes it.
        uvi.rewind()
        uvo = aipy.miriad.UV(uvofile, status='new')
        def rfi_mfunc(uv, preamble, data):
            bl = preamble[-1]
            i, j = aipy.miriad.bl2ij(bl)
            t = times.index(preamble[-2])
            if i == j: mask = numpy.logical_or(total_mask[t], int_mask[t])
            else: mask = total_mask[t]
            data = numpy.ma.array(data.data, mask=mask)
            return preamble, data

        # Apply the pipe to the data
        aipy.miriad.pipe_uv(uvi, uvo, mfunc=rfi_mfunc)
        del(uvi)
        del(uvo)
