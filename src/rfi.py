"""
Module for detecting and flaging RFI related effects.
"""

import numpy as np, optimize

def gaussian(amp, sig, off, x):
    """Generate gaussian value at x given amplitude, sigma, and x offset."""
    return amp * np.exp(-((x - off)/sig)**2)

def fit_gaussian(xs, ys):
    """Fit a gaussian (returning amplitude, sigma, offset) to a set of x
    values and their corresponding y values."""
    def fitfunc(args):
        amp, sig, off = args
        sum = ((ys - gaussian(amp, sig, off, xs))**2).sum()
        return sum
    return optimize.fmin(fitfunc, (ys.sum(), 50, ys.argmax()), disp=0)

def gen_rfi_thresh(data, nsig=2, cnt_per_bin=1000):
    """Generate a threshold at nsig times the standard deviation (above,below)
    the mean, outside of which data is considered rfi."""
    try: data = data.compressed()
    except(AttributeError): pass
    data = np.log10(np.abs(data))
    try: h, bvals = np.histogram(data, bins=data.size/cnt_per_bin)
    except(ValueError): return None, None
    if h.size == 0: return None, None
    # Fit a gaussian to histogram (better than just std-dev of data)
    amp, sig, off = fit_gaussian(np.arange(h.size), h)
    sig = abs(sig)
    hi_thresh = np.clip(np.round(off + nsig*sig), 0, len(bvals)-1)
    lo_thresh = np.clip(np.round(off - nsig*sig), 0, len(bvals)-1)
    return 10**bvals[hi_thresh], 10**bvals[lo_thresh]

def flag_by_int(preflagged_auto, nsig=1, raw=False):
    """Flag rfi for an autocorrelation.  Iteratively removes outliers and
    fits a smooth passband, then uses smooth passband to remove outliers
    (both positive and negative).  If 'raw' is specified, no smooth function
    is removed."""
    pwr_vs_t = np.ma.average(np.abs(preflagged_auto), axis=1)
    mask = pwr_vs_t.mask
    if raw:
        spikey_pwr_vs_t = pwr_vs_t
    else: 
        spikey_pwr_vs_t = np.abs(pwr_vs_t - remove_spikes(pwr_vs_t, mask.copy()))
    hi_thr, lo_thr = gen_rfi_thresh(spikey_pwr_vs_t, cnt_per_bin=20, nsig=nsig)
    spikey_pwr_vs_t = spikey_pwr_vs_t.filled(hi_thr)
    mask = spikey_pwr_vs_t >= hi_thr
    return mask.astype(np.int)

def remove_spikes(data, mask=None, order=6, iter=3, return_poly=False):
    """Iteratively fits a smooth function by removing outliers and fitting
    a polynomial, then using the polynomial to remove other outliers."""
    xs = np.arange(data.size)
    if mask is None or mask.size != data.size:
        mask = np.zeros(data.shape, dtype=np.bool)
    im = np.logical_not(mask)
    nxs = xs.compress(im)
    ndata = data.compress(im)
    if len(nxs) != 0: p = np.polyfit(nxs, ndata, deg=order)
    else: p = np.polyfit(xs, data, deg=order)
    if iter != 0:
        residual = np.abs(ndata - np.polyval(p, nxs))
        sig = np.sqrt(np.average(residual**2))
        mask.put(np.where(residual > sig), 1)
        p = remove_spikes(data, mask=mask, order=order,
            iter=iter-1, return_poly=True)
    if return_poly: return p
    else: return np.polyval(p, xs)

