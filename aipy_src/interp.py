"""
A pure-python, spline-like implementation of 1D interpolation.  Uses
an FIR filter (spline) to interpolate, and a polynomial to extrapolate
boundaries.

Author: Aaron Parsons
Date: 12/12/07
Revisions:
"""

__all__ = ['interpolate', 'default_filter']

import numpy as np

def subsample(a, factor):
    rv = np.zeros((a.size, factor), dtype=a.dtype)
    rv[:,0] = a
    wgts = np.zeros((a.size, factor), dtype=np.float)
    wgts[:,0] = 1.
    return rv.flatten(), wgts.flatten()

def polyextend(y, nsamples, degree=2):
    """Extend y on either end by nsamples with a polynomial of the specified
    degree.  Assumes y is sampled on a uniform x grid."""
    x = np.arange(-nsamples,nsamples)
    p = np.polyfit(x[-nsamples:],y[:nsamples],degree)
    y0 = np.polyval(p, x[:nsamples])
    x = np.arange(-nsamples+1,nsamples+1)
    p = np.polyfit(x[:nsamples],y[-nsamples:],degree)
    y1 = np.polyval(p, x[-nsamples:])
    return np.concatenate([y0, y, y1])
    
def default_filter(xs, freq=.25):
    """A basic smoothing filter using a Hamming window and a sinc function
    of the specified frequency.  Input a sample grid [-x,x] to get the
    FIR coefficients used for smoothing."""
    i = np.arange(xs.size, dtype=np.float) / xs.size * 2*np.pi
    return (1-np.cos(i)) * np.sinc(freq*xs)

def interpolate(ys, factor, filter=default_filter, order=4):
    """Oversample ys by the specified factor using a filter function to
    interpolated between samples.  The filter function will be used to
    construct an FIR filter for x values [-order,order] in steps of 1/order.
    Order should be even to make this an averaging filter.  Internally,
    ys is extended by 2nd degree polynomial in both directions to attempt
    a smooth boundary transition."""
    step = 1./factor
    ys = polyextend(ys, order)
    ys_ss, wgts = subsample(ys, factor)
    fir = filter(np.arange(-order,order+step,step))
    ys_ss = np.convolve(ys_ss, fir, mode='valid')
    wgts = np.convolve(wgts, fir, mode='valid')
    return ys_ss / wgts
   
