"""
A pure-python, spline-like implementation of 1D interpolation.  Uses
an FIR filter (spline) to interpolate, and a polynomial to extrapolate
boundaries.

Author: Aaron Parsons
Date: 12/12/07
Revisions:
"""

__all__ = ['interpolate', 'default_filter']

import numpy as n

def subsample(a, factor):
    rv = n.zeros((a.size, factor), dtype=a.dtype)
    rv[:,0] = a
    wgts = n.zeros((a.size, factor), dtype=n.float)
    wgts[:,0] = 1.
    return rv.flatten(), wgts.flatten()

def polyextend(y, nsamples, degree=2):
    """Extend y on either end by nsamples with a polynomial of the specified
    degree.  Assumes y is sampled on a uniform x grid."""
    x = n.arange(-nsamples,nsamples)
    p = n.polyfit(x[-nsamples:],y[:nsamples],degree)
    y0 = n.polyval(p, x[:nsamples])
    x = n.arange(-nsamples+1,nsamples+1)
    p = n.polyfit(x[:nsamples],y[-nsamples:],degree)
    y1 = n.polyval(p, x[-nsamples:])
    return n.concatenate([y0, y, y1])
    
def default_filter(xs, freq=.25):
    """A basic smoothing filter using a Hamming window and a sinc function
    of the specified frequency.  Input a sample grid [-x,x] to get the
    FIR coefficients used for smoothing."""
    i = n.arange(xs.size, dtype=n.float) / xs.size * 2*n.pi
    return (1-n.cos(i)) * n.sinc(freq*xs)

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
    fir = filter(n.arange(-order,order+step,step))
    ys_ss = n.convolve(ys_ss, fir, mode='valid')
    wgts = n.convolve(wgts, fir, mode='valid')
    return ys_ss / wgts
   
