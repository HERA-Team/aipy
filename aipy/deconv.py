"""
A module implementing various techniques for deconvolving an image by a
kernel.  Currently implemented are Clean, Least-Squares, and Maximum Entropy.

Author: Aaron Parsons
Date: 11/29/07
Revisions: None
"""

import numpy
# Find smallest representable # > 0 for setting clip level
clip_lev = numpy.finfo(numpy.float).tiny 

def recenter(a, c):
    """Slide the (0,0) point of matrix a to a new location tuple c.  This is
    useful for making an image centered on your screen after performing an
    inverse fft of uv data."""
    s = a.shape
    c = (c[0] % s[0], c[1] % s[1])
    a1 = numpy.concatenate([a[c[0]:], a[:c[0]]], axis=0)
    a2 = numpy.concatenate([a1[:,c[1]:], a1[:,:c[1]]], axis=1)
    return a2

def clean(im, ker, mdl=None, gain=.2, maxiter=10000, chkiter=100,verbose=False):
    """This is an implementation of the standard Hoegbom clean deconvolution 
    algorithm, which operates on the assumption that the image is composed of
    point sources.  This makes it a poor choice for images with distributed 
    flux.  The algorithm works by iteratively constructing a model.  Each 
    iteration, a point is added to this model at the location of the maximum 
    residual, with a fraction (specified by 'gain') of the magnitude.  The 
    convolution of that point is removed from the residual, and the process
    repeats.  Termination happens after 'maxiter' iterations, or when the
    clean loops starts increasing the magnitude of the residual.
    im: The image to be deconvolved.
    ker: The kernel to deconvolve by (must be same size as im).
    mdl: An a priori model of what the deconvolved image should look like.
    gain: The fraction of a residual used in each iteration.  If this is too
        low, clean takes unnecessarily long.  If it is too high, clean does
        a poor job of deconvolving.
    maxiter: The maximum number of iterations performed before terminating.
    chkiter: The number of iterations between when clean checks if the 
        residual is increasing.
    verbose: If true, prints some info on how things are progressing."""
    dim = im.shape[1]
    ker_pwr = ker.sum()
    G = ker_pwr / gain
    if mdl is None: mdl = numpy.zeros_like(im)
    # Get the starting residual
    dif = im - numpy.fft.ifft2(numpy.fft.fft(mdl) * numpy.fft.fft(ker)).real
    score = numpy.average(dif**2)
    prev_a = None
    n_mdl = mdl.copy()
    n_dif = dif.copy()
    mode = 0
    # Begin the clean loop
    for i in range(maxiter):
        # Rather than perform a convolution each loop, we'll just subtract
        # from the residual a scaled kernel centered on the point just added,
        # and to avoid recentering the kernel each time (because clean often
        # chooses the same point over and over), we'll buffer the previous one.
        a = numpy.argmax(n_dif)
        if a != prev_a:
            prev_a = a
            rec_ker = recenter(ker, (-a/dim, -a%dim))
        v = n_dif.flat[a] / G
        n_mdl.flat[a] += v
        n_dif -= v * rec_ker
        if i % chkiter == 0:
            # Check in on how clean is progressing.  Potenially terminate.
            n_score = numpy.average(abs(n_dif)**2)
            if verbose: print i, n_score
            if n_score > score:
                n_mdl, n_dif = mdl, dif
                break
            # If we're still doing fine, buffer this solution as the new best.
            score = n_score
            mdl = n_mdl.copy()
            dif = n_dif.copy()
    return n_mdl, n_dif

def lsq(im, ker, mdl=None, gain=.2, tol=1e-3, maxiter=1000, verbose=False):
    """Implements a simple least-square fitting procedure for deconvolving
    an image.  However, to save computing, the gradient of the fit at each 
    pixel with respect pixels in the image is approximated as being diagonal.
    In essence, this assumes that the convolution kernel is a delta-function.
    This assumption works for small kernels, but not so well for large ones.
    See Cornwell and Evans, 1984 "A Simple Maximum Entropy Deconvolution
    Algorithm" for more information about this approximation.  This
    deconvolution algorithm, unlike maximum entropy, makes no promises about
    maximizing smoothness while fitting to the expected noise and flux levels.
    That is, it can introduce structure for which there is no evidence in the
    original image.  Termination happens after 'maxiter' iterations, or
    when the score is changing by a fraction less than 'tol' between iterations.
    This makes the assumption that the true optimum has a smooth approach.
    im: The image to be deconvolved.
    ker: The kernel to deconvolve by (must be same size as im).
    mdl: An a priori model of what the deconvolved image should look like.
    gain: The fraction of the step size (calculated from the gradient) taken
        in each iteration.  If this is too low, the fit takes unnecessarily 
        long.  If it is too high, the fit process can oscillate.
    tol: See the termination information above.
    maxiter: The maximum number of iterations performed before terminating.
    verbose: If true, prints some info on how things are progressing."""
    if mdl is None: mdl = numpy.zeros_like(im)
    x = mdl.copy()
    # Estimate diagonal response of the "delta-function" kernel
    q = numpy.sqrt((ker**2).sum())
    # Buffer the FFT of the kernel so we don't have to recalculate it.
    inv_ker = numpy.fft.fft2(ker)
    # Function to calculate chi-square and gradient
    def f(x):
        x_conv_ker = numpy.fft.ifft2(numpy.fft.fft2(x) * inv_ker).real
        diff = im - x_conv_ker
        g_chi2 = -2*q*(diff)
        chi2 = diff**2
        return chi2, g_chi2
    prev_score = 0
    # Start the fit loop
    for i in range(maxiter):
        chi2, g_chi2 = f(x)
        score = numpy.average(chi2)
        term = abs(1 - prev_score/score)
        if verbose:
            slope = numpy.sqrt(numpy.average(g_chi2**2))
            print 'Step', i, ':', 'score', score, ', slope', slope, 
            print ', term', term
        # Terminate if change in score and slope is small
        if term < tol: break
        prev_score = score
        # For especially clean images, g_chi2 in some components can go to 0.
        # This check makes lsq a little slower for most images, though...
        d_x = numpy.where(abs(g_chi2) > 0, -(1/g_chi2) * chi2, 0)
        x += gain * d_x
        # Force positivity
        x = numpy.where(x < 0, 0, x)
    return x, (im -  numpy.fft.ifft2(numpy.fft.fft2(x) * inv_ker).real)

def maxent(im, ker, var0, mdl=None, gain=.1, tol=1e-3,
        maxiter=1000, verbose=False):
    """The maximum entropy deconvolution (MEM) (see Cornwell and Evans 1984
    "A Simple Maximum Entropy Deconvolution Algorithm" and Sault 1990
    "A Modification of the Cornwell and Evans Maximum Entropy Algorithm")
    is similar to least-squares deconvolution, but instead of simply
    minimizing the fit, maximum entropy seeks to do so only to within the
    specified variance var0, and then attempts to maximize the "smoothness"
    of the model.  This has several desirable effects, including uniqueness
    of solution, independence from flux distribution (all image scales are
    equally weighted in the model), and absence of spurious structure in the
    model.  Similar approximations are made in this implementation as in
    the least-squares implementation.
    im: The image to be deconvolved.
    ker: The kernel to deconvolve by (must be same size as im).
    var0: The estimated variance (noise power) in the image.  If none is
        provided, a quick lsq is used to estimate the variance of the residual.
    mdl: An a priori model of what the deconvolved image should look like.
        If none is provided, will default to a flat image with the same
        flux as the dirty image, normalized for the gain in the kernel.
    gain: The fraction of the step size (calculated from the gradient) taken
        in each iteration.  If this is too low, the fit takes unnecessarily 
        long.  If it is too high, the fit process can oscillate.
    tol: The termination criteria.  Lower = more optimized.  .01-.001 is normal.
    maxiter: The maximum number of iterations performed before terminating.
    verbose: If true, prints some info on how things are progressing."""
    d_i = im.flatten()
    q = numpy.sqrt((ker**2).sum())
    minus_two_q = -2*q
    two_q_sq = 2*q**2
    if mdl is None: mdl = numpy.ones_like(im) * numpy.average(im) / ker.sum()
    Nvar0 = d_i.size * var0
    inv_ker = numpy.fft.fft2(ker)
    m_i = mdl.flatten()
    def next_step(b_i, alpha, verbose=False):
        b_i.shape = im.shape
        b_i_conv_ker = numpy.fft.ifft2(numpy.fft.fft2(b_i) * inv_ker).real
        b_i.shape = m_i.shape
        diff = (im - b_i_conv_ker).flatten()
        chi2 = numpy.dot(diff,diff) - Nvar0
        g_chi2 = minus_two_q*(diff)
        g_J = (-numpy.log(b_i/m_i) - 1) - alpha * g_chi2
        gg_J = (-1/b_i) - alpha * two_q_sq
        # Define dot product using the metric of -gg_J^-1
        def dot(x, y): return (x*y/-gg_J).sum()
        score = dot(g_J, g_J) / dot(1,1)
        #score = abs(dot(g_J, g_J) / dot(1,1))
        d_alpha = (chi2 + dot(g_chi2,g_J)) / dot(g_chi2,g_chi2)
        # For especially clean images, gg_J in some components can go to 0.
        # This check makes lsq a little slower for most images, though...
        d_b_i = numpy.where(abs(gg_J) > 0, -1/gg_J * (g_J - d_alpha*g_chi2), 0)
        if verbose:
            print '    score', score, 'fit', numpy.dot(diff,diff)
            print '    alpha', alpha, 'd_alpha', d_alpha
        return d_b_i, d_alpha, score
    alpha = 0.
    b_i = m_i.copy()
    info = {'success':True, 'term':'maxiter', 'var0':var0, 'tol':tol}
    for i in range(maxiter):
        if verbose: print 'Step', i, ':'
        d_b_i, d_alpha, score = next_step(b_i, alpha, verbose=verbose)
        if score < tol and score > 0:
            info['term'] = 'tol'
            break
        elif score > 1e10 or numpy.isnan(score) or score <= 0:
            info.update({'term':'divergence', 'success':False})
            break
        b_i += gain * d_b_i
        alpha += gain * d_alpha
        b_i = numpy.where(b_i < clip_lev, clip_lev, b_i)
    b_i.shape = im.shape
    info.update({'res':im - numpy.fft.ifft2(numpy.fft.fft2(b_i) * inv_ker).real,
        'score': score, 'alpha': alpha, 'iter':i+1})
    return b_i, info

