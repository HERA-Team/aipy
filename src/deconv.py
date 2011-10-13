"""
A module implementing various techniques for deconvolving an image by a
kernel.  Currently implemented are Clean, Least-Squares, Maximum Entropy,
and Annealing.  Standard parameters to these functions are:
im = image to be deconvolved.
ker = kernel to deconvolve by (must be same size as im).
mdl = a priori model of what the deconvolved image should look like.
maxiter = maximum number of iterations performed before terminating.
tol = termination criterion, lower being more optimized.
verbose =  print info on how things are progressing.
lower = lower bound of pixel values in deconvolved image
upper = upper bound of pixel values in deconvolved image
"""

import numpy as n, sys, _deconv

# Find smallest representable # > 0 for setting clip level
lo_clip_lev = n.finfo(n.float).tiny 

def clean(im, ker, mdl=None, area=None, gain=.1, maxiter=10000, tol=1e-3, 
        stop_if_div=True, verbose=False):
    """This standard Hoegbom clean deconvolution algorithm operates on the 
    assumption that the image is composed of point sources.  This makes it a 
    poor choice for images with distributed flux.  In each iteration, a point 
    is added to the model at the location of the maximum residual, with a 
    fraction (specified by 'gain') of the magnitude.  The convolution of that 
    point is removed from the residual, and the process repeats.  Termination 
    happens after 'maxiter' iterations, or when the clean loops starts 
    increasing the magnitude of the residual.  This implementation can handle 
    1 and 2 dimensional data that is real valued or complex.
    gain: The fraction of a residual used in each iteration.  If this is too
        low, clean takes unnecessarily long.  If it is too high, clean does
        a poor job of deconvolving."""
    if mdl is None:
        #mdl = n.zeros_like(im)
        mdl = n.zeros(im.shape, dtype=im.dtype)
        res = im.copy()
    else:
        mdl = mdl.copy()
        if len(mdl.shape) == 1:
            res = im - n.fft.ifft(n.fft.fft(mdl) * \
                                  n.fft.fft(ker)).astype(im.dtype)
        elif len(mdl.shape) == 2:
            res = im - n.fft.ifft2(n.fft.fft2(mdl) * \
                                   n.fft.fft2(ker)).astype(im.dtype)
        else: raise ValueError('Number of dimensions != 1 or 2')
    if area is None:
        area = n.ones(im.shape, dtype=n.int)
    else:
        area = area.astype(n.int)
        
    iter = _deconv.clean(res, ker, mdl, area,
            gain=gain, maxiter=maxiter, tol=tol, 
            stop_if_div=int(stop_if_div), verbose=int(verbose))
    score = n.sqrt(n.average(n.abs(res)**2))
    info = {'success':iter > 0 and iter < maxiter, 'tol':tol}
    if iter < 0: info.update({'term':'divergence', 'iter':-iter})
    elif iter < maxiter: info.update({'term':'tol', 'iter':iter})
    else: info.update({'term':'maxiter', 'iter':iter})
    info.update({'res':res, 'score':score})
    if verbose:
        print 'Term Condition:', info['term']
        print 'Iterations:', info['iter']
        print 'Score:', info['score']
    return mdl, info

def recenter(a, c):
    """Slide the (0,0) point of matrix a to a new location tuple c."""
    s = a.shape
    c = (c[0] % s[0], c[1] % s[1])
    a1 = n.concatenate([a[c[0]:], a[:c[0]]], axis=0)
    a2 = n.concatenate([a1[:,c[1]:], a1[:,:c[1]]], axis=1)
    return a2

#def clean(im, ker, mdl=None, gain=.2, maxiter=10000, chkiter=100, tol=1e-3,
#        lower=lo_clip_lev, upper=n.Inf, verbose=False):
#    """This standard Hoegbom clean deconvolution algorithm operates on the 
#    assumption that the image is composed of point sources.  This makes it a 
#    poor choice for images with distributed flux.  In each iteration, a point 
#    is added to the model at the location of the maximum residual, with a 
#    fraction (specified by 'gain') of the magnitude.  The convolution of that 
#    point is removed from the residual, and the process repeats.  Termination 
#    happens after 'maxiter' iterations, or when the clean loops starts 
#    increasing the magnitude of the residual.
#    gain: The fraction of a residual used in each iteration.  If this is too
#        low, clean takes unnecessarily long.  If it is too high, clean does
#        a poor job of deconvolving.
#    chkiter: The number of iterations between when clean checks if the 
#        residual is increasing."""
#    dim = im.shape[1]
#    q = n.sqrt((ker**2).sum())
#    G = q / gain
#    if mdl is None: mdl = n.zeros_like(im)
#    # Get the starting residual
#    inv_ker = n.fft.fft2(ker)
#    dif = im - n.fft.ifft2(n.fft.fft2(mdl) * inv_ker).real
#    score = n.sqrt(n.average(dif**2))
#    a0 = None
#    n_mdl, n_dif = mdl.copy(), dif.copy()
#    info = {'success':True, 'term':'maxiter', 'tol':tol}
#    # Begin the clean loop
#    for i in range(maxiter):
#        # Rather than perform a convolution each loop, we'll just subtract
#        # from the residual a scaled kernel centered on the point just added,
#        # and to avoid recentering the kernel each time (because clean often
#        # chooses the same point over and over), we'll buffer the previous one.
#        a = n.argmax(n_dif)
#        if a != a0:
#            a0 = a
#            rec_ker = recenter(ker, (-a/dim, -a%dim))
#        v = n_dif.flat[a] / G
#        n_mdl.flat[a] += v
#        n_dif -= v * rec_ker
#        # Check in on how clean is progressing.  Potenially terminate.
#        if i % chkiter == 0:
#            # Mystery: why does clean to worse when exact dif is computed?
#            #n_mdl = n_mdl.clip(lower, upper)
#            #n_dif = im - n.fft.ifft2(n.fft.fft2(mdl) * inv_ker).real
#            n_score = n.sqrt(n.average(n_dif**2))
#            if verbose:
#                print 'Step %d:' % i, 'score %f,' % n_score, 'best %f' % score
#            if n_score - score > score * tol:
#                n_mdl, n_dif = mdl, dif
#                info['term'] = 'tol'
#                break
#            score,mdl,dif = n_score, n_mdl.copy(), n_dif.copy()
#    info.update({'res':n_dif, 'score': score, 'iter':i+1})
#    return n_mdl, info

def lsq(im, ker, mdl=None, gain=.1, tol=1e-3, maxiter=200, 
        lower=lo_clip_lev, upper=n.Inf, verbose=False):
    """This simple least-square fitting procedure for deconvolving an image 
    saves computing by assuming a diagonal pixel-pixel gradient of the fit.
    In essence, this assumes that the convolution kernel is a delta-function.
    This works for small kernels, but not so well for large ones.  See Cornwell 
    and Evans, 1984 "A Simple Maximum Entropy Deconvolution Algorithm" for more 
    information about this approximation.  Unlike maximum entropy, lsq makes 
    no promises about maximizing smoothness, but needs no information
    about noise levels.  Structure can be introduced for which there is no 
    evidence in the original image.  Termination happens when the fractional 
    score change is less than 'tol' between iterations.
    gain: The fraction of the step size (calculated from the gradient) taken
        in each iteration.  If this is too low, the fit takes unnecessarily 
        long.  If it is too high, the fit process can oscillate."""
    if mdl is None:
        #mdl = n.zeros_like(im)
        mdl = n.zeros(im.shape, dtype=im.dtype)
    x = mdl.copy()
    # Estimate gain of the kernel
    q = n.sqrt((ker**2).sum())
    ker_i = n.fft.fft2(ker)
    info = {'success':True, 'term':'maxiter', 'tol':tol}
    # Function to calculate chi-square and gradient
    def f(x):
        x_conv_ker = n.fft.ifft2(n.fft.fft2(x) * ker_i).real
        diff = im - x_conv_ker
        g_chi2 = -2*q*(diff)
        chi2 = diff**2
        return chi2, g_chi2
    score = 0
    # Start the fit loop
    for i in range(maxiter):
        chi2, g_chi2 = f(x)
        n_score = n.average(chi2)
        term = abs(1 - score/n_score)
        if verbose:
            slope = n.sqrt(n.average(g_chi2**2))
            print 'Step %d:' % i, 'score',  score, 
            print 'slope', slope, 'term', term
        if term < tol:
            info['term'] = 'tol'
            break
        score = n_score
        # For especially clean imgs, g_chi2 in some components can go to 0.
        # This check makes lsq a little slower for most images, though...
        d_x = n.where(abs(g_chi2) > 0, -(1/g_chi2) * chi2, 0)
        x = n.clip(x + gain * d_x, lower, upper)
    info.update({'res':im - n.fft.ifft2(n.fft.fft2(x) * ker_i).real,
        'score': score, 'iter':i+1})
    return x, info

def maxent(im, ker, var0, mdl=None, gain=.1, tol=1e-3, maxiter=200, 
        lower=lo_clip_lev, upper=n.Inf, verbose=False):
    """Maximum entropy deconvolution (MEM) (see Cornwell and Evans 1984
    "A Simple Maximum Entropy Deconvolution Algorithm" and Sault 1990
    "A Modification of the Cornwell and Evans Maximum Entropy Algorithm")
    is similar to lsq, but the fit is only optimized to within the specified
    variance (var0) and then "smoothness" is maximized.  This has several 
    desirable effects including uniqueness of solution, equal weighting of
    Fourier components, and absence of spurious structure.  The same 
    delta-kernel approximation (see lsq) is made here.
    var0: The estimated variance (noise power) in the image.  If none is
        provided, a quick lsq is used to estimate the variance of the residual.
    gain: The fraction of the step size (calculated from the gradient) taken
        in each iteration.  If this is too low, the fit takes unnecessarily 
        long.  If it is too high, the fit process can oscillate."""
    d_i = im.flatten()
    q = n.sqrt((ker**2).sum())
    minus_two_q = -2*q
    two_q_sq = 2*q**2
    if mdl is None:
         #mdl = n.ones_like(im) * n.average(im) / ker.sum() 
         mdl = n.ones(im.shape, dtype=im.dtype) * n.average(im) / ker.sum() 
    #if mdl is None: mdl = n.ones_like(im) * n.average(im) / q
    Nvar0 = d_i.size * var0
    inv_ker = n.fft.fft2(ker)
    m_i = mdl.flatten()
    def next_step(b_i, alpha, verbose=False):
        b_i.shape = im.shape
        b_i_conv_ker = n.fft.ifft2(n.fft.fft2(b_i) * inv_ker).real
        b_i.shape = m_i.shape
        diff = (im - b_i_conv_ker).flatten()
        chi2 = n.dot(diff,diff) - Nvar0
        g_chi2 = minus_two_q*(diff)
        g_J = (-n.log(b_i/m_i) - 1) - alpha * g_chi2
        gg_J = (-1/b_i) - alpha * two_q_sq
        # Define dot product using the metric of -gg_J^-1
        def dot(x, y): return (x*y/-gg_J).sum()
        score = dot(g_J, g_J) / dot(1,1)
        d_alpha = (chi2 + dot(g_chi2,g_J)) / dot(g_chi2,g_chi2)
        # For especially clean images, gg_J in some components can go to 0.
        # This check makes lsq a little slower for most images, though...
        d_b_i = n.where(abs(gg_J) > 0, -1/gg_J * (g_J - d_alpha*g_chi2), 0)
        if verbose:
            print '    score', score, 'fit', n.dot(diff,diff)
            print '    alpha', alpha, 'd_alpha', d_alpha
        return d_b_i, d_alpha, score
    alpha = 0.
    b_i = m_i.copy()
    info = {'success':True, 'term':'maxiter', 'var0':var0, 'tol':tol}
    for i in range(maxiter):
        if verbose: print 'Step %d:' % i
        d_b_i, d_alpha, score = next_step(b_i, alpha, verbose=verbose)
        if score < tol and score > 0:
            info['term'] = 'tol'
            break
        elif score > 1e10 or n.isnan(score) or score <= 0:
            info.update({'term':'divergence', 'success':False})
            break
        b_i = n.clip(b_i + gain * d_b_i, lower, upper)
        alpha += gain * d_alpha
    b_i.shape = im.shape
    info.update({'res':im - n.fft.ifft2(n.fft.fft2(b_i) * inv_ker).real,
        'score': score, 'alpha': alpha, 'iter':i+1})
    return b_i, info

def maxent_findvar(im, ker, var=None, f_var0=.6, mdl=None, gain=.1, tol=1e-3, 
        maxiter=200, lower=lo_clip_lev, upper=n.Inf, verbose=False, 
        maxiterok=False):
    """This frontend to maxent tries to find a variance for which maxent will
    converge.  If the starting variance (var) is not specified, it will be
    estimated as a fraction (f_var0) of the variance of the residual of a 
    lsq deconvolution, and then a search algorithm tests an ever-widening
    range around that value.  This function will search until it succeeds."""
    cl, info, cnt = None, None, -1
    if var is None:
        # Get a starting estimate of variance to use via residual of lsq
        junk, info = lsq(im, ker, mdl=mdl, gain=gain, tol=tol,
            maxiter=maxiter/4, lower=lower, upper=upper, verbose=False)
        var = n.var(info['res'])
        if verbose: print 'Using', f_var0, 'of LSQ estimate of var=', var
        var *= f_var0
    else:
        if verbose: print 'Using specified var=', var
    while cl is None:
        print cnt
        if cnt == -1: v = var
        else: v = var / (1.5**cnt)
        while cnt < 0 or v < var * (1.5**cnt):
            if verbose:
                print 'Trying var=', v,
                sys.stdout.flush()
            c, i = maxent(im, ker, v, mdl=mdl, gain=gain, tol=tol,
                maxiter=maxiter, lower=lower, upper=upper, verbose=False)
            if verbose:
                print 'success %d,' % i['success'],
                print 'term: %s,' % i['term'], 'score:' , i['score']
            # Check if fit converged
            if i['success'] and (maxiterok or i['term'] == 'tol'):
                cl, info = c, i
                break
            else:
                if not cl is None or cnt == -1: break
                v *= 1.2 ** (1./(2*(cnt+1)))
        cnt += 1
    if verbose: print 'Done with MEM.'
    return cl, info

def anneal(im, ker, mdl=None, maxiter=1000, lower=lo_clip_lev, upper=n.Inf,
        cooling=lambda i,x: 1e+1*(1-n.cos(i/50.))*(x**2), verbose=False):
    """Annealing takes a non-deterministic approach to deconvolution by
    randomly perturbing the model and selecting perturbations that improve the 
    residual.  By slowly reducing the temperature of the perturbations, 
    annealing attempts to settle into a global minimum.  Annealing is slower
    than lsq for a known gradient, but is less sensitive to gradient errors 
    (it can solve for wider kernels). Faster cooling speeds terminate more 
    quickly, but are less likely to find the global minimum.  This 
    implementation assigns a temperature to each pixel proportional to the 
    magnitude of the residual in that pixel and the global cooling speed.
    cooling: A function accepting (iteration,residuals) that returns a 
        vector of standard deviation for noise in the respective pixels.
        Picking the scaling of this function correctly is vital for annealing
        to work."""
    if mdl is None: mdl = n.zeros_like(im)
    q = n.sqrt((ker**2).sum())
    inv_ker = n.fft.fft2(ker)
    dif = im - n.fft.ifft2(n.fft.fft2(mdl) * inv_ker).real
    score = n.average(dif**2)
    #info = {'success':True, 'term':'maxiter', 'speed':speed}
    info = {'success':True, 'term':'maxiter'}
    for i in range(maxiter):
        delta = n.random.normal(scale=1., size=mdl.shape) * cooling(i, dif/q)
        n_mdl = n.clip(mdl + delta, lower, upper)
        n_dif = im - n.fft.ifft2(n.fft.fft2(n_mdl) * inv_ker).real
        n_score = n.average(n_dif**2)
        if verbose: print 'Step %d:' % i, n_score, score
        if n_score < score: mdl, dif, score = n_mdl, n_dif, n_score
    info.update({'res':dif, 'score': score, 'iter':i+1})
    return mdl, info
