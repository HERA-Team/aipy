#! /usr/bin/env python
"""
A script for fitting the delays, positions, and offsets of antennas given
starting parameters in "aipy.loc" and a list of sources (from "aipy.src").

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, os, optparse

def calc_gain(phsdata):
    img = n.fft.ifft(phsdata.filled(0))
    return img * fit_window(img.shape[0])

def fit_window(width):
    w = n.arange(width)
    w = n.where(w > width/2, width - w, w)
    return 10.**(-w)

o = optparse.OptionParser()
o.set_usage('fit_dly_pos.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, loc=True, 
    src=True, dec=True)
o.add_option('--fitants', dest='fitants',
    help='List of antennas to fit parameters for.')
o.add_option('--swath', dest='swath', default=0, type='int',
    help='Number of bins around center to wipe out when removing a src.')
o.add_option('--prms', dest='prms', 
    help='Comma delimited list of paramters to fit (x,y,z,delay,offset).')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
a.scripting.uv_selector(uv, opts.ant, opts.pol)
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
cat = a.scripting.parse_src(opts.src, force_cat=True)
fitants = map(int, opts.fitants.split(','))
prms = opts.prms.split(',')
del(uv)

# Generate list of parameters to fit
prm_dict = {}
for a in fitants: prm_dict[a] = prms
print prm_dict
start_prms = aa.get_params(prm_dict)
prm_list, key_list = a.fit.flatten_prms(start_prms)

first_fit = None    # Used to normalize fit values to the starting fit

def fit_func(prms):
    global first_fit
    prms = a.fit.reconstruct_prms(prms, key_list)
    a.fit.print_params(prms)
    aa.set_params(prms)
    score = 0
    cnt,curtime = 0,None
    src_phs = {}
    for uvfile in args:
        sys.stdout.write('.'), ; sys.stdout.flush()
        uv = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            # Use only every Nth integration, if decimation is specified.
            if curtime != t:
                curtime = t
                cnt = (cnt + 1) % opts.decimate
                for s in src_phs:
                    score += abs(src_phs[s]).sum()
                    src_phs = {}
                if cnt == 0:
                    aa.set_jultime(t)
                    cat.compute(aa)
            if cnt != 0: continue
            for s in cat:
                ds = n.where(f, 0, d)
                # Remove all other interfering sources from data
                for rs in cat:
                    if rs != s: ds = aa.rmsrc(ds, cat[rs], i, j)
                ds = n.where(f, 0, d)
                try:
                    ds = aa.phs2src(ds, cat[s], i, j)
                    src_phs[s] = src_phs.get(s,0) + calc_gain(ds)
                except(a.ant.PointingError): pass
    score = 1 / score
    if first_fit is None: first_fit = score
    print 'Score:', score, 
    print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
    print '-------------------------------------------------------------------'
    return score / first_fit

a.optimize.fmin(fit_func, prm_list)
print 'Done!'

