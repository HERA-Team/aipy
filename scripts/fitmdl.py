#! /usr/bin/env python
"""
A script for fitting the amplitudes and passbands of antennas given
starting parameters in "a.loc" and a list of sources (from "a.src").

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('fitmdl.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True,
    loc=True, src=True, dec=True)
o.add_option('--aprms', dest='aprms',
    help='Comma delimited list of paramters to fit independently for each antenna.')
o.add_option('--fitants', dest='fitants', default='all',
    help='Comma delimited list of antennas to fit. Default "all".')
o.add_option('--aaprms', dest='aaprms',
    help='Comma delimited list of AntennaArray (not Antenna) parameters to fit.  Unless added specifically in loc file, default AntennaArray has no parameters to fit.')
o.add_option('--shprms', dest='shprms',
    help='Comma delimited list of parameters to fit whose values will be shared for all antennas.')
o.add_option('--sprms', dest='sprms', 
    help='Source=param pairs for fitting source parameters.')
o.add_option('--sim_autos', dest='sim_autos', action='store_true',
    help='Use auto-correlations in fitting.  Default is to use only cross-correlations.')
o.add_option('-n', '--norm', dest='norm', action='store_true',
    help='Normalize residual to data strength.  Puts weaker sources on similar footing as strong sources.')
o.add_option('--maxiter', dest='maxiter', type='float', default=-1,
    help='Maximum # of iterations to run.  Default is infinite.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
a.scripting.uv_selector(uv, opts.ant, opts.pol)
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
cat = a.scripting.parse_srcs(opts.src, force_cat=True)
cat.set_params(a.loc.get_src_prms(opts.loc))
(uvw,t,(i,j)),d = uv.read()
aa.set_jultime(t)
cat.compute(aa)
if opts.fitants.startswith('all'): fitants = range(uv['nants'])
else: fitants = map(int, opts.fitants.split(','))
del(uv)

# Generate list of parameters to fit
if opts.aprms is None: aprms = []
else: aprms = opts.aprms.split(',')
if opts.shprms is None: shprms = []
else: shprms = opts.shprms.split(',')
if opts.aaprms is None: aaprms = []
else: aaprms = opts.aaprms.split(',')
prm_dict = {}
for ant in fitants: prm_dict[ant] = aprms[:]
prm_dict[fitants[0]].extend(shprms)
prm_dict['aa'] = aaprms
print prm_dict
start_prms = aa.get_params(prm_dict)
if opts.sprms != None:
    sprms = [w.split('=') for w in opts.sprms.split(',')]
    sprm_dict = {}
    for src,prm in sprms:
        if sprm_dict.has_key(src): sprm_dict[src].append(prm)
        else: sprm_dict[src] = [prm]
    start_prms.update(cat.get_params(sprm_dict))
else: sprms = []
prm_list, key_list = a.fit.flatten_prms(start_prms)

first_fit = None    # Used to normalize fit values to the starting fit
mfq = cat.get_mfreqs()

def fit_func(prms):
    global first_fit
    prms = a.fit.reconstruct_prms(prms, key_list)
    a.fit.print_params(prms)
    for ant in fitants:
        for prm in shprms:
            prms[ant][prm] = prms[fitants[0]][prm]
    aa.set_params(prms)
    cat.set_params(prms)
    a1,a2,th = cat.get_srcshapes()
    #if n.all(asz == 0): asz = None  # Making None bypasses a computation step
    score,nsamples = 0.,0.
    cnt,curtime = 0,None
    for uvfile in args:
        sys.stdout.write('.') ; sys.stdout.flush()
        uv = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            # Use only every Nth integration, if decimation is specified.
            if curtime != t:
                curtime = t
                cnt = (cnt + 1) % opts.decimate
                if cnt == 0:
                    aa.set_jultime(t)
                    cat.compute(aa)
                    eqs = cat.get_crds('eq', ncrd=3)
                    flx = cat.get_fluxes()
                    ind = cat.get_indices()
                    aa.sim_cache(eqs, flx, indices=ind,
                        mfreqs=mfq, srcshapes=(a1,a2,th))
            if cnt != 0: continue
            if not opts.sim_autos and i == j: continue
            d = d.take(chans)
            f = f.take(chans)
            sim_d = aa.sim(i, j, pol=a.miriad.pol2str[uv['pol']])
            sq,difsq = n.abs(d)**2, n.abs(d - sim_d)**2
            difsq = n.where(f, 0, difsq)
            if opts.norm:
                sq = n.ma.masked_equal(sq, 0)
                score += (difsq/sq).sum()
                nsamples += n.logical_not(f).sum()
            else:
                score += difsq.sum()
                nsamples += sq.sum()
    score = n.sqrt(score / nsamples)
    print
    if first_fit is None: first_fit = score
    print 'Score:', score, 
    print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
    print '-------------------------------------------------------------------'
    return score / first_fit

if opts.maxiter < 0: opts.maxiter = n.Inf
a.optimize.fmin(fit_func, prm_list,
    maxfun=n.Inf, maxiter=opts.maxiter, ftol=1e-100, xtol=1e-100)

'''
def anneal(func, x0, T_x, cooling=lambda i: 3*(n.cos(i/50.)+1), 
        maxiter=1000, verbose=True):
    score = func(x0)
    for i in range(maxiter):
        if verbose: print 'Step:', i, 'T fraction:', cooling(i)
        delta = n.random.normal(scale=1., size=x0.shape) * cooling(i) * T_x
        n_x = x0 + delta
        n_score = func(n_x)
        if n_score < score: x0, score = n_x, n_score
        if verbose:
            print 'Best score:', score
            print 'Best x:', x0
    return x0, score

prm_list = n.array(prm_list)
x0, score = anneal(fit_func, prm_list, .0001*prm_list)
print x0
print 'Score:', score
print 'Done!'
'''
