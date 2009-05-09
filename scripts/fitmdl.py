#! /usr/bin/env python
"""
A script for fitting parameters of a measuremet equation given 
starting parameters in a cal file and a list of sources.
"""

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('fitmdl.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True,
    cal=True, src=True, dec=True)
o.add_option('--snap', dest='snap', action='store_true',
    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('--fitants', dest='fitants', default='all',
    help='Comma delimited list of antennas to fit. Default "all".')
o.add_option('--aprms', dest='aprms',
    help='Comma delimited list of paramters to fit independently for each antenna.')
o.add_option('--shprms', dest='shprms',
    help='Comma delimited list of parameters to fit whose values will be shared for all antennas.')
o.add_option('--sprms', dest='sprms', 
    help='Source=param pairs for fitting source parameters.')
o.add_option('--aaprms', dest='aaprms',
    help='Comma delimited list of AntennaArray (not Antenna) parameters to fit.  Unless added specifically in cal file, default AntennaArray has no parameters to fit.')
o.add_option('-q', '--quiet', dest='quiet', action='store_true',
    help='Be less verbose.')
o.add_option('--maxiter', dest='maxiter', type='float', default=-1,
    help='Maximum # of iterations to run.  Default is infinite.')
o.add_option('--xtol', dest='xtol', type='float', default=1e-10,
    help='Fractional change sought in it parameters before convergence.  Default 1e-10.')
o.add_option('--ftol', dest='ftol', type='float', default=1e-10,
    help='Fractional tolerance sought in score before convergence.  Default 1e-10.')
o.add_option('--sim_autos', dest='sim_autos', action='store_true',
    help='Use auto-correlations in fitting.  Default is to use only cross-correlations.')

opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
opts.ant += ',cross'
a.scripting.uv_selector(uv, opts.ant, opts.pol)
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
srclist,cutoff = a.scripting.parse_srcs(opts.src)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff)
(uvw,t,(i,j)),d = uv.read()
aa.set_jultime(t)
cat.compute(aa)
if opts.fitants.startswith('all'): fitants = range(uv['nants'])
else: fitants = map(int, opts.fitants.split(','))
del(uv)
if opts.maxiter < 0: opts.maxiter = n.Inf

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
#print prm_dict
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
mfq = cat.get('mfreq')
dbuf = None

# The function to be optimized
def fit_func(prms, filelist, decimate, decphs):
    global first_fit, dbuf
    if first_fit == 0: return 0
    prms = a.fit.reconstruct_prms(prms, key_list)
    if not opts.quiet: a.fit.print_params(prms)
    for ant in fitants:
        for prm in shprms:
            prms[ant][prm] = prms[fitants[0]][prm]
    aa.set_params(prms)
    cat.set_params(prms)
    a1,a2,th = cat.get('srcshape')
    score,nsamples = 0.,0.
    # Cache data from file to avoid hitting disk excessively
    if dbuf is None:
        dbuf = {}
        for uvfile in filelist:
            sys.stdout.write('.') ; sys.stdout.flush()
            uv = a.miriad.UV(uvfile)
            a.scripting.uv_selector(uv, opts.ant, opts.pol)
            uv.select('decimate', decimate, decphs)
            for (uvw,t,(i,j)),d,f in uv.all(raw=True):
                if not dbuf.has_key(t): dbuf[t] = {}
                if not opts.sim_autos and i == j: continue
                bl = a.miriad.ij2bl(i,j)
                dbuf[t][bl] = (d.take(chans), f.take(chans), 
                        n.where(f, 0, n.abs(d)**2).sum(),
                        a.miriad.pol2str[uv['pol']])
    # Process data from cache
    for t in dbuf:
        aa.set_jultime(t)
        cat.compute(aa)
        eqs = cat.get_crds('eq', ncrd=3)
        flx = cat.get_jys()
        dra,ddec = cat.get('ionref')
        aa.sim_cache(eqs, flx, mfreqs=mfq, 
            ionrefs=(dra,ddec), srcshapes=(a1,a2,th))
        for bl in dbuf[t]:
            i,j = a.miriad.bl2ij(bl)
            d,f,nsamp,pol = dbuf[t][bl]
            sim_d = aa.sim(i, j, pol=pol)
            difsq = n.abs(d - sim_d)**2
            difsq = n.where(f, 0, difsq)
            score += difsq.sum()
            nsamples += nsamp
    if nsamples == 0:
        first_fit = 0.
        return 0.
    score = n.sqrt(score / nsamples)
    if first_fit is None: first_fit = score
    if not opts.quiet:
        print
        print 'Score:', score, 
        print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
        print '------------------------------------------------------------'
    return score / first_fit

# Call the optimizer
if not opts.snap:
    rv = a.optimize.fmin(
        fit_func, prm_list,
        args=(args, opts.decimate, opts.decphs),
        full_output=1, disp=0,
        maxfun=opts.maxiter, maxiter=n.Inf, 
        ftol=opts.ftol, xtol=opts.xtol
    )
    prms,score = rv[:2]
    prms = a.fit.reconstruct_prms(prms, key_list)
    print
    a.fit.print_params(prms)
    print 'Score:', score * first_fit, 
    print '(%2.2f%% of %f)' % (100 * score, first_fit)
    print '------------------------------------------------------------'
else:
    for uvfile in args:
        # Figure out what times are in the file
        uv = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        times = [0]
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            if times[-1] != t: times.append(t)
        times = times[1:]
        del(uv)
        decimate = len(times) * opts.decimate
        # Fit each time separately
        for cnt, t in enumerate(times):
            print 'Time:', t
            print 'Iter: %d / %d' % (cnt+1, len(times))
            first_fit,dbuf = None,None
            rv = a.optimize.fmin(
                fit_func, prm_list,
                args=([uvfile], decimate, opts.decimate*cnt + opts.decphs),
                full_output=1, disp=0,
                maxfun=opts.maxiter, maxiter=n.Inf, 
                ftol=opts.ftol, xtol=opts.xtol
            )
            prms,score = rv[:2]
            prms = a.fit.reconstruct_prms(prms, key_list)
            print
            a.fit.print_params(prms)
            print 'Score:', score * first_fit, 
            print '(%2.2f%% of %f)' % (100 * score, first_fit)
            print '------------------------------------------------------------'
