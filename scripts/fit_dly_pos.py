#! /usr/bin/env python
"""
A script for fitting the delays, positions, and offsets of antennas given
starting parameters in "aipy.loc" and a list of sources (from "aipy.src").

Author: Aaron Parsons
Date: 6/03/07
Revisions:
    12/11/07 arp    Ported to use new miriad file interface
"""

import aipy, numpy, sys, os
from optparse import OptionParser

def calc_gain(phsdata):
    img = numpy.fft.ifft(phsdata.filled(0))
    return img * fit_window(img.shape[0])

def fit_window(width):
    w = numpy.arange(width)
    w = numpy.where(w > width/2, width - w, w)
    return 10.**(-w)

p = OptionParser()
p.set_usage('fit_dly_pos.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-w', '--swath', dest='swath', default=0, type='int',
    help='Number of bins around center to wipe out when removing a src.')
p.add_option('-c', '--cat', dest='cat',
    help='A list of several sources (separated by commas) to use.')
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
p.add_option('-x', '--decimate', dest='decimate', default=1, type='int',
    help='Only use every Nth time step in fitting calculations.')
p.add_option('-p', '--prms', dest='prms', 
    help='Comma delimited list of paramters to fit (x,y,z,delay,offset).')
p.add_option('-a', '--ants', dest='ants', default='1,2,3',
    help='Comma delimited list of antennas to fit. Default 1,2,3')
opts, args = p.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
aa = aipy.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

srcs = opts.cat.split(',')
cat = aipy.src.get_catalog(srcs, type='sim')

prms = opts.prms.split(',')
ants = eval('['+opts.ants+']')
prm_dict = {}
for a in ants: prm_dict[a] = prms
print prm_dict
start_prms = aa.get_params(prm_dict)

prm_list, key_list = aipy.fit.flatten_prms(start_prms)
first_fit = None    # Used to normalize fit values to the starting fit

def fit_func(prms):
    global first_fit
    prms = aipy.fit.reconstruct_prms(prms, key_list)
    aipy.fit.print_params(prms)
    aa.set_params(prms)
    score = 0
    cnt = 0
    curtime = -1
    src_phs = {}
    #data = []
    for uvfile in args:
        print uvfile
        uvi = aipy.miriad.UV(uvfile)
        uvi.select('auto', 0, 0, include=False)
        for p,d in uvi.all():
            uvw, t, (i,j) = p
            # Use only every Nth integration, if decimation is specified.
            if curtime != t:
                curtime = t
                cnt = (cnt + 1) % opts.decimate
                for s in src_phs:
                    score += abs(src_phs[s]).sum()
            if cnt != 0: continue
            aa.set_jultime(t)
            cat.compute(aa)
            for s in cat:
                ds = d.filled(0)
                # Remove all other interfering sources from data
                for rs in cat:
                    if rs != s: ds = aa.rmsrc(ds, cat[rs], i, j)
                ds = numpy.ma.array(ds, mask=d.mask)
                try:
                    ds = aa.phs2src(ds, cat[s], i, j)
                    src_phs[s] = src_phs.get(s,0) + calc_gain(ds)
                except(aipy.ant.PointingError): pass
    score = 1 / score
    if first_fit is None:
        print 'Base score:', score
        first_fit = score
    print score / first_fit
    #import pylab
    #pylab.clf()
    #data = numpy.array(data)
    #data = numpy.fft.ifft(data, axis=1)
    #pylab.imshow(numpy.abs(data), aspect='auto')
    #pylab.show()
    return score / first_fit

aipy.optimize.fmin(fit_func, prm_list)
print 'Done!'

