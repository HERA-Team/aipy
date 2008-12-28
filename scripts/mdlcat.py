#! /usr/bin/env python
"""
Models visibilities for various catalog sources and creates a new Miriad UV
file containing either the simulated data, or the residual when the model
is removed from measured data.

Author: Aaron Parsons
Date: 03/06/08
Revisions:
"""
import sys, os, numpy, aipy, optparse

p = optparse.OptionParser()
p.set_usage('mdlcat.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-s', '--sim', dest='sim', action='store_true',
    help='Output a simulated dataset (rather than subtracting).')
p.add_option('-c', '--cat', dest='cat',
    help='A list of several sources (separated by commas) to remove.')
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
p.add_option('-b', '--use_bp', dest='use_bp', action='store_true',
    help='Use bandpass information (otherwise assumes files are bp calibrated.')
p.add_option('-n', '--noiselev', dest='noiselev', default=0., type='float',
    help='RMS amplitude of noise added to each UV sample of simulation.')
opts, args = p.parse_args(sys.argv[1:])

srcs = opts.cat.split(',')
cat = aipy.src.get_catalog(srcs, type='sim')

uv = aipy.miriad.UV(args[0])
aa = aipy.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'],
    use_bp=opts.use_bp)
(uvw,t,(i,j)),d = uv.read()
aa.set_jultime(t)
cat.compute(aa)
s_eqs,fluxes,indices,mfreqs = cat.get_vecs()
del(uv)

# A pipe for just outputting the model
curtime = None
def mdl(uv, p, d):
    global curtime
    uvw, t, (i,j) = p
    if i == j: return p, d
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        cat.compute(aa)
    d = aa.sim(i, j, s_eqs, fluxes, indices=indices, mfreqs=mfreqs, 
        pol=aipy.miriad.pol2str[uv['pol']])
    d = numpy.ma.array(d, mask=numpy.zeros_like(d))
    # Add on some noise for a more realistic experience
    noise_amp = numpy.random.random(d.shape) * opts.noiselev
    noise_phs = numpy.random.random(d.shape) * 2*numpy.pi * 1j
    noise = noise_amp * numpy.exp(noise_phs)
    return p, d + noise

# A pipe to use for removing the model
def rm(uv, p, d):
    global curtime
    uvw, t, (i,j) = p
    if i == j: return p, d
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        cat.compute(aa)
    sd = aa.sim(i, j, s_eqs, fluxes, indices=indices, mfreqs=mfreqs, 
        pol=aipy.miriad.pol2str[uv['pol']])
    return p, d - sd

for filename in args:
    print filename
    uvofile = filename + 's'
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    if opts.sim: uvo.pipe(uvi, mfunc=mdl)
    else: uvo.pipe(uvi, mfunc=rm)

