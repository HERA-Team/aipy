#! /usr/bin/env python
"""
Models visibilities predicted by a Healpix map, creating a new Miriad UV
file containing either the simulated data, or the residual when the model
is removed from measured data.

Author: Aaron Parsons
Date: 03/06/08
Revisions:
"""
import aipy as a, numpy as n, ephem as e, pylab as p, sys, optparse

o = optparse.OptionParser()
o.set_usage('mdlpix.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-m', '--map', dest='map',
    help='The Healpix map to use for simulation.')
o.add_option('-l', '--loc', dest='loc', 
    help='The array location to use for simulation (see loc.py).')
p.add_option('-s', '--sim', dest='sim', action='store_true',
    help='Output a simulated dataset (rather than subtracting).')
o.add_option('--iepoch', dest='iepoch', default=2000., type='float',
    help='The epoch of coordinates in the map.')
o.add_option('-x', '--decimate', dest='decimate', type='int', default=1,
    help='Only use every Nth sample.')
p.add_option('-n', '--noiselev', dest='noiselev', default=0., type='float',
    help='RMS amplitude of noise added to each UV sample of simulation.')
opts,args = o.parse_args(sys.argv[1:])


# Initialize an AntennaArray for simulation
uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
(uvw,t,(i,j)),d = uv.read()
aa.set_jultime(t)
del(uv)

# Open the HealpixMap containing the sky used for simulation
h = a.healpix.HealpixMap()
h.from_fits(opts.map)
fluxes = h.map
s_eqs = h.px2crd(n.arange(h.Npix()), crd_type='vec').tranpose()
m_precess = a.coord.convert_m('eq','eq', iepoch=opts.iepoch, oepoch=aa.epoch)
s_eqs = n.dot(m_precess, s_eqs)
#indices = ???
#mfreqs = ???

# A pipe for just outputting the model
curtime = None
def mdl(uv, p, d):
    global curtime
    uvw, t, (i,j) = p
    if i == j: return p, d
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
    d = aa.sim(i, j, s_eqs, fluxes, #indices=indices, mfreqs=mfreqs,
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
    sd = aa.sim(i, j, s_eqs, fluxes, #indices=indices, mfreqs=mfreqs,
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

