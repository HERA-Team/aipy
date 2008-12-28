#! /usr/bin/env python
"""
Models visibilities for various catalog sources and creates a new Miriad UV
file containing either the simulated data, or the residual when the model
is removed from measured data.

Author: Aaron Parsons
Date: 03/06/08
Revisions:
    04/30/08    arp Added option of using a Map (stored in a fits file) 
                    for simulation.
    05/07/08    arp Changed to sim_cache() and pipe(raw=True) for 2.5x speedup
"""
import numpy as n, aipy as a, optparse, os, sys, ephem

p = optparse.OptionParser()
p.set_usage('mdlvis.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-s', '--sim', dest='sim', action='store_true',
    help='Output a simulated dataset (rather than subtracting).')
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
p.add_option('-c', '--cat', dest='cat', 
    help='A list of several sources (separated by commas) to remove.')
p.add_option('-m', '--map', dest='map',
    help='The Healpix map to use for simulation.')
p.add_option('--iepoch', dest='iepoch', default=ephem.J2000, 
    help='The epoch of coordinates in the map. Default J2000.')
p.add_option('--freq', dest='freq', default=.150, type='float',
    help='Frequency of flux data in map.')
p.add_option('-n', '--noiselev', dest='noiselev', default=0., type='float',
    help='RMS amplitude of noise added to each UV sample of simulation.')
opts, args = p.parse_args(sys.argv[1:])

# Initialize AntennaArray
uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
p,d,f = uv.read(raw=True)
no_flags = n.zeros_like(f)
del(uv)

# Generate a model of the sky with point sources and a pixel map

mfq = []
asz = []
# Initialize point sources
if not opts.cat is None:
    srcs = opts.cat.split(',')
    cat = a.src.get_catalog(srcs)
    cat.set_params(a.loc.get_src_prms(opts.loc))
    mfq.append(cat.get_mfreqs())
    asz.append(cat.get_angsizes())

# Initialize pixel map
if not opts.map is None:
    h = a.map.Map(fromfits=opts.map)
    px = n.arange(h.npix())
    try: mflx, i_poly = h[px]
    except(ValueError): mflx = h[px]
    px = n.compress(mflx > 0, px)
    try: mflx, i_poly = h[px]
    except(ValueError):
        mflx = h[px]
        i_poly = [n.zeros_like(mflx)]
    mind = i_poly[0]    # Only implementing first index term for now
    mmfq = opts.freq * n.ones_like(mind)
    mfq.append(mmfq)
    x,y,z = h.px2crd(px, ncrd=3)
    m_eq = n.array((x,y,z))
    # Should pixels include information for resolving them?
    asz.append(n.zeros_like(mmfq))
mfq = n.concatenate(mfq)
asz = n.concatenate(asz)
if n.all(asz == 0): asz = None

# A pipe for just outputting the model
curtime = None
def mdl(uv, p, d, f):
    global curtime, eqs
    uvw, t, (i,j) = p
    if i == j: return p, d, f
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        eqs,flx,ind = [],[],[]
        if not opts.cat is None:
            cat.compute(aa)
            eqs.append(cat.get_crds('eq', ncrd=3))
            flx.append(cat.get_fluxes())
            ind.append(cat.get_indices())
        if not opts.map is None:
            m_precess = a.coord.convert_m('eq','eq',
                iepoch=opts.iepoch, oepoch=aa.epoch)
            eqs.append(n.dot(m_precess, m_eq))
            flx.append(mflx); ind.append(mind)
        eqs = n.concatenate(eqs, axis=-1)
        flx = n.concatenate(flx)
        ind = n.concatenate(ind)
        aa.sim_cache(eqs, flx, indices=ind, mfreqs=mfq, angsizes=asz)
    sd = aa.sim(i, j, pol=a.miriad.pol2str[uv['pol']])
    #sd = n.ma.array(sd, mask=n.zeros_like(sd))
    if opts.sim:
        d = sd
        f = no_flags
    else: d -= sd
    if opts.noiselev != 0:
        # Add on some noise for a more realistic experience
        noise_amp = n.random.random(d.shape) * opts.noiselev
        noise_phs = n.random.random(d.shape) * 2*n.pi * 1j
        noise = noise_amp * n.exp(noise_phs)
        d += noise
    return p, d, f

# Run mdl on all files
for filename in args:
    print filename
    uvofile = filename + 's'
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mdl, raw=True)

