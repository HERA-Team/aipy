#! /usr/bin/env python
"""
Models visibilities for various catalog sources and creates a new Miriad UV
file containing either the simulated data, or the residual when the model
is removed from measured data.

Author: Aaron Parsons
"""
import numpy as n, aipy as a, optparse, os, sys, ephem

o = optparse.OptionParser()
o.set_usage('mdlvis.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, loc=True, src=True)
o.add_option('--sim', dest='sim', action='store_true',
    help='Output a simulated dataset (rather than subtracting).')
o.add_option('-f', '--flag', dest='flag', action='store_true',
    help='If outputting a simulated data set, mimic the data flagging of the original dataset.')
o.add_option('-m', '--map', dest='map',
    help='The Healpix map to use for simulation input.')
o.add_option('--iepoch', dest='iepoch', default=ephem.J2000, 
    help='The epoch of coordinates in the map. Default J2000.')
o.add_option('--freq', dest='freq', default=.150, type='float',
    help='Frequency of flux data in map.')
o.add_option('-n', '--noiselev', dest='noiselev', default=0., type='float',
    help='RMS amplitude of noise added to each UV sample of simulation.')
o.add_option('--nchan', dest='nchan', default=256, type='float',
    help='Number of channels in simulated data if no input data to mimic.')
o.add_option('--sfreq', dest='sfreq', default=.75, type='float',
    help='Start frequency (GHz) in simulated data if no input data to mimic.')
o.add_option('--sdf', dest='sdf', default=.150/256, type='float',
    help='Channel spacing (GHz) in simulated data if no input data to mimic.')
o.add_option('--inttime', dest='inttime', default=10, type='float',
    help='Integration time (s) in simulated data if no input data to mimic.')
o.add_option('--startjd', dest='startjd', default=2454600., type='float',
    help='Julian Date to start observation if no input data to mimic.')
o.add_option('--endjd', dest='endjd', default=2454601., type='float',
    help='Julian Date to end observation if no input data to mimic.')
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
p,d,f = uv.read(raw=True)
no_flags = n.zeros_like(f)
del(uv)

# Generate a model of the sky with point sources and a pixel map

mfq,asz = [], []
# Initialize point sources
if not opts.src is None:
    cat = a.scripting.parse_srcs(opts.src, force_cat=True)
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
        if not opts.src is None:
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
    if opts.sim:
        d = sd
        if not opts.flag: f = no_flags
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

