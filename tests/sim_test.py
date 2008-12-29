#! /usr/bin/env python
import aipy as a, numpy as n, time, os, sys

uvi = a.miriad.UV('test.uv')
freqs = n.arange(uvi['nchan'], dtype=n.float) * uvi['sdf'] + uvi['sfreq']
loc = ('0:00', '0:00')
beam = a.sim.BeamFlat(freqs)
ants = [
    a.sim.Antenna(0,0,0,beam),
    a.sim.Antenna(100,0,0,beam),
    a.sim.Antenna(0,100,0,beam),
    a.sim.Antenna(100,100,0,beam),
]
aa = a.sim.AntennaArray(ants=ants, location=loc)
src = a.sim.RadioFixedBody('0','0', 1e4)
if os.path.exists('test.uvs'):
    print 'test.uvs exists... skipping'
    sys.exit(0)
uvo = a.miriad.UV('test.uvs', status='new')
uvo.init_from_uv(uvi)

curtime = None
cat = a.sim.SrcCatalog([src] * 4)
def mfunc(uv, p, data):
    global curtime
    uvw,t,(i,j) = p
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        src.compute(aa)
        s_eqs = cat.get_crds('eq')
        fluxes = cat.get_fluxes()
        indices = cat.get_indices()
        mfreqs = cat.get_mfreqs()
        aa.sim_cache(s_eqs=s_eqs,fluxes=fluxes,indices=indices,mfreqs=mfreqs)
    d = aa.sim(i, j)
    return p, n.ma.array(d, mask=data.mask)

t0 = time.time()
uvo.pipe(uvi, mfunc=mfunc)
print 'Simulation finished in %f seconds.' % (time.time() - t0)
del(uvo)
print 'Checking against knowngood.uvs ...'
uv1 = a.miriad.UV('knowngood.uvs')
uv2 = a.miriad.UV('test.uvs')

for (p1,d1),(p2,d2) in zip(uv1.all(), uv2.all()):
    dif = (d1.filled(0) - d2.filled(0)).sum()
    if dif > 1e-10:
        raise ValueError('Failed verification (%f > 1e-10).' % dif)
print 'Passed.'
