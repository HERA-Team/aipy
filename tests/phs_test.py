#! /usr/bin/env python
import aipy as a, numpy as n

freqs = n.arange(1.,2.,.1)
loc = ('0:00', '0:00')
beam = a.ant.Beam(freqs)
ants = [a.ant.Antenna(0,0,0,beam), a.ant.Antenna(100,0,0,beam)]
aa = a.ant.AntennaArray(ants=ants, location=loc)
aa.set_jultime(2454490)
src = a.ant.RadioFixedBody('0','0')
src.compute(aa)
print 'RA=%s' % src.ra, 'DEC=%s' % src.dec
phs = aa.gen_phs(src, 0, 1)
seq1 = a.coord.radec2eq((src.ra, src.dec))
print 'VEC=%s' % seq1
assert(n.all(aa.gen_phs(seq1, 0, 1) == phs))
seq2 = n.array([seq1] * 4)
for p in aa.gen_phs(seq2, 0, 1): assert(n.all(p == phs))
print 'Phases correct'
print 'Passed.'
