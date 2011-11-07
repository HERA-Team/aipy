import timeit

setup = '''
import numpy as n, aipy as a
freqs = n.arange(.1,.2,.0001)
beam = a.fit.Beam(freqs)
ants = [
    a.fit.Antenna(  0,   0,  0, beam),
    a.fit.Antenna(  0, 100,  0, beam),
    a.fit.Antenna(100,   0,  0, beam),
    a.fit.Antenna(100, 100,  0, beam),
]
aa = a.fit.AntennaArray(('45:00','90:00'), ants)
aa.set_jultime(2455400.1)
s_eqs = n.array([[0,1,0]]*100).transpose()
'''
expr = '''
aa.gen_phs(s_eqs,0,1,)
'''

t = timeit.Timer(expr, setup=setup)
print (t.timeit(number=10)/ 10) * 1e3, 'ms'
