#! /usr/bin/env python
"""
Apply calibration parameters to a data set.
"""

import aipy as a
import numpy as np
import optparse,sys,os

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
freqs = np.linspace(uv['sfreq'],uv['sfreq']+uv['nchan']*uv['sdf'],uv['nchan'])
aa = a.cal.get_aa(opts.cal,freqs)
del(uv)

def mfunc(uv,p,d):
    uvw,t,(i,j) = p
    pol = a.miriad.pol2str[uv['pol']]
    aa.set_active_pol(pol)
    G_ij = 1./aa.passband(i,j)
    if i != j:
        E_ij = np.exp(-1j*2*np.pi*aa.get_phs_offset(i,j))
    else:
        E_ij = np.ones_like(d)
    d_fix = d * G_ij * E_ij
    return p,d_fix

for infile in args:
    outfile = infile+'C'
    print infile,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping....'
        continue 

    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc,append2hist='APPLY_CAL: Applied calibration information from %s' % opts.cal)
