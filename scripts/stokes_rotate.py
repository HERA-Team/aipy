#!/usr/bin/env python
"""
Rotate linearly polarized data into Stokes' I,Q,U,V
"""
import aipy as a
import numpy as np
import optparse,sys,os

o = optparse.OptionParser()
o.set_usage('stokes_rotate.py *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])


for uvfile in args:
    
    infile = uvfile
    outfile = infile+'S'
    
    print infile,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping'
        continue

    uv = a.pol.UV(uvfile)
    DD = {}
    for (uvw,t,bl),d,f in uv.all(raw=True):
        plzn = uv.read_pol()
        if not bl in DD.keys(): DD[bl] = {}
        if not t in DD[bl].keys(): DD[bl][t] = {}
        if not plzn in DD[bl][t].keys():
            DD[bl][t][plzn] = np.ma.array(d,mask=f)
    del(uv)
    
    for bl in DD:
        for t in DD[bl]:
            DD[bl][t] = a.pol.xy2stokes(DD[bl][t])
    
    def mfunc(uv,p,d,f):
        uvw,t,bl = p
        print uvi['pol']
        plzn = uvi.read_pol()
        print bl,plzn 
        if plzn == 'xx':
            uvo.write_pol('I')
            print '-->',uvo.read_pol()
            return p,DD[bl][t]['I'],f
        if plzn == 'xy': 
            uvo.write_pol('Q')
            print '-->',uvo.read_pol()
            return p,DD[bl][t]['Q'],f
        if plzn == 'yx': 
            uvo.write_pol('U')
            print '-->',uvo.read_pol()
            return p,DD[bl][t]['U'],f
        if plzn == 'yy': 
            uvo.write_pol('V')
            print '-->',uvo.read_pol()
            return p,DD[bl][t]['V'],f
    
    uvi = a.pol.UV(infile)
    uvo = a.pol.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,raw=True,mfunc=mfunc,append2hist='XY --> STOKES \n')
