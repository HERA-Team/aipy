#!/usr/global/paper/bin/python
import aipy as a
import numpy as np
import optparse,sys,os

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])


for uvfile in args:
    
    infile = uvfile
    outfile = infile+'S'
    
    print infile,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping'
        continue

    uv = a.miriad.UV(uvfile)
    DD = {}
    for (uvw,t,bl),d,f in uv.all(raw=True):
        plzn = a.miriad.pol2str[uv['pol']]
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
        plzn = a.miriad.pol2str[uv['pol']]
        if plzn == 'xx':
            uvo['pol'] = a.miriad.str2pol['I']
            return p,DD[bl][t]['I'],f
        if plzn == 'xy': 
            uvo['pol'] = a.miriad.str2pol['Q']
            return p,DD[bl][t]['Q'],f
        if plzn == 'yx': 
            uvo['pol'] = a.miriad.str2pol['U']
            return p,DD[bl][t]['U'],f
        if plzn == 'yy': 
            uvo['pol'] = a.miriad.str2pol['V']
            return p,DD[bl][t]['V'],f
    
    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,raw=True,mfunc=mfunc,append2hist='XY --> STOKES \n')
