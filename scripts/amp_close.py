#!/usr/bin/env python

"""
Get a first-pass gain calibration by enforcing amplitude closure.
"""

import aipy as a
import numpy as np
from matplotlib import pylab as pl
import optparse,sys,os
from time import time

o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('amp_close.py [options] *.uv')
o.add_option('-C','--cal',dest='cal',default=None,help='Cal file to use')
o.add_option('-f','--flag',dest='flag',type='str',help='Antennae to flag. Ex. 1,2,3')
o.add_option('-p','--plots',dest='plots',default=None,action='store_true',help='Plots? [default=no]')
o.add_option('-s','--src',dest='src',default='None',help='Flux Calibration Source [default flux = 1]')
o.add_option('--minuv',dest='minuv',type='float',default=30,help='Minimum baseline length to consider [default=30lambda]')
o.add_option('--wgt',dest='wgt',type='str',default='equal',help='Weighting scheme, options are equal, radial, var, snr [default=equal]')
opts,args = o.parse_args(sys.argv[1:])

t0 = time()

#Parse Arguments:
FlagAnt = []
if opts.flag:
    for antstr in opts.flag.split(','):
        FlagAnt.append(int(antstr))
def uvlen(A):
    return np.sqrt(np.dot(A,A))
aa = a.cal.get_aa(opts.cal,np.array([0.15]))

#Find source Flux:
if opts.src == None: flx = 1.
else: 
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src,'helm,misc')
    cat = a.cal.get_catalog(opts.cal,srclist,cutoff,catalogs)
    cat.compute(aa)
    flx_list = cat.get_jys()
    flx = 0.
    print "Setting Jansky levels to the following sources:"
    for i,src in enumerate(srclist):
        print ' '*5,src,flx_list[i].squeeze(),'Jys'
        flx += flx_list[i].squeeze()

#Read in Data:
DD,cnt,bl_len,AntNos = {},{},{},[]

for uvfile in args:
    print 'Reading',uvfile
    uv = a.miriad.UV(uvfile)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        pol = a.miriad.pol2str[uv['pol']]
        #Selection Criteria:
        if i == j: continue
        if np.sum(f) == uv['nchan']: continue
        if pol[0] != pol[1]: continue
        if i in FlagAnt or j in FlagAnt: continue
        if not (i,j) in bl_len.keys(): bl_len[(i,j)] = uvlen(aa.get_baseline(i,j))
        if bl_len[(i,j)]*0.15 <= opts.minuv: continue
        #Populate arrays:
        if not pol in DD.keys(): 
            DD[pol] = {}
            cnt[pol] = {}
        if not (i,j) in DD[pol].keys(): 
            DD[pol][(i,j)] = np.zeros(uv['nchan'],dtype=np.complex)
            cnt[pol][(i,j)] = np.zeros(uv['nchan'])
        DD[pol][(i,j)] += np.abs(d)
        cnt[pol][(i,j)] += np.where(f,0,1)
        if not i in AntNos: AntNos.append(i)
        if not j in AntNos: AntNos.append(j)
pols = DD.keys()
AntNos.sort()
pols.sort()
dnu = uv['sdf']
del(uv)

Npol,Nant = len(DD.keys()),len(AntNos)
Nbl = int(Nant*(Nant-1)/2)

#Integrate and all that shit...
G_ij,W_ij = {},{}
for pol in pols:
    G_ij[pol] = {}
    W_ij[pol] = {}
    for bl in DD[pol]:
        DD[pol][bl] = np.where(cnt[pol][bl] == 0, 0, DD[pol][bl]/cnt[pol][bl])
        G_ij[pol][bl] = np.log(np.mean(DD[pol][bl]))
        if opts.wgt == 'equal' : W_ij[pol][bl] = 1.
        if opts.wgt == 'radial': W_ij[pol][bl] = bl_len[bl]
        if opts.wgt == 'var'   : W_ij[pol][bl] = 1./np.log(np.var(DD[pol][bl]))
        if opts.wgt == 'snr'   : W_ij[pol][bl] = np.log(np.mean(DD[pol][bl]))/np.var(DD[pol][bl])

#Do the linear algebra
G = {}
for pol in pols:
    A = np.zeros((Nbl,Nant))
    W = np.zeros((Nbl,Nbl))
    b = np.zeros(Nbl)
    
    #populatin the countryside...
    index = 0
    for i,ant1 in enumerate(AntNos):
        for j,ant2 in enumerate(AntNos):
            if not (ant1,ant2) in G_ij[pol].keys(): continue
            b[index] = G_ij[pol][(ant1,ant2)]
            A[index,i],A[index,j] = 1.,1.
            W[index,index] = W_ij[pol][(ant1,ant2)]
            index += 1
    W /= np.sum(W)

    G[pol] = 10**np.dot(np.linalg.inv(np.dot(np.dot(A.T,W),A)),np.dot(np.dot(A.T,W),b))
    G[pol] /= flx*np.mean(G[pol])

print "'amps' : {"
for i,ant in enumerate(AntNos):
    ampstr = str(ant)+' : {'
    for pol in pols:
        ampstr += "'"+pol[0]+"' : "+str(G[pol][i])
        if pol != pols[-1]: ampstr += ', '
    ampstr += ' },'
    print ampstr
print '},'

t1 = time()
print 'Computation time =',t1-t0,'s'

if opts.plots:
    print 'Generating plots... this may take a while...'
    AntNos = np.array(AntNos)
    figcnt = 0

    for pol in pols:
        pl.figure(figcnt)
        DDD = np.zeros((Nant,Nant))
        for i,ant1 in enumerate(AntNos):
            for j,ant2 in enumerate(AntNos):
                if not (ant1,ant2) in G_ij[pol].keys(): continue
                DDD[i,j] = 10**G_ij[pol][(ant1,ant2)]/flx
        DDD = np.ma.array(DDD,mask=np.where(DDD==0,1,0))
        pl.imshow(DDD,aspect='auto',interpolation='nearest')
        pl.colorbar()
        pl.title('Gains by baseline, polarization %s'%pol[0])
        pl.draw()
        figcnt += 1
    
    if opts.wgt != 'equal':
        for pol in pols:
            pl.figure(figcnt)
            WWW = np.zeros((Nant,Nant))
            for i,ant1 in enumerate(AntNos):
                for j,ant2 in enumerate(AntNos):
                    if not (ant1,ant2) in W_ij[pol].keys(): continue
                    WWW[i,j] = W_ij[pol][(ant1,ant2)]
            WWW = np.ma.array(WWW/np.sum(WWW),mask=np.where(WWW==0,1,0))
            pl.imshow(WWW,aspect='auto',interpolation='nearest')
            pl.colorbar()
            pl.title('Weights for each baseline, polarization %s'%pol[0])
            pl.draw()
            figcnt += 1 

    pl.figure(figcnt)
    for pol in pols:
        pl.plot(AntNos,G[pol],'o',label=pol[0])
    pl.xlabel('Antenna Number')
    pl.ylabel('Gain (units)')
    pl.legend()
    pl.draw()
    figcnt += 1

    pl.show()

