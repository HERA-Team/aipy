#!/usr/bin/env python

"""
Do one iteration of a delay-solving algorithm.
"""

import aipy as a
import numpy as np
from matplotlib import pylab as pl
import optparse,sys,os
from time import time
from glob import glob

o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('find_delay.py [options] *.uv.src')
o.add_option('-C','--cal',dest='cal',default=None,help='Cal file to use')
o.add_option('-f','--flag',dest='flag',type='str',help='Antennae to flag. Ex. 1,2,3')
o.add_option('--fixant',dest='fixant',type='str',default=None,help='Antennae whose delays you want to fix at calfile value')
o.add_option('--refant',dest='refant',default=0,help='Reference Antenna Number (default=0)')
o.add_option('-p','--plots',dest='plots',default=None,action='store_true',help='Plots?')
opts,args = o.parse_args(sys.argv[1:])

t0 = time()

#Parse Arguments
FlagAnt = []
if opts.flag:
    for antstr in opts.flag.split(','):
        FlagAnt.append(int(antstr))

exec('from %s import prms' % opts.cal)
        

#Read in Data
DD = {}
AntNos = []
for uvfile in args:
    print 'Reading',uvfile
    uv = a.miriad.UV(uvfile)
    for (uvw,t,(i,j)),d in uv.all():
        if i == j: continue
        pol = a.miriad.pol2str[uv['pol']]
        if pol == 'xy' : continue
        if pol == 'yx' : continue
        if i in FlagAnt: continue
        if j in FlagAnt: continue
        if not pol in DD.keys(): DD[pol] = {}
        if not (i,j) in DD[pol].keys(): DD[pol][(i,j)] = np.zeros_like(d)
        #Construct flagging kernel
        flags = np.logical_not(d.mask).astype(np.float)
        gain = np.sqrt(np.average(flags**2))
        ker = np.fft.ifft(flags)
        d = d.filled(0)
        #Delay transform
        d = np.fft.ifft(d)
        d = np.ma.array(d)
        d = np.ma.concatenate([d[d.shape[0]/2:],d[:d.shape[0]/2]],axis=0)
        DD[pol][(i,j)] += np.abs(d)
        if not i in AntNos: AntNos.append(i)
        if not j in AntNos: AntNos.append(j)
delays = np.arange(-0.5/uv['sdf'],0.5/uv['sdf'],1/(uv['sdf']*uv['nchan']))
pols = DD.keys()
AntNos.sort()
pols.sort()
del(uv)

#Pick out maximum values DO A BETTER JOB OF THIS!!!
Phi_ij = {}
for pol in pols:
    if not pol in Phi_ij.keys(): Phi_ij[pol] = {}
    for ij in DD[pol].keys():
        for i,x in enumerate(DD[pol][ij]):
            if x == np.max(DD[pol][ij]): Phi_ij[pol][ij] = delays[i]
del(DD)

#Construct matrices to do linear algebra and shit.
Npol = len(Phi_ij.keys())
Nant = len(AntNos)
Nbl = int(Nant*(Nant-1)/2)


#Set delays you want to fix.
if opts.refant == None: refant = AntNos[0]
else: refant = opts.refant

if not opts.fixdelay: FixDelay = {}
if opts.fixdelay == 'all':
    FixDelay = {}
    for pol in pols:
        try: FixDelay[pol] = [prms['delays'][i][pol] for i in AntNos]
        except(KeyError): FixDelay[pol] = [prms['delays'][i] for i in AntNos]
        FixDelayIndex[pi] = range(len(AntNos))
else: 
    FixDelay = {}
    FixDelayIndex = []
    for pol in pols:
        FixDelay[pol] = []
        for i,ant in AntNos:
            if ant in opts.fixdelay.split(','):
                try: FixDelay[pol].append(prms['delays'][ant][pol])
                except(KeyError): FixDelay[pol].append(prms['delays'][ant])
                if not i in FixDelayIndex: FixDelayIndex.append(i)
            
Tau = {}
CC = {}
for pol in pols:
    
    A,b = {},{}
    A = np.zeros((Nbl+len(FixDelay[pol])+1,Nant))
    b = np.zeros(Nbl+len(FixDelay[pol])+1)

    #populatin the countryside, populatin the villages!
    index = 0
    for i,ant1 in enumerate(AntNos):
        if ant1 == refant:
            A[-1,i] = 1.
            b[-1] = 0.
        for j,ant2 in enumerate(AntNos):
            if j <= i: continue
            b[index] = Phi_ij[pol][(ant1,ant2)]
            A[index,i],A[index,j] = 1.,-1.    
            index += 1
    #Set Constraints
    for i,Di in enumerate(FixDelayIndex):
        A[Nbl+i,Di] = 1.
        b[Nbl+i] = FixDelay[pol][i]
    Tau[pol] = np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,b))
    #These are the delays (in ns)!!!

    #Compute closure phases
    CC[pol] = {}
    for i in AntNos:
        for j in AntNos:
            if j <= i: continue
            for k in AntNos:
                if k <= j: continue
                CC[pol][(i,j,k)] = Phi_ij[pol][(i,j)] + Phi_ij[pol][(j,k)] - Phi_ij[pol][(i,k)]

print 'Computing time =', time() - t0,'s for',Nbl,'baselines.'

#Print delays:
for pol in pols:
    for i,d in enumerate(Tau[pol]):
        badstr = ''
        if np.abs(d-np.average(Tau[pol])) >= 2.*np.std(Tau[pol]): badstr = '!!!'
        print 'Antenna',AntNos[i],pol[0],':',d,'ns',badstr


#Generate plots:
if opts.plots:
    #Delay vs. Antenna
    for pol in pols:
        pl.plot(AntNos,Tau[pol],'o',label=pol)
        pl.xlabel('Antenna Number')
        pl.ylabel('Delay (ns)')
    pl.legend()
    pl.show()

    #Histogram of all Closures
    c = {}
    for pol in CC:
        if not pol in c: c[pol] = {}
        for (i,j,k) in CC[pol]:
            if not 'all' in c[pol]: c[pol]['all'] = []
            if not i in c[pol]: c[pol][i] = []
            if not j in c[pol]: c[pol][j] = []
            if not k in c[pol]: c[pol][k] = []
            c[pol]['all'].append(CC[pol][(i,j,k)])
            c[pol][i].append(CC[pol][(i,j,k)])
            c[pol][j].append(CC[pol][(i,j,k)])
            c[pol][k].append(CC[pol][(i,j,k)])

    for i,pol in enumerate(c.keys()):
        for ant in c[pol]:
            pl.hist(c[pol][ant],bins=50)
            pl.title('Antenna %s, %s' % (ant,pol))
            pl.xlabel('Closure phase (ns)')
            pl.ylabel('Number')
            pl.show()
            
