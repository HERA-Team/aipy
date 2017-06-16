#!/usr/bin/env python

"""
Do a one iteration delay-solving algorithm.
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
o.add_option('-C','--cal',dest='cal',type='str',help='Cal file to use.')
o.add_option('-f','--flag',dest='flag',type='str',help='Antennae to flag. Ex. 1,2,3')
o.add_option('-p','--plots',dest='plots',default=None,action='store_true',help='Plots?')
opts,args = o.parse_args(sys.argv[1:])

t0 = time()
if opts.cal != None:
    exec('from %s import prms'%opts.cal)
    old_delays = prms['delays'].copy()
    for ant in old_delays:
        if type(old_delays[ant]) is not dict:
            old_delays[ant] = {'x': prms['delays'][ant]}

#Parse Arguments
FlagAnt = []
if opts.flag:
    for antstr in opts.flag.split(','):
        FlagAnt.append(int(antstr))

#Define a function to give me a cleaned delay.
def clean_delay(d):
    d = d.squeeze()
    d = np.fft.ifft(d)
    samples = np.logical_not(f).astype(float)
    ker = np.fft.ifft(samples)
    gain = np.sqrt(np.mean(samples**2))
    d,info = a.deconv.clean(d,ker,tol=1e-6)
    d += info['res']/gain
    d = np.fft.fftshift(d)
    return d

def length(A):
    return np.sqrt(np.dot(A,A))

aa = a.cal.get_aa(opts.cal,np.array([0.15]))

#Read in Data
DD = {}
bl_len = {}
AntNos = []
for uvfile in args:
    print 'Reading',uvfile
    uv = a.miriad.UV(uvfile)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if i == j: continue
        if not (i,j) in bl_len.keys(): bl_len[(i,j)] = length(aa.get_baseline(i,j))
        if bl_len[(i,j)]*0.15 <= 30.: continue
        pol = a.miriad.pol2str[uv['pol']]
        if pol == 'xy' : continue
        if pol == 'yx' : continue
        if i in FlagAnt: continue
        if j in FlagAnt: continue
        if not pol in DD.keys(): 
            DD[pol] = {}
        if not (i,j) in DD[pol].keys(): DD[pol][(i,j)] = np.zeros(uv['nchan'],dtype=np.complex)
        #Delay transform
        d = np.ma.array(d)
        DD[pol][(i,j)] += d
        if not i in AntNos: AntNos.append(i)
        if not j in AntNos: AntNos.append(j)
pols = DD.keys()
AntNos.sort()
pols.sort()
delays = np.linspace(-0.5/uv['sdf'],0.5/uv['sdf'],uv['nchan'])
del(uv)

Npol = len(DD.keys())
Nant = len(AntNos)
Nbl = int(Nant*(Nant-1)/2)

#Find the peak of each delay spectrum
D_ij = {}
W_ij = {}
for pol in pols:
    D_ij[pol] = {}
    W_ij[pol] = {}
    for bl in DD[pol]:
        DD[pol][bl] = np.abs(clean_delay(DD[pol][bl]))
        DD[pol][bl] /= DD[pol][bl].max()
        #multiply by some prior.
        for i,x in enumerate(DD[pol][bl]):
            if x == 1: D_ij[pol][bl] = delays[i]
        #W_ij[pol][bl] = 1./(np.std(DD[pol][bl])**2) 
        W_ij[pol][bl] = 1. 

#Construct matrices to do linear algebra and shit.

Tau = {}
for pol in pols:
    A = np.zeros((Nbl+1,Nant))
    W = np.zeros((Nbl+1,Nbl+1))
    b = np.zeros(Nbl+1)
        
    #populatin the countryside, populatin the villages!
    index = 0
    for i,ant1 in enumerate(AntNos):
        for j,ant2 in enumerate(AntNos):
            if not (ant1,ant2) in D_ij[pol].keys(): continue
            #if j <= i: continue
            b[index] = D_ij[pol][(ant1,ant2)]
            A[index,i],A[index,j] = 1.,-1.    
            W[index,index] = W_ij[pol][(ant1,ant2)]
            index += 1
    #Set Constraints
    W /= np.sum(W)
        
    Tau[pol] = np.dot(np.linalg.inv(np.dot(np.dot(A.T,W),A)),np.dot(np.dot(A.T,W),b))
    Tau[pol] -= np.mean(Tau[pol])

#Print things so you can just copy/paste into a cal file.
if not opts.cal is None:
    for i,ant in enumerate(AntNos):
        calstr = str(ant)+' : { '
        for pol in pols:
            calstr +="'"+pol[0]+"'"+' : '
            Tau[pol][i] += old_delays[ant][pol[0]]
            calstr += str(Tau[pol][i])
            if pol != pols[-1]: calstr += ', '
        calstr += ' }'
        if i != len(AntNos)-1: calstr += ','
        print calstr
else:
    for i,ant in enumerate(AntNos):
        calstr = str(ant)+' : { '
        for pol in pols:
            calstr +="'"+pol[0]+"'"+' : '
            calstr += str(Tau[pol][i])
            if pol != pols[-1]: calstr += ', '
        calstr += ' }'
        if i != len(AntNos)-1: calstr += ','
        print calstr

print 'Computing time =', time() - t0,'s for',Nbl,'baselines.'


if opts.plots:
    print 'Generating plots.... this may take a while...'
    figcnt = 0
    AntNos = np.array(AntNos) 
    
    def issubplotedge(rows,cols,index):
        inds = np.arange(1,cols*rows+1,1)
        inds.shape = (rows,cols)
        left = not np.argwhere(inds==index).squeeze()[1]
        bottom = not (np.argwhere(inds==index).squeeze()[0]+1)%rows
        return left,bottom
    
    for pol in pols:
        pl.figure(figcnt)
        
        nplots = len(DD[pol].keys())
        rows = Nant
        cols = Nant
        bls = DD[pol].keys()
        bls.sort()

        for plno,bl in enumerate(bls):
            #Make the plots
            i,j = np.argwhere(AntNos==bl[0]).squeeze(),np.argwhere(AntNos==bl[1]).squeeze()
            plindex = i*Nant+j
            
            ax = pl.subplot(rows,cols,plindex)
           
            pl.vlines(D_ij[pol][bl]/bl_len[bl],0,1.1,color='k')
            pl.vlines((Tau[pol][i]-Tau[pol][j])/bl_len[bl],0,1.1,color='r')
            pl.vlines((-1,1),0,1.1,linestyles='dotted',color='k')
            pl.plot(delays/bl_len[bl],DD[pol][bl],'b')
           
            pl.xlim([-1.5,1.5])
            pl.ylim([0,1.1])
            pl.yticks([])
            pl.xticks([])

            if i == 0: pl.text(0.5,1.2,str(AntNos[j]),transform=ax.transAxes)
            if j == Nant-1: pl.text(1.2,0.5,str(AntNos[i]),transform=ax.transAxes)
    
        #Give a sample plot
        pl.subplot(337)
        bl = bls[0]
        i,j = np.argwhere(AntNos==bl[0]).squeeze(),np.argwhere(AntNos==bl[1]).squeeze()
        pl.vlines(D_ij[pol][bl]/bl_len[bl],0,1.1,color='k')
        pl.vlines((Tau[pol][i]-Tau[pol][j])/bl_len[bl],0,1.1,color='r')
        pl.vlines((-1,1),0,1.1,linestyles='dotted',color='k')
        pl.plot(delays/bl_len[bl],DD[pol][bl],'b')
        pl.title('Sample baseline, %d_%d'%bl)
        pl.xlim([-1.5,1.5])
        pl.ylim([0,1.1])
        pl.xlabel('Delay, baseline lenghth.')
        pl.ylabel('Amplitude')

        pl.subplots_adjust(wspace=0,hspace=0,top=0.95,bottom=0.05)
        pl.draw()
        figcnt += 1
        
    #Plot a matrix of all delays
    for pol in pols:
        pl.figure(figcnt)
        DDD = np.zeros((Nant,Nant))
        for i,ant1 in enumerate(AntNos):
            for j,ant2 in enumerate(AntNos):
                if not (ant1,ant2) in D_ij[pol].keys(): continue
                DDD[i,j] = D_ij[pol][(ant1,ant2)]
        DDD = np.ma.array(DDD,mask=np.where(DDD==0,1,0))
        pl.matshow(DDD,aspect='auto')
        pl.colorbar()
        pl.draw()
        figcnt += 1
        
    
    #Delay vs. Antenna
    pl.figure(figcnt)
    for pol in pols:
        pl.plot(AntNos,Tau[pol],'o',label=pol[0])
        pl.xlabel('Antenna Number')
        pl.ylabel('Delay (ns)')
    pl.legend()
    pl.draw()

    pl.show()
