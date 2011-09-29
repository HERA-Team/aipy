#!/usr/bin/python
"""
Parse large data files by antenna number. Keep antenna 0 in all of them as a phase reference.
"""

import aipy as a
import sys,os,optparse

o = optparse.OptionParser()
o.set_usage('pull_ants.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-a','--ants',dest='ants',type='str',default=None,help='List of Antennae to pull.')
opts,args = o.parse_args(sys.argv[1:])

if opts.ants == None:
    print 'You forgot to tell me which antennae to pull. Try using a -a ant1,ant2,ant3 option or something, dude'
    sys.exit()

for filename in args:
    outfile = filename+'A'
    print filename,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping'
        continue
    
    ants2use = opts.ants.split(',') 
    
    def mfunc(uv,p,d):
        uvw,t,(i,j) = p
        if i in ants2use and j in ants2use: return p,d
        else: return None,None
    
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    histstr = 'PULL ANTS: ants='%opts.refant
    for ant in ants2use:
        histstr += str(ant)+','
        if ant == ants2use[-1]: histstr += '\n'
    uvo.pipe(uvi,mfunc=mfunc,append2hist=histstr)
    del uvo,uvi
