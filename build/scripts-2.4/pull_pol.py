#!/usr/bin/python
"""
Parse large data files by antenna number. Keep antenna 0 in all of them as a phase reference.
"""

import aipy as a
import sys,os,optparse

o = optparse.OptionParser()
o.set_usage('pull_ants.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-p','--pols',dest='pols',type='str',default=8,help='Which polarizations do you like?')
opts,args = o.parse_args(sys.argv[1:])

for filename in args:

    outfile = filename+'P'
    print filename,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping'
        continue

    pols2use = opts.pols.split(',')

    def mfunc(uv,p,d):
        pol = a.miriad.pol2str[uv['pol']] 
        if pol in pols2use: return p,d
        else: return None,None

    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    histstr = 'PULL POLS: pols='
    for pol in pols2use:
        histstr += str(pol)+','
        if pol == pols2use[-1]: histstr += '\n'
    uvo.pipe(uvi,mfunc=mfunc,append2hist=histstr)
    del uvo,uvi
