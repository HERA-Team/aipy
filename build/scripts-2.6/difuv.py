#!/usr/global/paper/bin/python
"""
Subtract the data in the second UV file from the data in the first UV file.
I imagine they better have corresponding integrations/data order.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, optparse

o = optparse.OptionParser()
o.set_usage('difuv.py [options] file1.uv file2.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

assert(len(args) == 2)
uv1 = a.miriad.UV(args[0])
uv2 = a.miriad.UV(args[1])

def mfunc(uv, p, d, f):
    p2,d2,f2 = uv2.read(raw=True)
    f = n.logical_or(f,f2)
    d = d - d2
    return p, d, f

uvo = a.miriad.UV(args[0] + 'd', status='new')
uvo.init_from_uv(uv1)
uvo.pipe(uv1, mfunc=mfunc, raw=True, 
    append2hist='DIFUV: subtracted file %s\n' % args[1])
