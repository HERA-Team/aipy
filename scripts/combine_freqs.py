#! /usr/bin/env python
"""
A script for reducing the number of channels in a UV data set by coherently
adding adjacent channels together.

Author: Aaron Parsons
Date: 6/03/07
Revisions:
    12/11/07    arp Ported to new miriad file interface
"""

import aipy, sys, os, numpy
from optparse import OptionParser

p = OptionParser()
p.set_usage('combine_freqs.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-n', '--nchan', dest='nchan', default=256, type='int',
    help='Reduce the number of channels in a spectrum to this number.')
p.add_option('-c', '--careful_flag', dest='careful_flag', action='store_true',
    help='Flag resultant bin if any component bins are flagged (otherwise, flags only when there are more flagged bins than unflagged bins).')
p.add_option('-u', '--unify', dest='unify', action='store_true',
    help='Output to a single UV file.')
opts, args = p.parse_args(sys.argv[1:])

uvo = None
for uvfile in args:
    print uvfile
    uvi = aipy.miriad.UV(uvfile)
    sfreq = uvi['sfreq']
    sdf = uvi['sdf']
    nchan = uvi['nchan']
    newsfreq = sfreq + ((nchan/opts.nchan) * sdf) / 2
    newsdf = sdf * (nchan/opts.nchan)
    
    if uvo is None:
        uvofile = uvfile+'m'
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
        uvo = aipy.miriad.UV(uvofile, status='new')
        if nchan != opts.nchan:
            uvo.init_from_uv(uvi, override={'nchan':opts.nchan, 
                'sfreq':newsfreq, 'sdf':newsdf, 'nschan':opts.nchan, 
                'freq':newsfreq, 'nchan0':opts.nchan, },)
        else: uvo.init_from_uv(uvi)

    if nchan != opts.nchan:
        def f(uv, p, d):
            d.shape = (opts.nchan, nchan/opts.nchan)
            m = d.mask.sum(axis=1)
            if opts.careful_flag: m = numpy.where(m > 0, 1, 0)
            else: m = numpy.where(m >= nchan/opts.nchan/2, 1, 0)
            d = numpy.ma.average(d, axis=1)
            d = numpy.ma.array(d.data, mask=m, dtype=numpy.complex64)
            return p, d
        uvo.pipe(uvi, mfunc=f, append2hist='COMB_FREQ: nchan=%d, careful=%s, unify=%s\n' % (opts.nchan, opts.careful_flag, opts.unify))
    else:
        uvo.pipe(uvi, append2hist='Miniaturized...\n')
    if not opts.unify:
        del(uvo)
        uvo = None
    
