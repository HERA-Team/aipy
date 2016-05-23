#! /usr/bin/env python
"""
A script for reducing the number of channels in a UV data set by coherently
adding adjacent channels together.

Author: Aaron Parsons
Date: 6/03/07
Revisions:
    12/11/07    arp Ported to new miriad file interface
"""

import aipy as a, sys, os, numpy as np, optparse

o = optparse.OptionParser()
o.set_usage('combine_freqs.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-n', '--nchan', dest='nchan', default=256, type='int',
    help='Reduce the number of channels in a spectrum to this number.')
o.add_option('-c', '--careful_flag', dest='careful_flag', action='store_true',
    help='Flag resultant bin if any component bins are flagged (otherwise, flags only when there are more flagged bins than unflagged bins).')
o.add_option('-d', '--dont_flag', dest='dont_flag', action='store_true',
    help='Only flag resultant bin if every component bin is flagged (otherwise, uses any data available).')
o.add_option('-u', '--unify', dest='unify', action='store_true',
    help='Output to a single UV file.')
opts, args = o.parse_args(sys.argv[1:])

uvo = None
for uvfile in args:
    print uvfile,'->',uvfile+'m'
    uvi = a.miriad.UV(uvfile)
    sfreq,sdf,nchan = uvi['sfreq'], uvi['sdf'], uvi['nchan']
    newsfreq = sfreq + ((nchan/opts.nchan) * sdf) / 2
    newsdf = sdf * (nchan/opts.nchan)
    
    if uvo is None:
        uvofile = uvfile+'m'
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
        uvo = a.miriad.UV(uvofile, status='new')
        if nchan != opts.nchan:
            uvo.init_from_uv(uvi, override={'nchan':opts.nchan, 
                'sfreq':newsfreq, 'sdf':newsdf, 'nschan':opts.nchan, 
                'freq':newsfreq, 'nchan0':opts.nchan, },)
        else: uvo.init_from_uv(uvi)

    if nchan != opts.nchan:
        def mfunc(uv, p, d, f):
            d = np.where(f, 0, d)
            d.shape = (opts.nchan, nchan/opts.nchan)
            d = d.sum(axis=1)
            f.shape = (opts.nchan, nchan/opts.nchan)
            f = np.logical_not(f).astype(np.int).sum(axis=1)
            d /= f.clip(1,np.Inf)
            if opts.careful_flag: f = np.where(f < nchan/opts.nchan, 1, 0)
            elif opts.dont_flag: f = np.where(f < 1, 1, 0)
            else: f = np.where(f <= nchan/opts.nchan/2, 1, 0)
            return p, np.where(f, 0, d), f
        uvo.pipe(uvi, mfunc=mfunc, raw=True,
            append2hist='COMB_FREQ: nchan=%d careful=%s dont=%s unify=%s\n' % \
                (opts.nchan, opts.careful_flag, opts.dont_flag, opts.unify))
    else:
        uvo.pipe(uvi, 
            append2hist='COMB_FREQ: nchan=%d careful=%s dont=%s unify=%s\n' % \
                (opts.nchan, opts.careful_flag, opts.dont_flag, opts.unify))
    if not opts.unify:
        del(uvo)
        uvo = None
    
