#! /usr/bin/env python
"""
Removes stable crosstalk from PAPER data.  Uses manually set parameters to
filter delays and fringes that are unlikely to correspond to celestial
sources.  Has same shortcomings as xtalk2.py, and doesn't do as good a
job of customizing filters to antenna distribution, but requires less a priori
knowledge of array.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, os, sys, optparse

o = optparse.OptionParser()
o.set_usage('xtalk1.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--sky_max_fringe', dest='sky_max_fringe', 
    type=float, default=.1,
    help='The highest fringe rate (in Hz) expected from celestial sources.')
o.add_option('--cross_fringe_center', dest='cross_fringe_center',
    type=float, default=0,
    help='The center fringe rate for crosstalk.')
o.add_option('--cross_fringe_width', dest='cross_fringe_width',
    type=float, default=0,
    help='The width of the fringe rate for crosstalk.')
o.add_option('--sky_max_delay', dest='sky_max_delay', 
    type=float, default=600,
    help='The highest delay (in ns) expected from celestial sources.')
o.add_option('--cross_delay_center', dest='cross_delay_center', 
    type=float, default=12,
    help='The center delay (in ns) which is corrupted by xtalk.')
o.add_option('--cross_delay_width', dest='cross_delay_width', 
    type=float, default=140,
    help='The width of the delay (in ns) which is corrupted by xtalk.')
opts, args = p.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
inttime,sdf,nchan = uv['inttime'], uv['sdf'], uv['nchan']
bin_dly = 1 / (sdf * nchan)
max_dly_bin = n.round(opts.sky_max_delay / bin_dly)
d_center = n.round(opts.cross_delay_center / bin_dly)
d_width = n.round(opts.cross_delay_width / bin_dly)
d0 = d_center - d_width
d1 = d_center + d_width
print d0, d1, max_dly_bin
assert(d0 * d1 < 0)     # d0 and d1 have to be of opposite signs
assert(max_dly_bin > 0)
del(uv)

# Group files into batches of 4
args.sort()
for g in range(0, len(args), 4):
    phs_off = {}
    cnt = {}
    files = args[g:g+4]
    print 'Processing group', files
    flag = 0
    for uvfile in files:
        uvofile = uvfile+'x'
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping group.'
            flag = 1
            break
    if flag: continue
    for uvfile in files:
        print 'Reading', uvfile
        uvi = a.miriad.UV(uvfile)
        uvi.select('auto', 0, 0, include=0)
        # Gather all data
        for (uvw,t,(i,j)),d in uvi.all():
            bl = a.miriad.ij2bl(i,j)
            c = nchan - n.sum(d.mask)
            d = n.fft.ifft(d.filled(0))
            try:
                phs_off[bl].append(d)
                cnt[bl] += c
            except(KeyError):
                phs_off[bl] = [d]
                cnt[bl] = c
        del(uvi)

    bin_dly = 1. / (len(phs_off[phs_off.keys()[0]]) * inttime)
    fng_width = n.round(opts.cross_fringe_width / bin_dly)
    max_fng_bin = n.round(opts.sky_max_fringe / bin_dly)
    fng_center = n.round(opts.cross_fringe_center / bin_dly)
    f0 = fng_center - fng_width
    f1 = fng_center + fng_width
    print f0, f1, max_fng_bin
    assert(f0 * f1 >= 0)   # Ensure fw0 and fw1 have the same sign
    assert(max_fng_bin > 0)

    # Average data by channel over entire file
    for bl in phs_off:
        d = n.array(phs_off[bl])
        d = n.fft.fft(d, axis=0)
        d[:,max_dly_bin:-max_dly_bin] = 0
        d[max_fng_bin:-max_fng_bin] = 0
        if f0 == f1:
            d[f0,d0:] = 0
            d[f0,:d1] = 0
        else:
            d[f0:f1,d0:] = 0
            d[f0:f1,:d1] = 0
        d = n.fft.ifft(d, axis=0)
        d = n.fft.fft(d, axis=1)
        phs_off[bl] = d

    cnt = {}
    for bl in phs_off: cnt[bl] = 0

    # Generate a pipe for removing average phase bias from data
    def phs_corr_mfunc(uv, p, d, f):
        uvw,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        if i == j: return p, d, f
        cnt[bl] += 1
        return p, phs_off[bl][cnt[bl],:], f

    for uvfile in files:
        # Apply the pipe to the data
        print 'Working on', uvfile
        uvofile = uvfile+'x'
        uvi = a.miriad.UV(uvfile)
        uvo = a.miriad.UV(uvofile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=phs_corr_mfunc, raw=True,
            append2hist='XTALK1: removed crosstalk\n')
