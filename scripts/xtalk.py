#! /usr/bin/env python
"""
A script for removing stable crosstalk from UV files.

Author: Aaron Parsons
Date: 6/03/07
Revisions: None
"""

import aipy.miriad, numpy, os, sys
from optparse import OptionParser

p = OptionParser()
p.set_usage('xtalk.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-f', '--sky_max_fringe', dest='sky_max_fringe', 
    type=float, default=.1,
    help='The highest fringe rate (in Hz) expected from celestial sources.')
p.add_option('-a', '--cross_fringe_center', dest='cross_fringe_center',
    type=float, default=0,
    help='The center fringe rate for crosstalk.')
p.add_option('-b', '--cross_fringe_width', dest='cross_fringe_width',
    type=float, default=0,
    help='The width of the fringe rate for crosstalk.')
p.add_option('-d', '--sky_max_delay', dest='sky_max_delay', 
    type=float, default=600,
    help='The highest delay (in ns) expected from celestial sources.')
p.add_option('-x', '--cross_delay_center', dest='cross_delay_center', 
    type=float, default=12,
    help='The center delay (in ns) which is corrupted by xtalk.')
p.add_option('-y', '--cross_delay_width', dest='cross_delay_width', 
    type=float, default=140,
    help='The width of the delay (in ns) which is corrupted by xtalk.')

opts, args = p.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
inttime = uv['inttime']
sdf = uv['sdf']
nchan = uv['nchan']
bin_dly = 1 / (sdf * nchan)
max_dly_bin = numpy.round(opts.sky_max_delay / bin_dly)
d_center = numpy.round(opts.cross_delay_center / bin_dly)
d_width = numpy.round(opts.cross_delay_width / bin_dly)
d0 = d_center - d_width
d1 = d_center + d_width
print d0, d1, max_dly_bin
assert(d0 * d1 < 0)     # d0 and d1 have to be of opposite signs
assert(max_dly_bin > 0)
del(uv)

# Group files into 1/10 of a julian day for processing in batches.
groups = {}
for uvfile in args:
    uv = aipy.miriad.UV(uvfile)
    p, d = uv.read_data()
    t = p[-2]
    t = int(t*10) / 2
    if not groups.has_key(t): groups[t] = [uvfile]
    else: groups[t].append(uvfile)
    del(uv)

for g in groups:
    phs_off = {}
    cnt = {}
    files = groups[g]
    files.sort()
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
        uvi = aipy.miriad.UV(uvfile)
        uvi.select_data('auto', 0, 0, include_it=0)
        # Gather all data
        while True:
            p, d = uvi.read_data()
            if d.size == 0: break
            bl = p[-1]
            c = nchan - numpy.sum(d.mask)
            d = numpy.fft.ifft(d.filled(0))
            try:
                phs_off[bl].append(d)
                cnt[bl] += c
            except(KeyError):
                phs_off[bl] = [d]
                cnt[bl] = c
        del(uvi)

    bin_dly = 1. / (len(phs_off[phs_off.keys()[0]]) * inttime)
    fng_width = numpy.round(opts.cross_fringe_width / bin_dly)
    max_fng_bin = numpy.round(opts.sky_max_fringe / bin_dly)
    fng_center = numpy.round(opts.cross_fringe_center / bin_dly)
    f0 = fng_center - fng_width
    f1 = fng_center + fng_width
    print f0, f1, max_fng_bin
    assert(f0 * f1 >= 0)   # Ensure fw0 and fw1 have the same sign
    assert(max_fng_bin > 0)

    # Average data by channel over entire file
    for bl in phs_off:
        d = numpy.array(phs_off[bl])
        d = numpy.fft.fft(d, axis=0)
        d[:,max_dly_bin:-max_dly_bin] = 0
        d[max_fng_bin:-max_fng_bin] = 0
        if f0 == f1:
            d[f0,d0:] = 0
            d[f0,:d1] = 0
        else:
            d[f0:f1,d0:] = 0
            d[f0:f1,:d1] = 0
        d = numpy.fft.ifft(d, axis=0)
        d = numpy.fft.fft(d, axis=1)
        phs_off[bl] = d

    cnt = {}
    for bl in phs_off: cnt[bl] = 0

    # Generate a pipe for removing average phase bias from data
    def phs_corr_mfunc(uv, preamble, data):
        bl = preamble[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return preamble, data
        data = numpy.ma.array(phs_off[bl][cnt[bl],:], mask=data.mask)
        cnt[bl] += 1
        return preamble, data

    for uvfile in files:
        # Apply the pipe to the data
        print 'Working on', uvfile
        uvofile = uvfile+'x'
        uvi = aipy.miriad.UV(uvfile)
        uvo = aipy.miriad.UV(uvofile, status='new')
        aipy.miriad.pipe_uv(uvi, uvo, mfunc=phs_corr_mfunc)
        del(uvo)
        del(uvi)
