#! /usr/bin/env python
"""
A script for removing a source using a fringe-rate/delay transform.

Author: Aaron Parsons
Date: 6/03/07
Revisions:
    12/11/07 arp    Ported to use new miriad file interface
"""

import aipy, numpy, os, sys
from optparse import OptionParser

p = OptionParser()
p.set_usage('rm_src.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-f', '--fng_w', dest='fng_w', 
    type=int, default=1,
    help='The number of fringe bins to null.')
p.add_option('-d', '--dly_w', dest='dly_w',
    type=float, default=5,
    help='The number of delay bins to null.')
p.add_option('-s', '--src', dest='src',
    help='The source to remove.')
p.add_option('-e', '--extract', dest='extract', action='store_true',
    help='Extract the source instead of removing it.')
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
opts, args = p.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
aa = aipy.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

src = aipy.src.get_src(opts.src, type='ant')

for uvfile in args:
    print 'Working on', uvfile
    phs_dat = {}
    cnt = {}
    uvofile = uvfile+'.'+opts.src
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = aipy.miriad.UV(uvfile)

    # Gather all data
    for (uvw,t,(i,j)),d in uvi.all():
        if i == j: continue
        aa.set_jultime(t)
        src.compute(aa)
        bl = aa.ij2bl(i,j)
        try:
            d = aa.phs2src(d, src, i, j)
            cnt[bl] = cnt.get(bl, 0) + d.mask.sum()
            d = numpy.fft.ifft(d.filled(0))
        except(aipy.ant.PointingError): d = numpy.zeros_like(d.data)
        try: phs_dat[bl].append(d)
        except(KeyError): phs_dat[bl] = [d]

    # Perform fringe rate transform
    for bl in phs_dat:
        d = numpy.array(phs_dat[bl])
        d /= float(d.size - cnt.get(bl,0)) / d.size
        d = numpy.fft.fft(d, axis=0)
        x1, x2 = opts.fng_w, -opts.fng_w+1
        if x2 == 0: x2 = d.shape[0]
        y1, y2 = opts.dly_w, -opts.dly_w
        if y2 == 0: y2 = d.shape[1]
        #print x1, x2, y1, y2
        #import pylab
        #pylab.imshow(numpy.log10(numpy.abs(d)))
        #pylab.show()
        d[x1:x2,:] = 0
        d[:,y1:y2] = 0
        d = numpy.fft.ifft(d, axis=0)
        d = numpy.fft.fft(d, axis=1)
        phs_dat[bl] = d

    cnt = {}
    for bl in phs_dat: cnt[bl] = 0

    # Generate a pipe for removing average phase bias from data
    def rm_mfunc(uv, p, d):
        uvw, t, (i,j) = p
        if i == j: return p, d
        bl = aa.ij2bl(i,j)
        aa.set_jultime(t)
        src.compute(aa)
        data = phs_dat[bl][cnt[bl],:]
        try: data = aa.unphs2src(data, src, i, j)
        except(aipy.ant.PointingError):
            if opts.extract: d *= 0
            data = 0
        cnt[bl] += 1
        if opts.extract: return p, numpy.ma.array(data, mask=d.mask)
        else: return p, d - data

    # Apply the pipe to the data
    uvi.rewind()
    uvo = aipy.miriad.UV(uvofile, status='new')
    # Apply the pipe to the data
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=rm_mfunc)
    del(uvi); del(uvo)
