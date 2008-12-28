#! /usr/bin/env python
"""
A script for removing stable crosstalk from UV files.

Author: Aaron Parsons
Date: 6/03/07
Revisions: None
"""

import aipy.miriad, numpy, os, sys, aipy.ant, aipy.src
from optparse import OptionParser

def uv2freqs(uv):
    #sfreq = uv['sfreq']
    sfreq = 0.1126
    #sdf = uv['sdf']
    sdf = 0.00234 / 2
    return numpy.arange(uv['nchan'], dtype=numpy.float) * sdf + sfreq

def uv2aa(uv):
    freqs = uv2freqs(uv)
    antpos = uv['antpos']
    antpos.shape = (3, uv['nants'])
    antpos = antpos.transpose()
    try: delays = uv['delay']
    except(KeyError):
        delays = numpy.zeros((uv['nants'],), dtype=numpy.float64)
    ants = []
    for n, d in enumerate(delays):
        b = aipy.ant.Beam(freqs)
        ants.append(aipy.ant.Antenna(antpos[n,0], antpos[n,1], antpos[n,2],
            b, delay=d))
    location = (uv['latitud'], uv['longitu'])
    return aipy.ant.AntennaArray(ants, location)


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
p.add_option('-o', '--override', dest='override', action='store_true',
    help='Override antenna positions and delays.')


opt, args = p.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
if opt.override:
    print 'Overriding antenna parameters in UV file.'
    freqs = uv2freqs(uv)
    b = aipy.ant.Beam(freqs)
    ants = [
        aipy.ant.Antenna(    0.,    0.,    0.,b,delay=+0.0, offset=+0.0),
        aipy.ant.Antenna(-100.7, 139.1,-197.7,b,delay=-1.62,offset=+0.16),
        aipy.ant.Antenna( -68.7, 383.4,-132.6,b,delay=-1.41,offset=-0.20),
        aipy.ant.Antenna(  59.4, 465.3, 120.7,b,delay=-1.67,offset=+0.28),
    ]
    location = (-0.466646451963, 2.03621421241) # Boolardy
    aa = aipy.ant.AntennaArray(ants, location)
else:
    print 'Using antenna parameters from UV file.'
    aa = uv2aa(uv)

src = aipy.src.get_src(opt.src, type='ant')

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
    phs_dat = {}
    cnt = {}
    files = groups[g]
    files.sort()
    print 'Processing group', files
    flag = 0
    for uvfile in files:
        uvofile = uvfile+'.'+opt.src
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
            aa.set_jultime(p[-2])
            src.compute(aa)
            try:
                d = aa.phs2src(d, src, bl)
                cnt[bl] = cnt.get(bl, 0) + d.mask.sum()
                d = numpy.fft.ifft(d.filled(0))
            except(aipy.ant.PointingError): d = numpy.zeros_like(d.data)
            try: phs_dat[bl].append(d)
            except(KeyError): phs_dat[bl] = [d]
        del(uvi)

    # Perform fringe rate transform
    for bl in phs_dat:
        d = numpy.array(phs_dat[bl])
        d /= float(d.size - cnt.get(bl,0)) / d.size
        d = numpy.fft.fft(d, axis=0)
        x1, x2 = opt.fng_w, -opt.fng_w+1
        if x2 == 0: x2 = d.shape[0]
        y1, y2 = opt.dly_w, -opt.dly_w
        if y2 == 0: y2 = d.shape[1]
        d[x1:x2,:] = 0
        d[:,y1:y2] = 0
        d = numpy.fft.ifft(d, axis=0)
        d = numpy.fft.fft(d, axis=1)
        phs_dat[bl] = d

    cnt = {}
    for bl in phs_dat: cnt[bl] = 0

    # Generate a pipe for removing average phase bias from data
    def rm_mfunc(uv, p, d):
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return p, d
        aa.set_jultime(p[-2])
        src.compute(aa)
        data = phs_dat[bl][cnt[bl],:]
        try: data = aa.unphs2src(data, src, bl)
        except(aipy.ant.PointingError): data = 0
        cnt[bl] += 1
        return p, d - data

    for uvfile in files:
        # Apply the pipe to the data
        print 'Working on', uvfile
        uvofile = uvfile+'.'+opt.src
        uvi = aipy.miriad.UV(uvfile)
        uvo = aipy.miriad.UV(uvofile, status='new')
        aipy.miriad.pipe_uv(uvi, uvo, mfunc=rm_mfunc)
        del(uvo)
        del(uvi)
