#!/usr/bin/python
"""
A script for filtering data in delay and fringe-rate domain.

Author: Aaron Parsons
Date: 10/10/07
Revisions:
    12/11/07 arp    Ported to use new miriad file interface
"""
import os, sys, aipy, numpy

__version__ = '0.0.1'

def gen_skypass_delay(aa, sdf, nchan, max_bl_frac=1.):
    bin_dly = 1. / (sdf * nchan)
    filters = {}
    for i,j in [aa.bl2ij(bl) for bl in aa.bl_order]:
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)
        max_bl = max_bl_frac * numpy.sqrt(numpy.dot(max_bl, max_bl))
        dly = aa.get_delay(i,j)
        uthresh, lthresh = (dly + max_bl)/bin_dly + 1, (dly - max_bl)/bin_dly
        uthresh, lthresh = int(numpy.round(uthresh)), int(numpy.round(lthresh))
        f = numpy.ones((nchan,), dtype=numpy.float)
        f[uthresh:lthresh] = 0
        filters[bl] = f
    return filters

def gen_skypass_fringe(aa, inttime, N_int, nchan, margin=1.2):
    freqs = aa.ants[0].beam.freqs
    bin_dly = 1. / (inttime * N_int)
    dth_dt = 2*numpy.pi / aipy.const.sidereal_day
    filters = {}
    for i,j in [aa.bl2ij(bl) for bl in aa.bl_order]:
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)[:2]
        max_bl = numpy.sqrt(numpy.dot(max_bl, max_bl)) * freqs
        max_fr = margin * dth_dt * max_bl
        uthresh, lthresh = max_fr / bin_dly + 1, -max_fr / bin_dly
        uthresh = numpy.round(uthresh).astype(numpy.int)
        lthresh = numpy.round(lthresh).astype(numpy.int)
        f = numpy.ones((N_int, nchan), dtype=numpy.float)
        for n in range(nchan):
            f[uthresh[n]:lthresh[n],n] = 0
        filters[bl] = f
    return filters
    
from optparse import OptionParser

p = OptionParser()
p.set_usage('xtalk2.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-l', '--loc', dest='loc', default='pwa303',
    help='Use location-specific info for this location (default pwa303).')
opts, args = p.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
nchan, inttime = uv['nchan'], uv['inttime']
sdf, sfreq = uv['sdf'], uv['sfreq']
aa = aipy.loc.get_aa(opts.loc, sdf, sfreq, nchan, use_bp=False)
skypass_delay = gen_skypass_delay(aa, sdf, nchan, max_bl_frac=1.)
del(uv)

for uvfile in args:
    print 'Working on', uvfile
    uvofile = uvfile+'X'
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = aipy.miriad.UV(uvfile)

    # Gather all data from file
    data = {}
    times = []
    for (uvw,t,(i,j)),d in uvi.all():
        bl = aa.ij2bl(i,j)
        if len(times) == 0 or times[-1] != t: times.append(t)
        if i != j:
            # Filter out delays that don't correspond to the sky
            d = numpy.fft.fft(numpy.fft.ifft(d.filled(0)) * skypass_delay[bl])
        try: data[bl].append(d)
        except(KeyError): data[bl] = [d]
    # Filter out fringes that don't correspond to the sky
    skypass_fringe = gen_skypass_fringe(aa, inttime, len(times), nchan)
    for bl in data:
        i, j = aa.bl2ij(bl)
        if i == j: continue
        d = numpy.array(data[bl])
        d = numpy.fft.fft(d, axis=0) * skypass_fringe[bl]
        # Also kill 0 fringe (narrows FOV, but kills a lot of xtalk and rfi)
        d[0] = 0
        d = numpy.fft.ifft(d, axis=0)
        data[bl] = d

    # Generate a pipe for swapping in new data
    def rfi_mfunc(uv, p, d):
        uvw, t, (i,j) = p
        bl = aipy.miriad.ij2bl(i,j)
        t = times.index(t)
        d = numpy.ma.array(data[bl][t], mask=d.mask)
        return p, d

    uvi.rewind()
    uvo = aipy.miriad.UV(uvofile, status='new')
    # Apply the pipe to the data
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=rfi_mfunc,
        append2hist='XTALK2: version %s\n' % __version__)
