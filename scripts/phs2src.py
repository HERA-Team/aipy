#! /usr/bin/env python

import aipy.ants, numpy

def uv2freqs(uv):
    #sfreq = uv['sfreq']
    sfreq = 0.1126
    #sdf = uv['sdf']
    sdf = 0.00234 / 2
    return numpy.arange(uv['nchan'], dtype=numpy.float) * sdf + sfreq

def uv2aa(uv, antpos=None, delays=None):
    if antpos is None:
        antpos = uv['antpos']
        antpos.shape = (3, uv['nants'])
        antpos = antpos.transpose()
    else: a = antpos
    if delays is None:
        try: delays = uv['delay']
        except(KeyError):
            delays = numpy.zeros((uv['nants'],), dtype=numpy.float64)
    ants = []
    for n, d in enumerate(delays):
        print n, a[n], d
        ants.append(aipy.ants.Antenna(antpos[n,0], antpos[n,1], antpos[n,2], 
            delay=d))
    location = (uv['latitud'], uv['longitu'])
    freqs = uv2freqs(uv)
    return aipy.ants.AntennaArray(ants, location, freqs=freqs)

if __name__ == '__main__':
    import sys, aipy.miriad, os, aipy.srcs
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('phs2src.py [options] *.uv')
    p.set_description(__doc__)
    p.add_option('-s', '--source', dest='source',
        help='The source to phase to.  Options are sun and cyg.')
    p.add_option('-a', '--antpos', dest='antpos', action='store_true',
        help='Override antenna positions.')
    p.add_option('-d', '--delays', dest='delays', action='store_true',
        help='Override antenna delays.')
    p.add_option('-r', '--rmsrc', dest='rmsrc', action='store_true',
        help='Phase to src, remove, then unphase.')
    p.add_option('-w', '--swath', dest='swath', default=0, type='int',
        help='Number of bins around center to wipe out when removing a src.')
    opts, args = p.parse_args(sys.argv[1:])
    
    uv = aipy.miriad.UV(args[0])

    antpos = numpy.array([
        [  0.,     0.,     0.],
        [ -90.51,  142.55, -203.70], #[-98.,   138.,  -185.],
        [ -60.72,  377.92, -135.90], #[-70.,   382.,  -112.],
        [  56.77,  460.56,  117.51], #[ 58.,   463.,   114.],
    ])
    delays = [0., -7.93, -5.90, 3.09]

    if not opts.antpos: antpos = None
    if not opts.delays: delays = None
    aa = uv2aa(uv, antpos=antpos, delays=delays)
    del(uv)

    src = aipy.srcs.srcs[opts.source]

    # A pipe to use for phasing to a source
    def phs(uv, p, d):
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return p, d
        aa.set_jultime(p[-2])
        d = aa.phs2src(d, src, bl)
        return p, d

    # A pipe to use for removing a source
    def rm(uv, p, d):
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return p, d
        aa.set_jultime(p[-2])
        d = aa.rmsrc(d, src, bl, swath=opts.swath)
        return p, d

    if opts.rmsrc: f = rm
    else: f = phs

    for filename in args:
        print filename
        uvofile = filename + '.' + opts.source
        if os.path.exists(uvofile):
            print 'File exists: skipping'
            continue
        uvi = aipy.miriad.UV(filename)
        uvo = aipy.miriad.UV(uvofile, status='new')
        aipy.miriad.pipe_uv(uvi, uvo, mfunc=f, exclude=['bandpass', 'leakage'])
        del(uvi); del(uvo)

