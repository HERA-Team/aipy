#! /usr/bin/env python

import aipy.ant, numpy

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

if __name__ == '__main__':
    import sys, aipy.miriad, os, aipy.src
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('phs2src.py [options] *.uv')
    p.set_description(__doc__)
    p.add_option('-s', '--source', dest='source',
        help='The source to phase to. Can be several sources separated by commas if -r is set.')
    p.add_option('-o', '--override', dest='override', action='store_true',
        help='Override antenna positions and delays.')
    p.add_option('-r', '--rmsrc', dest='rmsrc', action='store_true',
        help='Phase to src, remove, then unphase.')
    p.add_option('-w', '--swath', dest='swath', default=0, type='int',
        help='Number of bins around center to wipe out when removing a src.')
    opts, args = p.parse_args(sys.argv[1:])
    
    uv = aipy.miriad.UV(args[0])
    if opts.override:
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
    norm_bandpass = {}
    while True:
        p, d = uv.read_data()
        if d.size == 0: break
        i,j = aa.bl2ij(p[-1])
        if i == j:
            norm_bandpass[i] = norm_bandpass.get(i,0) + d.filled(0)
    del(uv)
    for i in norm_bandpass:
        norm_bandpass[i] /= norm_bandpass[i].sum()
        norm_bandpass[i] = numpy.sqrt(norm_bandpass[i])

    cat = aipy.src.get_catalog(opts.source.split(','), type='ant')
    
    # A pipe to use for phasing to a source
    def phs(uv, p, d):
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return p, d
        aa.set_jultime(p[-2])
        cat.compute(aa)
        d = aa.phs2src(d, cat.values()[0], bl)
        return p, d

    # A pipe to use for removing a source
    def rm(uv, p, d):
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return p, d
        aa.set_jultime(p[-2])
        cat.compute(aa)
        nbp = norm_bandpass[i] * norm_bandpass[j]
        d = aa.rmsrc(d, cat.values(), bl, swath=opts.swath, norm_bandpass=nbp)
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
        aipy.miriad.pipe_uv(uvi, uvo, mfunc=f)
        del(uvi); del(uvo)

