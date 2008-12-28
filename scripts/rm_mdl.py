#! /usr/bin/env python

import aipy.sim, aipy.srcs, aipy.eor, numpy

def uv2freqs(uv):
    #sfreq = uv['sfreq']
    sfreq = 0.1126
    #sdf = uv['sdf']
    sdf = 0.00234 / 2
    return numpy.arange(uv['nchan'], dtype=numpy.float) * sdf + sfreq

def uv2aa(uv):
    freqs = uv2freqs(uv)
    p = numpy.array([  3.89595979e+04,  -4.74699190e+04,   2.52022159e+04,
                      -7.61402543e+03,   1.43157780e+03,  -1.71509536e+02,
                       1.27843283e+01,  -5.42007429e-01,   1.00051341e-02])
    b = aipy.eor.EorBeam(freqs)
    antennas = [
        aipy.fit.FitSimAntenna(
            0., 0., 0., delay=0.,
            amp=1.00,
            beam=b, gain_poly=p),
        aipy.fit.FitSimAntenna(
            #-101.08, 138.81, -177.51, delay=-7.88,
            #-100.55, 138.74, -177.18, delay=-8.41,
            # -98.60, 138.72, -174.29, delay=-11.09,
             -100.25, 138.71, -175.54, delay=-9.19,
            amp=1.04,
            beam=b, gain_poly=p),
        aipy.fit.FitSimAntenna(
             #-63.85,  383.32, -130.94, delay=-6.08,
             #-63.84,  383.49, -130.33, delay=-6.15,
             #-65.30,  383.05, -131.19, delay=-4.67,
              -66.82,  383.26, -132.35, delay=-2.91,
            amp=1.02,
            beam=b, gain_poly=p),
        aipy.fit.FitSimAntenna(
              #64.46,  465.39,  100.15, delay=+1.94,
              #64.08,  465.44,  101.01, delay=+2.08,
              #64.35,  464.97,  100.63, delay=+1.62,
               60.94,  465.09,   98.27, delay=+5.61,
            amp=0.95,
            beam=b, gain_poly=p),
    ]
    location = (uv['latitud'], uv['longitu'])
    return aipy.sim.SimAntennaArray(antennas, location)

def uv2srcs(uv):
    freqs = uv2freqs(uv)
    srcs = []
    for s in aipy.srcs.srcs: 
        ra, dec, strength, meas_freq, index = aipy.srcs.srcs[s]
        if s == 'sun':
            src = aipy.sim.SimRadioSun(strength, freqs,
                spec_index=index, meas_freq=meas_freq)
        else:
            src = aipy.sim.SimRadioFixedBody(ra, dec, strength, freqs,
                spec_index=index, meas_freq=meas_freq, name=s)
        if s in ['sun', 'cyg']: srcs.append(src) 
        #srcs.append(src)
    srclist = aipy.sim.SourceList(srcs)
    return srclist

if __name__ == '__main__':
    import sys, aipy.miriad, os, aipy.srcs
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('phs2src.py [options] *.uv')
    p.set_description(__doc__)
    p.add_option('-s', '--sim', dest='sim', action='store_true',
        help='Output a simulated dataset (rather than subtracting).')
    opts, args = p.parse_args(sys.argv[1:])
    
    uv = aipy.miriad.UV(args[0])

    aa = uv2aa(uv)
    srclist = uv2srcs(uv)
    del(uv)

    # A pipe for just outputting the model
    def mdl(uv, p, d):
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return p, d
        aa.set_jultime(p[-2])
        #d -= aa.sim_data(srclist.sources, bl, uv['pol'])
        d = aa.sim_data(srclist.sources, bl, stokes=-6)
        d = numpy.ma.array(d, mask=numpy.zeros_like(d))
        return p, d

    # A pipe to use for removing the model
    def rm(uv, p, d):
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if i == j: return p, d
        aa.set_jultime(p[-2])
        #d -= aa.sim_data(srclist.sources, bl, uv['pol'])
        d -= aa.sim_data(srclist.sources, bl, stokes=-6)
        return p, d

    if opts.sim: f = mdl
    else: f = rm

    for filename in args:
        print filename
        uvofile = filename + 's'
        if os.path.exists(uvofile):
            print 'File exists: skipping'
            continue
        uvi = aipy.miriad.UV(filename)
        uvo = aipy.miriad.UV(uvofile, status='new')
        aipy.miriad.pipe_uv(uvi, uvo, mfunc=f, exclude=['bandpass', 'leakage'])
        del(uvi); del(uvo)

