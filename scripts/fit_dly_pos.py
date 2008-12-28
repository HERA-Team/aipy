#! /usr/bin/env python

import aipy.fit, numpy

def uv2freqs(uv):
    #sfreq = uv['sfreq']
    sfreq = .1126
    #sdf = uv['sdf']
    sdf = .00234 / 2
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
        ants.append(aipy.fit.Antenna(antpos[n,0], antpos[n,1], antpos[n,2],
            b, delay=d))
    location = (uv['latitud'], uv['longitu'])
    return aipy.fit.AntennaArray(ants, location)

def calc_gain(phsdata):
    img = numpy.fft.ifft(phsdata.filled(0))
    return img * fit_window(img.shape[0])

def fit_window(width):
    w = numpy.arange(width)
    w = numpy.where(w > width/2, width - w, w)
    return 10.**(-w)

if __name__ == '__main__':
    import sys, aipy.miriad, os, aipy.src
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('fit_dly_pos.py [options] *.uv')
    p.set_description(__doc__)
    p.add_option('-o', '--override', dest='override', action='store_true',
        help='Override antenna positions and delays.')
    p.add_option('-w', '--swath', dest='swath', default=0, type='int',
        help='Number of bins around center to wipe out when removing a src.')
    p.add_option('-x', '--decimate', dest='decimate', default=1, type='int',
        help='Only use every Nth time step in fitting calculations.')
    opts, args = p.parse_args(sys.argv[1:])

    uv = aipy.miriad.UV(args[0])
    if opts.override:
        print 'Overriding antenna parameters in UV file.'
        freqs = uv2freqs(uv)
        b = aipy.ant.Beam(freqs)
        ants = [
            aipy.fit.Antenna(    0.,    0.,    0., b),
            aipy.fit.Antenna(-100.9, 138.8,-196.4, b, delay=-1.6, offset=+0.16),
            aipy.fit.Antenna( -68.8, 383.0,-131.0, b, delay=-1.4, offset=-0.20),
            aipy.fit.Antenna(  60.3, 464.8, 120.9, b, delay=-1.6, offset=+0.30),
        ]
        location = (-0.466646451963, 2.03621421241) # Boolardy
        aa = aipy.fit.AntennaArray(ants, location)
    else:
        print 'Using antenna parameters from UV file.'
        aa = uv2aa(uv)
    del(uv)

    cat = aipy.src.get_catalog(['Sun', 'cyg', 'crab', 'vir', 'cen'])
    start_prms = aa.get_params({
        1:['x','y','z','delay','offset'], 
        2:['x','y','z','delay','offset'], 
        3:['x','y','z','delay','offset'], 
    })
    prm_list, key_list = aipy.fit.flatten_prms(start_prms)
    first_fit = None    # Used to normalize fit values to the starting fit

    def fit_func(prms):
        global first_fit
        prms = aipy.fit.reconstruct_prms(prms, key_list)
        aipy.fit.print_params(prms)
        aa.set_params(prms)
        score = 0
        cnt = 0
        curtime = -1
        src_phs = {}
        #data = []
        for uvfile in args:
            print uvfile
            uvi = aipy.miriad.UV(uvfile)
            uvi.select_data('auto', 0, 0, include_it=False)
            while True:
                p, d = uvi.read_data()
                if d.size == 0: break
                t, bl = p[-2:]
                # Use only every Nth integration, if decimation is specified.
                if curtime != t:
                    curtime = t
                    cnt = (cnt + 1) % opts.decimate
                    for s in src_phs:
                        score += abs(src_phs[s]).sum()
                if cnt != 0: continue
                aa.set_jultime(t)
                cat.compute(aa)
                for s in cat:
                    ds = d.filled(0)
                    # Remove all other interfering sources from data
                    for rs in cat:
                        if rs != s: ds = aa.rmsrc(ds, cat[rs], bl)
                    ds = numpy.ma.array(ds, mask=d.mask)
                    try:
                        ds = aa.phs2src(ds, cat[s], bl)
                        src_phs[s] = src_phs.get(s,0) + calc_gain(ds)
                    except(aipy.ant.PointingError): pass
            del(uvi)
        score = 1 / score
        if first_fit is None: first_fit = score
        print score / first_fit
        #import pylab
        #pylab.clf()
        #data = numpy.array(data)
        #data = numpy.fft.ifft(data, axis=1)
        #pylab.imshow(numpy.abs(data), aspect='auto')
        #pylab.show()
        return score / first_fit

    import scipy.optimize
    scipy.optimize.fmin(fit_func, prm_list)
    print 'Done!'

