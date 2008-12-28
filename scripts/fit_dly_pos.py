#! /usr/bin/env python

import aipy.fit, numpy

def uv2freqs(uv):
    #sfreq = uv['sfreq']
    sfreq = .1126
    #sdf = uv['sdf']
    sdf = .00234 / 2
    return numpy.arange(uv['nchan'], dtype=numpy.float) * sdf + sfreq

def uv2aa(uv, antpos=None, delays=None):
    if antpos is None:
        antpos = uv['antpos']
        antpos.shape = (3, uv['nants'])
        antpos = antpos.transpose()
    if delays is None:
        try: delays = uv['delay']
        except(KeyError):
            delays = numpy.zeros((uv['nants'],), dtype=numpy.float64)
    ants = []
    for n, d in enumerate(delays):
        print n, antpos[n], d
        ants.append(aipy.fit.FitAntenna(antpos[n,0], antpos[n,1], antpos[n,2], 
            delay=d))
    location = (uv['latitud'], uv['longitu'])
    freqs = uv2freqs(uv)
    return aipy.fit.FitAntennaArray(ants, location, freqs=freqs)

def aggressive_rmsrc(data, src, aa, i, j=None, swath=2):
    phsdata = aa.phs2src(data, src, i, j)
    if numpy.ma.isMA(phsdata): phsdata = phsdata.filled(0)
    img = numpy.fft.ifft(phsdata)
    img[:swath+1] = 0; img[-swath:] = 0
    phsdata = numpy.fft.fft(img)
    d = aa.unphs2src(phsdata, src, i, j)
    return d

def calc_gain(phsdata):
    img = numpy.abs(numpy.fft.ifft(phsdata.filled(0)))
    return numpy.sum(img * fit_window(img.shape[0]))

def fit_window(width):
    w = numpy.arange(width)
    w = numpy.where(w > width/2, width - w, w)
    return 10.**(-w)

if __name__ == '__main__':
    import sys, aipy.miriad, aipy.ants, aipy.srcs

    uv = aipy.miriad.UV(sys.argv[-1])

    antpos = numpy.array([
        [  0.,     0.,     0.],
        [ -90.51,  142.55, -203.70], #[-98.,   138.,  -185.],
        [ -60.72,  377.92, -135.90], #[-70.,   382.,  -112.],
        [  56.77,  460.56,  117.51], #[ 58.,   463.,   114.],
    ])
    delays = [0., -7.93, -5.90, 3.09]
    aa = uv2aa(uv, antpos=antpos, delays=delays)
    del(uv)

    srcs = ['sun', 'cyg', 'cen', 'vir']
    #start_prms = aa.get_params({1:'*', 2:'*', 3:'*'})
    start_prms = aa.get_params({
        1:['delay','pos'], 
        2:['delay','pos'], 
        3:['delay','pos']
    })
    prm_list, key_list = aipy.fit.flatten_prms(start_prms)

    def fit_func(args):
        prms = aipy.fit.reconstruct_prms(args, key_list)
        aipy.fit.print_params(prms)
        aa.set_params(prms)
        rv = 0
        data = []
        for uvfile in sys.argv[1:]:
            print uvfile,
            not_up = {}
            uvi = aipy.miriad.UV(uvfile)
            uvi.select_data('auto', 0, 0, include_it=False)
            #uvi.select_data('antennae', 1, 3, include_it=True)
            #uvi.select_data('antennae', 1, 4, include_it=True)
            #uvi.select_data('antennae', 3, 4, include_it=True)
            while True:
                p, d = uvi.read_data()
                if d.size == 0: break
                bl = p[-1]
                aa.set_jultime(p[-2])
                for s in srcs:
                    try:
                        phs_d = aa.phs2src(d, aipy.srcs.srcs[s], bl)
                        rv += calc_gain(phs_d)
                    except(aipy.ants.PointingError):
                        not_up[s] = 0
            not_up = not_up.keys()
            if len(not_up) == 0: print ''
            else: print 'not up:', not_up
            del(uvi)
        print rv
        #import pylab
        #pylab.clf()
        #data = numpy.array(data)
        #data = numpy.fft.ifft(data, axis=1)
        #pylab.imshow(numpy.abs(data), aspect='auto')
        #pylab.show()
        return -rv

    import scipy.optimize
    scipy.optimize.fmin(fit_func, prm_list)
    print 'Done!'

