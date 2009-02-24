#!/usr/bin/python
"""
Removes crosstalk from PAPER data by filtering delays and fringes that 
are unlikely to correspond to celestial sources.  Has shortcoming that
flagged data cause power to be scattered by filtering process.

Author: Aaron Parsons
"""
import os, sys, aipy as a, numpy as n

def gen_skypass_delay(aa, sdf, nchan, max_bl_frac=1.):
    bin_dly = 1. / (sdf * nchan)
    filters = {}
    for i,j in [aa.bl2ij(bl) for bl in aa.bl_order]:
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)
        max_bl = max_bl_frac * n.sqrt(n.dot(max_bl, max_bl))
        dly,off = aa.phsoff[-2:]
        uthresh, lthresh = (dly + max_bl)/bin_dly + 1, (dly - max_bl)/bin_dly
        uthresh, lthresh = int(n.round(uthresh)), int(n.round(lthresh))
        f = n.ones((nchan,), dtype=n.float)
        f[uthresh:lthresh] = 0
        filters[bl] = f
    return filters

def gen_skypass_fringe(aa, inttime, N_int, nchan, margin=1.2):
    freqs = aa[0].beam.freqs
    bin_dly = 1. / (inttime * N_int)
    dth_dt = 2*n.pi / a.const.sidereal_day
    filters = {}
    for i,j in [aa.bl2ij(bl) for bl in aa.bl_order]:
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)[:2]
        max_bl = n.sqrt(n.dot(max_bl, max_bl)) * freqs
        max_fr = margin * dth_dt * max_bl
        uthresh, lthresh = max_fr / bin_dly + 1, -max_fr / bin_dly
        uthresh = n.round(uthresh).astype(n.int)
        lthresh = n.round(lthresh).astype(n.int)
        f = n.ones((N_int, nchan), dtype=n.float)
        for c in range(nchan):
            f[uthresh[c]:lthresh[c],c] = 0
        filters[bl] = f
    return filters
    
if __name__ == '__main__':
    import optparse

    o = optparse.OptionParser()
    o.set_usage('xtalk2.py [options] *.uv')
    o.set_description(__doc__)
    a.scripting.add_standard_options(o, loc=True)
    o.add_option('--clean', dest='clean', type='float', default=1e-3,
        help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
    opts,args = o.parse_args(sys.argv[1:])

    # Parse command-line options
    uv = a.miriad.UV(args[0])
    nchan, inttime = uv['nchan'], uv['inttime']
    sdf, sfreq = uv['sdf'], uv['sfreq']
    aa = a.loc.get_aa(opts.loc, sdf, sfreq, nchan)
    skypass_delay = gen_skypass_delay(aa, sdf, nchan, max_bl_frac=1.)
    del(uv)

    for uvfile in args:
        print 'Working on', uvfile
        uvofile = uvfile+'X'
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
        uvi = a.miriad.UV(uvfile)

        # Gather all data from file
        data = {}
        print '    Performing delay transform and cleaning...'
        for (uvw,t,(i,j)),d,f in uvi.all(raw=True):
            bl = aa.ij2bl(i,j)
            if i == j: continue
            # Filter out delays that don't correspond to the sky
            flags = n.logical_not(f).astype(n.float)
            gain = n.sqrt(n.average(flags**2))
            ker = n.fft.ifft(flags)
            d = n.where(f, 0, d)
            d = n.fft.ifft(d)
            if not n.all(d == 0):
                d, info = a.deconv.clean1d(d, ker, tol=opts.clean)
                d += info['res'] / gain
            d *= skypass_delay[bl]
            d = n.fft.fft(d)
            try: data[bl].append(d)
            except(KeyError): data[bl] = [d]

        print '    Performing fringe rate transform and cleaning...'
        # Filter out fringes that don't correspond to the sky
        ntimes = len(data.values()[0])
        skypass_fringe = gen_skypass_fringe(aa, inttime, ntimes, nchan)
        for bl in data:
            d = n.array(data[bl])
            flags = n.where(d[:,0] != 0, 1., 0.)
            gain = n.sqrt(n.average(flags**2))
            ker = n.fft.ifft(flags)
            d = n.fft.ifft(d, axis=0)
            for chan in range(d.shape[1]):
                d[:,chan],info = a.deconv.clean1d(d[:,chan],ker,tol=opts.clean)
                d[:,chan] += info['res'] / gain
            d *= skypass_fringe[bl]
            # Also kill 0 fringe (narrows FOV, but kills a lot of xtalk and rfi)
            d[0] = 0
            d = n.fft.fft(d, axis=0)
            data[bl] = d

        cnt = {}
        for bl in data: cnt[bl] = 0

        # Generate a pipe for swapping in new data
        def mfunc(uv, p, d, f):
            uvw, t, (i,j) = p
            if i == j: return p,d,f
            bl = a.miriad.ij2bl(i,j)
            d = data[bl][cnt[bl],:]
            cnt[bl] += 1
            return p, d, f

        # Apply the pipe to the data
        uvi.rewind()
        uvo = a.miriad.UV(uvofile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=mfunc, raw=True,
            append2hist='XTALK2: clean=%f\n' % (opts.clean))
