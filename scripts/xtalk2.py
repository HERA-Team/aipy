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
        dly = aa.get_delay(i,j)
        uthresh, lthresh = (dly + max_bl)/bin_dly + 1, (dly - max_bl)/bin_dly
        uthresh, lthresh = int(n.round(uthresh)), int(n.round(lthresh))
        f = n.ones((nchan,), dtype=n.float)
        f[uthresh:lthresh] = 0
        filters[bl] = f
    return filters

def gen_skypass_fringe(aa, inttime, N_int, nchan, margin=1.2):
    freqs = aa.ants[0].beam.freqs
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
        for n in range(nchan):
            f[uthresh[n]:lthresh[n],n] = 0
        filters[bl] = f
    return filters
    
if __name__ == '__main__':
    import optparse

    o = optparse.OptionParser()
    o.set_usage('xtalk2.py [options] *.uv')
    o.set_description(__doc__)
    a.scripting.add_standard_options(o, loc=True)
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
        times = []
        for (uvw,t,(i,j)),d in uvi.all():
            bl = aa.ij2bl(i,j)
            if len(times) == 0 or times[-1] != t: times.append(t)
            if i != j:
                # Filter out delays that don't correspond to the sky
                d = n.fft.fft(n.fft.ifft(d.filled(0)) * skypass_delay[bl])
            try: data[bl].append(d)
            except(KeyError): data[bl] = [d]
        # Filter out fringes that don't correspond to the sky
        skypass_fringe = gen_skypass_fringe(aa, inttime, len(times), nchan)
        for bl in data:
            i, j = aa.bl2ij(bl)
            if i == j: continue
            d = n.array(data[bl])
            d = n.fft.fft(d, axis=0) * skypass_fringe[bl]
            # Also kill 0 fringe (narrows FOV, but kills a lot of xtalk and rfi)
            d[0] = 0
            d = n.fft.ifft(d, axis=0)
            data[bl] = d

        # Generate a pipe for swapping in new data
        def rfi_mfunc(uv, p, d, f):
            uvw, t, (i,j) = p
            bl = a.miriad.ij2bl(i,j)
            t = times.index(t)
            return p, data[bl][t], f

        uvi.rewind()
        uvo = a.miriad.UV(uvofile, status='new')
        # Apply the pipe to the data
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=rfi_mfunc, raw=True,
            append2hist='XTALK2: removed crosstalk\n')
