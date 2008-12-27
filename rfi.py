"""
Module for automatic detection and flagging of rfi.

Author: Aaron Parsons
Date: 11/10/2006
Revisions: 
    11/12/2006  arp     Added numpy.abs/numpy.absolute workaround.
"""
import numpy, miriad

# Add a workaround for some versions of numpy having 'absolute' instead of 'abs'
try: numpyabs = numpy.abs
except(AttributeError): numpyabs = numpy.absolute

def derivative(a, axis=0):
    """Take the derivative of 'a' along the specified axis."""
    o = range(len(a.shape))
    o[axis], o[0] = 0, axis
    a = a.transpose(o)
    if numpy.ma.isMaskedArray(a): rv = numpy.zeros_like(a.filled())
    else: rv = numpy.zeros_like(a)
    rv[:-1] = a[1:] - a[:-1]
    a = a.transpose(o)
    rv = rv.transpose(o)
    return rv

def range2list(*ranges):
    """Generate a list of values for specified range pairs, or single values.
    Includes endpoints."""
    L = []
    for r in ranges:
        if len(r) == 2: L += range(r[0], r[1] + 1)
        elif len(r) == 1: L.append(r[0])
        else:
            raise ValueError('Each range must have 1 or 2 values (%d given).' \
                % (len(r)))
    return L

def find_left_bounded_thresh(hist, min_flag_frac=.125, max_flag_frac=.5):
    cum, mid, thq = 0, 0, 0
    tot = hist.sum()
    for i, h in enumerate(hist):
        cum += h
        if mid == 0 and cum > max_flag_frac: mid = i
        if cum > 1 - min_flag_frac:
            thq = i
            break
    amax1 = numpy.argmax(h[:mid])
    amax2 = numpy.argmax(h[mid:]) + mid
    rmin = numpy.argmin(h[amax1:amax2]) + amax1
    rmin = max(rmin, mid)
    rmin = min(rmin, thq)
    return rmin

def gen_rfi_mask(data, also=[], bins=None, use_deriv=True):
    """Faster (and I think better) way of generating a per-frequency rfi mask.
    Don't know about robustness for non-PAPER array data."""
    if bins is None: bins = data.shape[-1] / 20
    adat = numpyabs(data)
    if use_deriv:
        d = derivative(adat)
        dd = numpyabs(derivative(d))
    else:
        dd = adat
    h, bvals = numpy.histogram(numpy.log10(dd+1), bins=bins)
    bvals = 10**bvals
    h[0] = 0
    tot = h.sum()
    cum = 0; mid = 0; thq = 0
    for i in range(bins):
        cum += h[i]
        if mid == 0 and cum > tot / 2: mid = i
        if thq == 0 and cum > tot * 7 / 8:
            thq = i
            break
    mid += 1
    try:
        amax1 = numpy.argmax(h[:mid])
        amax2 = numpy.argmax(h[mid:]) + mid
        rmin = numpy.argmin(h[amax1:amax2]) + amax1
        rmin = max(rmin, mid)
        rmin = min(rmin, thq)
        thresh = bvals[rmin]
    except: thresh = 0
    print amax1, amax2, rmin
    import pylab
    pylab.plot(h)
    pylab.show()
    m = numpy.where(dd > thresh, 1, 0)
    for a in also: m[a] = 1
    return m.astype(numpy.bool)

def gen_post_sim_mfunc(simulator, remask=False):
    dbuf = {}
    # Collect the data we want from all the uvfiles
    for uvf in simulator.uvfiles:
        uv = miriad.UV(uvf)
        uv.configure_preamble('time/baseline/pol')
        while True:
            preamble, data = uv.read_data()
            if data.size == 0: break
            t, bl, p = preamble
            if not simulator.is_active(bl, p): continue
            i, j = simulator.gen_ij(bl)
            data = data.take(simulator.chans)
            V_f = simulator.sim_data(bl, stokes=p)
            k = '%d-%d,%d' % (i,j,p)
            data -= V_f
            if dbuf.has_key(k): dbuf[k].append(data)
            else: dbuf[k] = [data]
        del(uv)
    # For each baseline,pol,channel generate an rfi mask for the time series
    for k in dbuf:
        d = numpy.array(dbuf[k])
        print d.shape
        for i in range(d.shape[1]):
            d[:,i] = gen_rfi_mask(d[:,i])
        dbuf[k] = d.astype(numpy.bool)
    # Create the mfunc which will mask data flagged by gen_rfi_mask
    def mfunc(uv, preamble, data):
        u, v, w, t, bl = preamble
        p = uv['pol']
        if not simulator.is_active(bl, p): return preamble, data
        i, j = simulator.gen_ij(bl)
        if i == j: return preamble, data
        k = '%d-%d,%d' % (i,j,p)
        m = dbuf[k][0]
        dbuf[k] = dbuf[k][1:]
        if not remask: m |= data.mask.take(simulator.chans)
        data.mask.put(simulator.chans, m)
        return preamble, data
    # Return this mfunc to do work elsewhere
    return mfunc

#def gen_freq_mfunc(also=[], bins=None, remask=False):
#    """Create an mfunc suitable for miriad.pipe_uv use which masks rfi frequency
#    bins based on premise that data + rfi is a bimodal distribution, and that 
#    there is a local minimum between them.  Will not flag more than 1/2 the 
#    data, or less than 1/4, unless the data is totally rfi dominated, in
#    which case it masks the entire spectrum.  Uses the union of all 
#    autocorrelations masks to mask all baselines (autos + crosses) for a
#    given integration.  This may be a little aggressive, but it works."""
#    def mfunc(uv, p, d):
#        # Generate an rfi mask which is the 'or' of all autocorr rfi masks.
#        m = numpy.zeros(dbuf[0].shape)
#        for i, p in enumerate(pbuf):
#            bl = int(p[-1])
#            if (bl>>8) & 255 == bl & 255:
#                m = numpy.logical_or(m, 
#                    gen_rfi_mask(dbuf[i].data, also, bins))
#        dout = []
#        for d in dbuf:
#            if remask: dout.append(numpy.ma.array(d.data, mask=m))
#            else:
#                new_m = numpy.logical_or(m, d.mask)
#                dout.append(numpy.ma.array(d.data, mask=new_m))
#        return pbuf, dout
#    return mfunc

def gen_rfi_int_thresh(pwrs, bins=10):
    """Given a list of total powers for many integrations, figures out
    a threshold for excising rfi-dominated integrations.  This threshold is
    halfway (in log space) between the pwr mode (most populated histogram bin)
    and the maximum recorded value."""
    h, bvals = numpy.histogram(numpy.log10(pwrs + 1), bins=bins)
    bvals = 10**bvals
    amax = (numpy.argmax(h) + (len(h) -1)) / 2
    thresh = bvals[amax]
    return thresh

def gen_int_mfunc(uvi, remask=False):
    """Create an mfunc suitable for miriad.pipe_uv use which masks rfi
    integrations based on total recorded power.  This works, but a better
    way may be to look at the hairiness (squiggles from FFT overflows) which
    ruin all the data in a spectrum and warrents throwing out the integration.
    """
    pwr_dict = {}
    try:
        while True:
            preamble, data = uvi.read_data()
            if data.size == 0: break
            bl = int(preamble[-1])
            pol = int(uvi['pol'])
            k = '%d,p%d' % (bl, pol)
            p = numpyabs(data.data).sum()
            try: pwr_dict[k].append(p)
            except(KeyError): pwr_dict[k] = [p]
    except(KeyboardInterrupt): pass
    thresh = {}
    for k in pwr_dict:
        pwr_dict[k] = numpy.array(pwr_dict[k])
        thresh[k] = gen_rfi_int_thresh(pwr_dict[k])
    uvi.rewind()
    def mfunc(pbuf, dbuf, vbuf, cbuf=None):
        dout = []
        for p, d, v in zip(pbuf, dbuf, vbuf):
            bl = int(p[-1])
            pol = int(v['pol'])
            k = '%d,p%d' % (bl, pol)
            pwr = numpyabs(data).sum()
            if pwr > thresh[k]: 
                dout.append(numpy.ma.array(d.data, mask=1))
            elif remask:
                dout.append(numpy.ma.array(d.data, mask=0))
            else: dout.append(d)
        return pbuf, dout
    return mfunc
        
#  ____            _       _     ___       _             __
# / ___|  ___ _ __(_)_ __ | |_  |_ _|_ __ | |_ ___ _ __ / _| __ _  ___ ___
# \___ \ / __| '__| | '_ \| __|  | || '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
#  ___) | (__| |  | | |_) | |_   | || | | | ||  __/ |  |  _| (_| | (_|  __/
# |____/ \___|_|  |_| .__/ \__| |___|_| |_|\__\___|_|  |_|  \__,_|\___\___|
#                   |_|

if __name__ == '__main__':
    import sys, os, miriad
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('rfi.py [options] *.uv')
    p.set_description(__doc__)
    p.add_option('-r', '--remask', dest='remask', action='store_true',
        help='Clobber any existing masks the data might have.')
    p.add_option('-f', '--freqonly', dest='freqonly', action='store_true',
        help='Flag only by frequency (and not by integration).')
    p.add_option('-i', '--intonly', dest='intonly', action='store_true',
        help='Flag only by integration (and not by frequency).')
    p.add_option('-s', '--silent', dest='silent', action='store_true',
        help="Don't say a word.")
    opts, args = p.parse_args(sys.argv[1:])

    if opts.silent:
        def output(s): pass
    else:
        def output(s): print(s)

    if not opts.intonly:
        also = range2list((0,69), (235,255), (183,194))
        freq_mfunc = gen_freq_mfunc(also, remask=opts.remask)

    for f in args:
        if os.path.exists(f + '.xrfi'):
            output('Skipping %s: xrfi file already exists.' % (f))
            continue
        output('Working on file: ' + f)
        uvi = miriad.UV(f)
        uvo = miriad.UV(f + '.xrfi', 'new')
        if not opts.freqonly:
            output('    Creating mfunc')
            int_mfunc = gen_int_mfunc(uvi, remask=opts.remask)
        if opts.freqonly:
            mfunc = freq_mfunc
            hist_lines = 'RFI: Flagged by frequency only.\n'
            output('    Flagging data by freq')
        elif opts.intonly:
            mfunc = int_mfunc
            hist_lines = 'RFI: Flagged by integration only.\n'
            output('    Flagging data by int')
        else:
            def mfunc(pbuf, dbuf, vbuf, cbuf):
                pbuf, dbuf = freq_mfunc(pbuf, dbuf, vbuf, cbuf)
                pbuf, dbuf = int_mfunc(pbuf, dbuf, vbuf, cbuf)
                return pbuf, dbuf
            hist_lines = 'RFI: Flagged by frequency and integration.\n'
            output('    Flagging data by freq and int')
        miriad.pipe_uv(uvi, uvo, mfunc, append2hist=hist_lines)
        del(uvi); del(uvo)
        output('    Done.')

