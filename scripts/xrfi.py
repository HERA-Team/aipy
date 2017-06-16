#!/usr/bin/python
"""
Detect and flag RFI related effects in UV files.  Uses statistical thresholding
to identify outliers in power and narrowness of features.  Does an improved
job of identifying low-level interference if sky model and crosstalk are
removed from the data first.
"""

import numpy as np, aipy as a, os, sys, pickle, optparse

o = optparse.OptionParser()
o.set_usage('xrfi.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-c', '--chan', dest='chan', 
    help='Manually flag channels before xrfi processing.  Options are "<chan1 #>,..." (a list of channels), or "<chan1 #>_<chan2 #>" (a range of channels).  Default is None.')
o.add_option('-n', '--nsig', dest='nsig', default=2., type='float',
    help='Number of standard deviations above mean to flag.  Default 2.')
o.add_option('-m', '--flagmode', dest='flagmode', default='both',
    help='Can be val,int,both,none for flagging by value only, integration only, both, or only manually flagged channels.  Default both.')
o.add_option('-s', '--share', dest='share', action='store_true',
    help='Flag a channel, integration if any pol/baseline flags it (share flags).')
o.add_option('--ch_thresh', dest='ch_thresh',type='float',default=.33,
    help='Fraction of the data in a channel which, if flagged, will result in the entire channel being flagged.  Default .33')
o.add_option('--int_thresh', dest='int_thresh',type='float',default=.99,
    help='Fraction of the data in an integration which, if flagged, will result in the entire integration being flagged.  Default .99')
o.add_option('--raw', dest='raw', action='store_true',
    help='Flag by integration without removing a smooth function.')
o.add_option('--reflag', dest='reflag', action='store_true',
    help='Ignore any previous flagging.')
opts, args = o.parse_args(sys.argv[1:])

def share_mask(mask):
    '''Flag a chan,integration if any pol/baseline has flagged it.'''
    m = {}
    for pol in mask:
      for bl in mask[pol]:
        for t in mask[pol][bl]:
            m[t] = m.get(t,0) | mask[pol][bl]
    for pol in mask:
      for bl in mask[pol]:
        for t in mask[pol][bl]:
            mask[pol][bl][t] = m[t]
    return mask

# Parse command-line options
uv = a.miriad.UV(args[0])
if not opts.chan is None:
    chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
else:
    chans = []
    opts.chan = 'None'
del(uv)

for uvfile in args:
    uvofile = uvfile+'r'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)
    (uvw,jd,(i,j)),d,f = uvi.read(raw=True)
    uvi.rewind()
    # Gather all data and each time step
    window = None
    data,mask,times = {}, {}, []
    for (uvw,t,(i,j)), d, f in uvi.all(raw=True):
        if len(times) == 0 or times[-1] != t: times.append(t)
        pol = uvi['pol']
        if not pol in data:
            data[pol] = {}
            mask[pol] = {}
        bl = a.miriad.ij2bl(i,j)
        if not bl in data[pol]:
            data[pol][bl] = {}
            mask[pol][bl] = {}
        if opts.reflag: f = np.zeros_like(f)
        f = np.where(np.abs(d) == 0, 1, f)
        if i != j:
            if window is None: 
                window = d.size/2 - abs(np.arange(d.size) - d.size/2)
            d = np.fft.fft(np.fft.ifft(d) * window)
        # Manually flagged channels
        f[chans] = 1
        mask[pol][bl][t] = f
        data[pol][bl][t] = d

    # Generate a single mask for all baselines which masks data if any
    # baseline has an outlier at that freq/time.  Data flagged
    # strictly on basis of nsigma above mean.
    if opts.flagmode in ['val', 'both']:
        for pol in data:
          for bl in data[pol]:
            i, j = a.miriad.bl2ij(bl)
            if i == j: continue
            data_times = data[pol][bl].keys()
            d = np.ma.array([data[pol][bl][t] for t in data_times],
                mask=[mask[pol][bl][t] for t in data_times])
            hi_thr, lo_thr = a.rfi.gen_rfi_thresh(d, nsig=opts.nsig)
            m = np.where(np.abs(d) > hi_thr,1,0)
            for i, t in enumerate(data_times): mask[pol][bl][t] |= m[i]
            # If more than ch_thresh of the data in a channel or 
            # integration is flagged, flag the whole thing
            ch_cnt = np.array([mask[pol][bl][t] for t in data_times]).sum(axis=0)
            ch_msk = np.where(ch_cnt > ch_cnt.max()*opts.ch_thresh,1,0)
            for t in mask[pol][bl]:
                if mask[pol][bl][t].sum() > mask[pol][bl][t].size \
                        * opts.int_thresh:
                    mask[pol][bl][t] |= 1
                else: mask[pol][bl][t] |= ch_msk

    # Use autocorrelations to flag entire integrations which have
    # anomalous powers.  At least 2 antennas must agree for a 
    # integration to get flagged.
    if opts.flagmode in ['int', 'both']:
        for pol in data:
          for bl in data[pol]:
            i, j = a.miriad.bl2ij(bl)
            if i != j: continue
            data_times = data[pol][bl].keys()
            data_times.sort()
            d = np.ma.array([data[pol][bl][t] for t in data_times],
                mask=[mask[pol][bl][t] for t in data_times])
            bad_ints = a.rfi.flag_by_int(d, nsig=opts.nsig, raw=opts.raw)
            for i in np.where(bad_ints)[0]:
                t = data_times[i]
                mask[pol][bl][t] |= 1

    if opts.share: mask = share_mask(mask)
    # Generate a pipe for applying the mask to data as it comes in.
    def rfi_mfunc(uv, preamble, data, flags):
        uvw, t, (i,j) = preamble
        f = mask[uv['pol']][a.miriad.ij2bl(i,j)][t]
        return preamble, np.where(f, 0, data), f

    uvi.rewind()
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=rfi_mfunc, raw=True, append2hist='XRFI: nsig=%f chans=%s mode=%s ch_thresh=%f int_thresh=%f raw=%s reflag=%s share=%s\n' %  (opts.nsig, opts.chan, opts.flagmode, opts.ch_thresh, opts.int_thresh, opts.raw, opts.reflag, opts.share))
