#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, os

o = optparse.OptionParser()
o.set_usage('xrfi_simple.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-c', '--chan', dest='chan',
    help='Manually flag channels before xrfi processing.  Options are "<chan1 #>,..." (a list of channels), or "<chan1 #>_<chan2 #>" (a range of channels).  Default is None.')
o.add_option('-n', '--nsig', dest='nsig', default=2., type='float',
    help='Number of standard deviations above mean to flag if neither --dt nor --df are specified.  Default 2.')
o.add_option('--df', dest='df', type='float', 
    help='Number of standard deviations above mean to flag, after taking derivative of frequency axis')
o.add_option('--dt', dest='dt', type='float',
    help='Number of standard deviations above mean to flag, after taking derivative of time axis')
o.add_option('--combine', dest='combine', action='store_true',
    help='Use the same mask for all baselines/pols (and use thresh to decide how many concidences it takes to flag all data.')
o.add_option('-t', '--thresh', dest='thresh', default=1, type='int',
    help='Number of flagging coincidences (baselines/pols) required to flag a time/chan.')
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
if not opts.chan is None:
    chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
else:
    chans = []
    opts.chan = 'None'
del(uv)


for uvfile in args:
    uvofile = uvfile+'R'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)
    #uvi.select('auto', -1, -1, include=False)
    # Gather all data and each time step
    data,mask,times = {}, {}, []
    for (uvw,t,(i,j)), d, f in uvi.all(raw=True):
        if len(times) == 0 or times[-1] != t: times.append(t)
        bl = a.miriad.ij2bl(i,j)
        pol = uvi['pol']
        if not pol in data:
            data[pol] = {}
            mask[pol] = {}
        if not bl in data[pol]:
            data[pol][bl] = {}
            mask[pol][bl] = {}
        # Manually flag channels
        f[chans] = 1
        mask[pol][bl][t] = f
        data[pol][bl][t] = d

    # Generate statistical mask
    for pol in data:
      for bl in data[pol]:
        i, j = a.miriad.bl2ij(bl)
        data_times = data[pol][bl].keys()
        data_times.sort()
        d = n.array([data[pol][bl][t] for t in data_times])
        m = n.array([mask[pol][bl][t] for t in data_times])
        if opts.df != None:
            ddf = d[:,1:-1] - .5 * (d[:,:-2] + d[:,2:])
            ddf2 = n.abs(ddf)**2
            sig = n.sqrt(n.median(ddf2, axis=1))
            sig.shape = (sig.size,1)
            m[:,0] |= 1; m[:,-1] |= 1
            m[:,1:-1] |= n.where(ddf2/sig**2 > opts.df**2, 1, 0)
        if opts.dt != None:
            ddt = d[1:-1,:] - .5 * (d[:-2,:] + d[2:,:])
            ddt2 = n.abs(ddt)**2
            sig = n.sqrt(n.median(ddt2, axis=0))
            sig.shape = (1,sig.size)
            m[0,:] |= 1; m[-1,:] |= 1
            m[1:-1,:] |= n.where(ddt2/sig**2 > opts.dt**2, 1, 0)
        if opts.df == None and opts.dt == None:
            ad = n.abs(d)
            med = n.median(ad)
            sig = n.sqrt(n.median(n.abs(ad-med)**2))
            m |= n.where(ad > med + opts.nsig * sig, 1, 0)
        for i, t in enumerate(data_times): mask[pol][bl][t] |= m[i]
    if opts.combine:
        new_mask = {}
        for pol in mask:
          for bl in mask[pol]:
            for t in mask[pol][bl]:
                new_mask[t] = new_mask.get(t,0)+mask[pol][bl][t].astype(n.int)
        for t in new_mask:
            m = n.where(new_mask[t] > opts.thresh, 1, 0)
            for pol in mask:
              for bl in mask[pol]:
                mask[pol][bl][t] = m

    # Generate a pipe for applying the mask to data as it comes in.
    def rfi_mfunc(uv, preamble, data, flags):
        uvw, t, (i,j) = preamble
        bl = a.miriad.ij2bl(i,j)
        m = mask[uv['pol']][bl][t]
        return preamble, n.where(m, 0, data), m

    del(uvi)
    uvi = a.miriad.UV(uvfile)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=rfi_mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')



