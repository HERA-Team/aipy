#! /usr/bin/env python
import aipy as a, numpy as np
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
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
o.add_option('--to_npz', 
    help='Instead of applying mask to data, store it as npz of this name.  May only be used along with --combine.')
o.add_option('--from_npz', 
    help='Apply mask to data from this npz file (generated with --to_npz).  May only be used along with --combine.')
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
if opts.to_npz or opts.from_npz: assert(opts.combine)


for uvfile in args:
    uvofile = uvfile+'R'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    if opts.from_npz:
        print '    Reading flags from', opts.from_npz
        m = np.load(opts.from_npz)
        mask = {'xx':{257:{}}} # Just use dummy values here to mimic structure of mask dictionary
        for cnt,t in enumerate(m['times']):
            mask['xx'][257][t] = m[str(cnt)]
    else:
        uvi = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uvi, opts.ant)
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
            d = np.array([data[pol][bl][t] for t in data_times])
            m = np.array([mask[pol][bl][t] for t in data_times])
            if opts.df != None:
                ddf = d[:,1:-1] - .5 * (d[:,:-2] + d[:,2:])
                ddf2 = np.abs(ddf)**2
                sig = np.sqrt(np.median(ddf2, axis=1))
                sig.shape = (sig.size,1)
                m[:,0] |= True; m[:,-1] |= True
                m[:,1:-1] |= np.where(ddf2/sig**2 > opts.df**2, True, False)
            if opts.dt != None:
                ddt = d[1:-1,:] - .5 * (d[:-2,:] + d[2:,:])
                ddt2 = np.abs(ddt)**2
                sig = np.sqrt(np.median(ddt2, axis=0))
                sig.shape = (1,sig.size)
                m[0,:] |= True; m[-1,:] |= True
                m[1:-1,:] |= np.where(ddt2/sig**2 > opts.dt**2, True, False)
            if opts.df == None and opts.dt == None:
                ad = np.abs(d)
                med = np.median(ad)
                sig = np.sqrt(np.median(np.abs(ad-med)**2))
                m |= np.where(ad > med + opts.nsig * sig, True, False)
            for i, t in enumerate(data_times): mask[pol][bl][t] |= m[i]
        if opts.combine:
            new_mask = {}
            for pol in mask:
              for bl in mask[pol]:
                for t in mask[pol][bl]:
                    new_mask[t] = new_mask.get(t,0)+mask[pol][bl][t].astype(np.int)
            for t in new_mask:
                m = np.where(new_mask[t] >= opts.thresh, 1, 0)
                for pol in mask:
                  for bl in mask[pol]:
                    mask[pol][bl][t] = m
        del(uvi)

    if opts.to_npz:
        print '    Writing flags to', opts.to_npz
        m = {}
        _m = mask.values()[0].values()[0]
        times = np.array(_m.keys())
        for cnt,t in enumerate(times): m[str(cnt)] = _m[t]
        m['times'] = times
        np.savez(opts.to_npz, **m)
    else:
        # Generate a pipe for applying the mask to data as it comes in.
        def rfi_mfunc(uv, preamble, data, flags):
            uvw, t, (i,j) = preamble
            bl = a.miriad.ij2bl(i,j)
            if opts.combine:
                try: m = mask.values()[0].values()[0][t]
                except(KeyError): m = np.ones_like(flags) # default to flagging
            else: m = mask[uv['pol']][bl][t]
            return preamble, np.where(m, 0, data), m

        uvi = a.miriad.UV(uvfile)
        uvo = a.miriad.UV(uvofile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=rfi_mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')



