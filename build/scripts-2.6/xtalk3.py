#!/usr/global/paper/bin/python
"""
Removes crosstalk from PAPER data by modeling it over the course of a day
as a static per-channel cross-coupling that rises and falls uniformly across
the band with changing input amplitude.  Steps for crosstalk removal are:
run "xtalk3.py -o" on 1 day of UV files (which should be flagged for RFI, and
if possible, have a sky model removed) to 
generate *.xtalk files.  Then run "xtalk3.py -r" on same UV files to 
reprocess *.xtalk files, separating them into static shape/uniform gain 
terms that are stored in *.xtalk.rep files.  Finally, run "xtalk3.py -i" on 
UV files from the same JD (but which need not have RFI flagged or a sky model
removed) to use the *.xtalk.rep model to remove crosstalk from the visibility 
data.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, sys, os, optparse, pickle

o = optparse.OptionParser()
o.set_usage('xtalk3.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-i', '--infile', dest='infile', action='store_true',
    help='Apply xtalk calibrations generated with the -o option.')
o.add_option('-o', '--outfile', dest='outfile', action='store_true',
    help='Rather than apply the calibrations to the data, store them in a file (named by JD) to apply to a different file with the same JD.')
o.add_option('-r', '--reprocess', dest='reprocess', action='store_true',
    help='Reprocess xtalk files for the specified UV files to generate an average crosstalk profile and a time-dependent gain for scaling that profile.  The results will be saved to *.xtalk.rep files which will have precedence over *.xtalk files.')
o.add_option('-c', '--chan', dest='chan', default='160_720',
    help='Channels range to use for normalization when reprocessing.')
opts,args = o.parse_args(sys.argv[1:])

assert(not opts.infile or not opts.outfile)

if opts.reprocess:
    chans = map(int, opts.chan.split('_'))
    xsum, cnt, gain, times = {}, {}, {}, []
    for filename in args:
        uv = a.miriad.UV(filename)
        (uvw,jd,(i,j)),d,f = uv.read(raw=True)
        xfile = '%f.xtalk' % jd
        print xfile
        times.append(jd)
        if not os.path.exists(xfile):
            print xfile, 'does not exist.  Skipping...'
            continue
        f = open(xfile); xtalk = pickle.load(f); f.close()
        for bl in xtalk:
            dat = n.array(xtalk[bl])
            adat = n.ma.masked_equal(n.abs(dat), 0)
            if not gain.has_key(bl): gain[bl] = []
            gain[bl].append(n.average(adat[chans[0]:chans[1]]))
            dat /= gain[bl][-1]
            xsum[bl] = xsum.get(bl, 0) + dat
            cnt[bl] = cnt.get(bl, 0) + n.logical_not(adat.mask).astype(n.int)
    for bl in xsum: xsum[bl] /= n.where(cnt[bl] == 0, 1, cnt[bl])
    for c, jd in enumerate(times):
        repfile = '%f.xtalk.rep' % jd
        for bl in xsum: xtalk[bl] = list(gain[bl][c] * xsum[bl])
        print 'Writing', repfile
        f = open(repfile, 'w')
        pickle.dump(xtalk, f)
        f.close()
    import sys; sys.exit(0)

guess, cnt, xtalk = {}, {}, {}
for filename in args:
    print filename,'->',filename+'x'
    if not opts.outfile and os.path.exists(filename+'x'):
        print filename+'x', 'exists.  Skipping...'
        continue
    uv = a.miriad.UV(filename)
    uv.select('auto',0, 0, include=False)
    (uvw,jd,(i,j)),d,f = uv.read(raw=True)
    uv.rewind()
    if opts.infile:
        xfile = '%f.xtalk.rep' % jd
        if not os.path.exists(xfile): xfile = '%f.xtalk' % jd
        if not os.path.exists(xfile):
            print xfile, 'does not exist.  Skipping...'
            continue
        print '    using', xfile
        f = open(xfile); xtalk = pickle.load(f); f.close()
        for bl in xtalk: xtalk[bl] = n.array(xtalk[bl])
    else:
        guess, cnt, xtalk = {}, {}, {}
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            if not guess.has_key(bl): guess[bl], cnt[bl] = 0, 0
            guess[bl] += n.where(f, 0, d)
            cnt[bl] += n.logical_not(f)
        del(uv)
        for bl in guess: xtalk[bl] = guess[bl] / n.clip(cnt[bl], 1, n.Inf)
    if opts.outfile:
        for bl in xtalk: xtalk[bl] = list(xtalk[bl])
        print 'Writing %f.xtalk' % jd
        f = open('%f.xtalk' % jd, 'w')
        pickle.dump(xtalk, f)
        f.close()
    else:
        def mfunc(uv, p, d, f):
            uvw,t,(i,j) = p
            bl = a.miriad.ij2bl(i,j)
            if xtalk.has_key(bl): return p, d - xtalk[bl], f
            else: return p, d, f
        uvi = a.miriad.UV(filename)
        uvo = a.miriad.UV(filename+'x', status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=mfunc,
            append2hist='XTALK3: removed crosstalk\n', raw=True)
