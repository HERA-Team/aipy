#! /usr/bin/env python
import aipy.miriad, sys, pylab, numpy, math
from optparse import OptionParser

p = OptionParser()
p.set_usage('plot_uv.py [options] *.uv')
p.add_option('-p', '--phase', dest='phase', action='store_true',
    help='Plot phase.')
p.add_option('-l', '--linear', dest='linear', action='store_true',
    help='Plot on linear scale.')
p.add_option('-a', '--autos', dest='autos', action='store_true',
    help='Plot only auto-correlations.')
p.add_option('-c', '--cross', dest='cross', action='store_true',
    help='Plot only cross-correlations.')
p.add_option('-i', '--ant_i', dest='ant_i', default=-1, type='int',
    help='Only plot baselines including this antenna.')
p.add_option('-j', '--ant_j', dest='ant_j', default=-1, type='int',
    help='Only plot baselines which include both ant_i and this antenna.')
p.set_description(__doc__)

opts, args = p.parse_args(sys.argv[1:])

data = {}
for uvfile in args:
    print 'Reading', uvfile
    uv = aipy.miriad.UV(uvfile)
    while True:
        p, d = uv.read_data()
        if d.size == 0: break
        bl = p[-1]
        i, j = aipy.miriad.bl2ij(bl)
        if opts.ant_i != -1:
            if opts.ant_j != -1:
                if opts.ant_i != i or opts.ant_j != j: continue
            if opts.ant_i != i and opts.ant_i != j: continue
        if opts.autos and i != j: continue
        if opts.cross and i == j: continue
        if not data.has_key(bl): data[bl] = []
        d.shape = (1,) + d.shape
        if opts.phase: data[bl].append(numpy.angle(d.filled(0)))
        elif opts.linear: data[bl].append(numpy.abs(d))
        else: data[bl].append(numpy.ma.log10(numpy.abs(d)+1e-6))
    del(uv)

ks = data.keys()
ks.sort()
m2 = int(math.sqrt(len(ks)))
m1 = int(math.ceil(float(len(ks)) / m2))
for n, k in enumerate(ks):
    d = numpy.ma.concatenate(data[k], axis=0)
    pylab.subplot(m1, m2, n+1)
    pylab.imshow(d, aspect='auto')
    pylab.title(str(aipy.miriad.bl2ij(k)))
pylab.show()
