#! /usr/bin/env python
import aipy as a, ephem as e, numpy as np, sys, optparse, os

o = optparse.OptionParser()
o.set_usage('modmap.py [options]')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)
o.add_option('-m', '--map', dest='map', default='out.fits',
    help='The location to save the map.')
o.add_option('-i', '--in_map', dest='in_map', 
    help='Existing map to add sources to.  Treated as temperature, i.e. if resolution is degraded, pixels are averaged.  Use mscale to re-scale if you want Jy per pixel.')
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Frequency to use for source strengths.  Default .150 GHz.')
o.add_option('-j', '--juldate', dest='juldate', type='float',
    help='Julian date used for locating moving sources.')
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
o.add_option('--nindices', dest='nindices', type='int',default=0,
    help='Number of terms of the spectral index polynomial to save in map.  Default is 0.')
o.add_option('--mscale', dest='mscale', type='float',default=1,
    help='Scaling to apply to data in map.  Useful for changing T into Jy per pixel, and reverse.')
o.add_option('--sscale', dest='sscale', type='float',default=1,
    help='Scaling to apply to sources.  Useful for changing Jy to T.  It may be useful to know that the size of a map pixel is 4*pi / (12*nside**2).')
o.add_option('--iepoch', dest='iepoch', default=e.J2000,
    help='The epoch of coordinates in the input map.  Default J2000.')
o.add_option('--oepoch', dest='oepoch', default=e.J2000, 
    help='The epoch of coordinates in the output map.  Default J2000.')
o.add_option('--isys', dest='isys', default='eq',
    help='Coordinates (eq,ec,ga) of input map.  Default is eq')
o.add_option('--osys', dest='osys', default='eq',
    help='Coordinates (eq,ec,ga) of output map.  Default is eq')
o.add_option('--dtype', dest='dtype', 
    help='Data type (float, double) of output map.  Default is same as input map, or double if no input map is specified.')
opts,args = o.parse_args(sys.argv[1:])

assert(opts.in_map != None or opts.nside != None)
assert(opts.dtype == None or opts.dtype in ['float', 'double'])
if opts.dtype != None:
    if opts.dtype.startswith('fl'): opts.dtype = np.float32
    else: opts.dtype = np.double

if not opts.in_map is None:
    imap = a.map.Map(fromfits=opts.in_map)
    if opts.nside is None and opts.dtype is None:
        print 'Using data from %s' % opts.in_map
        h = imap
    else:
        if opts.dtype is None: opts.dtype = imap.get_dtype()
        if opts.nside is None: opts.nside = imap.nside()
        h = a.map.Map(nside=opts.nside,nindices=opts.nindices,dtype=opts.dtype)
        print 'Importing data from %s' % opts.in_map
        h.from_map(a.map.Map(fromfits=opts.in_map))
    h.map.map *= opts.mscale
    px = np.arange(h.npix())
    m = a.coord.convert_m(opts.osys, opts.isys,
        iepoch=opts.oepoch, oepoch=opts.iepoch)
    x,y,z = np.dot(m,h.px2crd(px, ncrd=3))
    h.set_interpol(True)
    h.map[px] = h.map[x,y,z]
    h.wgt[px] = h.wgt[x,y,z]
    for i in h.ind: i[px] = i[x,y,z]
    h.set_interpol(False)
else:
    if opts.dtype is None: opts.dtype = np.double
    h = a.map.Map(nside=opts.nside, nindices=opts.nindices, dtype=opts.dtype)
    print 'Starting a new map.'
print 'NSIDE:', h.nside()
h.reset_wgt()

if opts.src is None:
    print 'Saving to', opts.map
    h.to_fits(opts.map)
    sys.exit(0)

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
if not opts.cal is None:
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
    aa = a.cal.get_aa(opts.cal, .1, opts.freq, 1)
    aa.set_jultime(opts.juldate)
    cat.compute(aa)
else:
    cat = a.src.get_catalog(srclist, cutoff, catalogs)

m = a.coord.convert_m('eq', opts.osys, oepoch=opts.oepoch)
ths,phis = h.px2crd(np.arange(h.npix()), ncrd=2)
ras,decs = phis, np.pi/2 - ths
afreq = np.array([opts.freq])
for srcname in cat:
    src = cat[srcname]
    eq = e.Equatorial(src._ra, src._dec, epoch=e.J2000)
    eq = e.Equatorial(eq, epoch=opts.oepoch)
    ra,dec = eq.get()
    # Account for size of source
    a1,a2,th = src.srcshape
    dras,ddecs = ras - ra, decs - dec
    print '--------------------------------------------------'
    src.update_jys(afreq)
    strength = src.get_jys()[0] * opts.sscale
    print 'Adding', srcname, 'with strength %f Jy' % strength,
    print 'and index', src.index
    print 'Source shape: a1=%f, a2=%f, th=%f' % (a1, a2, th)
    da1 = dras*np.cos(th) - ddecs*np.sin(th)
    da2 = dras*np.sin(th) + ddecs*np.cos(th)
    delta = (da1/a1)**2 + (da2/a2)**2
    px = np.where(delta <= 1)[0]
    if len(px) == 0:
        print 'Treating as point source.'
        eq = a.coord.radec2eq((ra,dec))
        x,y,z = np.dot(m, eq)
        px = h.crd2px(np.array([x]),np.array([y]),np.array([z]))
    str_per_px = strength / len(px)
    print 'Putting %f in each of %d pixels' % (str_per_px, len(px))
    if opts.nindices == 0:
        for p in px:
            curflux = h[p]
            h.put(p,1,str_per_px+curflux)
    else:
        for p in px:
            curflux, i = h[p]
            h.put(p,1,str_per_px+curflux,[src.index])

print 'Saving to', opts.map
h.to_fits(opts.map)

