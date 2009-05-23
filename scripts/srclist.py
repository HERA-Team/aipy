#! /usr/bin/env python
"""
A script for listing sources (and optionally source parameters) from
catalogs.
"""
import aipy as a, ephem, sys, optparse

o = optparse.OptionParser()
o.set_usage('srclist.py [options]')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)
o.add_option('-P','--prms', dest='prms', 
    help='A comma-delimited list of parameters to print for each source.  Can be "*" to list all parameters.  If no parameters are specified, a list of the selected sources is printed.')
o.add_option('-x','--exclude', dest='exclude', 
    help='A string that is parsed similarly to "srcs", but excludes sources that may otherwise have been selected by the "srcs" parameter.')
o.add_option('-j','--juldate', dest='juldate', type='float',
    help='A Julian Date to use for placing sources.')
o.add_option('--ra',dest='ra_rng', 
    help='A range RA1_RA2 of right-ascensions to select for.  Default: None.')
o.add_option('--dec',dest='dec_rng',
    help='A range DEC1_DEC2 of declinations to select for.  Default: None.')
o.add_option('--divstr', dest='divstr', default=' ',
    help='Divider string to use between source names when printing.  Default is " ".')
opts,args = o.parse_args(sys.argv[1:])

srclist,cutoff = a.scripting.parse_srcs(opts.src)

if opts.cal != None:
    cat = a.cal.get_catalog(opts.cal, srcs=srclist, cutoff=cutoff)
else:
    cat = a.src.get_catalog(srcs=srclist, cutoff=cutoff)

if opts.exclude != None:
    xlist,xoff = a.scripting.parse_srcs(opts.exclude)
    if opts.cal != None:
        xcat = a.cal.get_catalog(opts.cal, srcs=xlist, cutoff=xoff)
    else:
        xcat = a.src.get_catalog(srcs=xlist, cutoff=xoff)
else: xcat = {}

o = ephem.Observer()
if opts.juldate is None: o.date = ephem.J2000
else: o.date = a.phs.juldate2ephem(opts.juldate)
o.epoch = o.date

for s in cat.keys():
    try: a.phs.RadioFixedBody.compute(cat[s], o)
    except(TypeError):
        if opts.juldate is None: del(cat[s])
        else: a.phs.RadioSpecial.compute(cat[s],o)

srcs = cat.keys()
srcs = [s for s in srcs if s not in xcat]

if opts.ra_rng != None:
    ra1,ra2 = map(ephem.hours, opts.ra_rng.split('_'))
    if ra1 < ra2:
        srcs = [s for s in srcs if (cat[s]._ra > ra1 and cat[s]._ra < ra2)]
    else:
        srcs = [s for s in srcs if (cat[s]._ra > ra1 or cat[s]._ra < ra2)]

if opts.dec_rng != None:
    dec1,dec2 = map(ephem.degrees, opts.dec_rng.split('_'))
    if dec1 < dec2:
        srcs = [s for s in srcs if (cat[s]._dec > dec1 and cat[s]._dec < dec2)]
    else:
        srcs = [s for s in srcs if (cat[s]._dec > dec1 or cat[s]._dec < dec2)]

srcs.sort()
if opts.prms == None:
    print opts.divstr.join(srcs)
else:
    prms = opts.prms.split(',')
    for s in srcs:
        p = cat.get_params({s:prms})
        a.fit.print_params(p)
