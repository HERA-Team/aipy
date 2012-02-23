#! /usr/bin/env python
"""
A script for listing sources (and optionally source parameters) from
catalogs.
"""
import aipy as a, ephem, sys, optparse,logging

o = optparse.OptionParser()
o.set_usage('srclist.py [options]')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)
o.add_option('-P','--prms', dest='prms', 
    help='A comma-delimited list of parameters to print for each source.  Can be "*" to list all parameters.  If no parameters are specified, a list of the selected sources is printed.')
o.add_option('-x','--exclude', dest='exclude', 
    help='A string that is parsed similarly to "srcs", but excludes sources that may otherwise have been selected by the "srcs" parameter.')
o.add_option('-c','--cen', dest='centers', 
    help='A string that is parsed similarly to "srcs", but is used along with --sep to select sources near specified locations.')
o.add_option('-j','--juldate', dest='juldate', type='float',
    help='A Julian Date to use for placing sources.')
o.add_option('--ra',dest='ra_rng', 
    help='A range RA1_RA2 of right-ascensions to select for.  Default: None.')
o.add_option('--dec',dest='dec_rng',
    help='A range DEC1_DEC2 of declinations to select for.  Default: None.')
o.add_option('--sep',dest='sep', type='float', 
    help='Include areas within the specified angular separation (in degrees) of any sources listed in --src.')
o.add_option('--divstr', dest='divstr', default=' ',
    help='Divider string to use between source names when printing.  Default is " ".')
o.add_option('-v',dest='verb',action='store_true',
    help="Print more")
o.add_option('--fitprms',action='store_true',
    help="Print as prm inputs for fitmdl")
opts,args = o.parse_args(sys.argv[1:])

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
if opts.verb:
    logging.basicConfig(level=logging.DEBUG)
else:
     logging.basicConfig(level=logging.WARNING)
log = logging.getLogger('srclist')

if opts.cal != None:
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    cat = a.src.get_catalog(srclist, cutoff, catalogs)

if opts.exclude != None:
    xlist,xoff,catalogs = a.scripting.parse_srcs(opts.exclude, opts.cat)
    if opts.cal != None:
        xcat = a.cal.get_catalog(opts.cal, xlist, xoff, catalogs)
    else:
        xcat = a.src.get_catalog(xlist, xoff, catalogs)
else: xcat = {}

if opts.centers != None:
    assert(opts.sep != None)
    clist,coff,catalogs = a.scripting.parse_srcs(opts.centers, opts.cat)
    if opts.cal != None:
        ccat = a.cal.get_catalog(opts.cal, clist, coff, catalogs)
    else:
        ccat = a.src.get_catalog(clist, coff, catalogs)
else: ccat = {}
    
if opts.juldate is None: date = ephem.J2000
else: 
    date = a.phs.juldate2ephem(opts.juldate)
    aa = a.scripting.get_null_aa()
    aa.set_jultime(opts.juldate)

for c in [cat, xcat, ccat]:
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: c[s].compute(aa)

srcs = cat.keys()
srcs = [s for s in srcs if s not in xcat]
if opts.sep != None:
    nsrcs = []
    for s1 in srcs:
        for s2 in ccat.keys():
            if ephem.separation(cat[s1], ccat[s2]) <= opts.sep * a.img.deg2rad:
                nsrcs.append(s1)
                break
    srcs = nsrcs
    

if opts.ra_rng != None:
    ra1,ra2 = map(ephem.hours, opts.ra_rng.split('_'))
    if ra1 < ra2:
        srcs = [s for s in srcs if (cat[s].ra > ra1 and cat[s].ra < ra2)]
    else:
        srcs = [s for s in srcs if (cat[s].ra > ra1 or cat[s].ra < ra2)]

if opts.dec_rng != None:
    dec1,dec2 = map(ephem.degrees, opts.dec_rng.split('_'))
    if dec1 < dec2:
        srcs = [s for s in srcs if (cat[s].dec > dec1 and cat[s].dec < dec2)]
    else:
        srcs = [s for s in srcs if (cat[s].dec > dec1 or cat[s].dec < dec2)]

# We're done selecting sources now.  Time to print information
srcs.sort()
if opts.prms == None:
    print opts.divstr.join(srcs)
else:
    if not opts.fitprms:
        prms = opts.prms.split(',')
        for s in srcs:
            p = cat.get_params({s:prms})
            if p[s].has_key('ra'):  p[s]['ra'] = ephem.hours(p[s]['ra'])
            if p[s].has_key('dec'): p[s]['dec'] = ephem.degrees(p[s]['dec'])
            a.fit.print_params(p)
    else:
        prms = opts.prms.split(',')
        outstring = ''
        scount =0
        snum = len(srcs)
        for s in srcs:
            scount +=1
            outstring += s+"=("
            p = cat.get_params({s:prms})
            pcnt = len(p[s].keys())
            for i,prm in enumerate(p[s].keys()):
                if i<pcnt-1: outstring += prm+"/"
                else: outstring += prm
            if scount<snum-1:outstring += "),"
            else: outstring += ")"
        print outstring
