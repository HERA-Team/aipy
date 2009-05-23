#! /usr/bin/env python
"""
A script for listing sources (and optionally source parameters) from
catalogs.
"""
import aipy as a, sys, optparse

o = optparse.OptionParser()
o.set_usage('srclist.py [options]')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)
o.add_option('-P','--prms', dest='prms', 
    help='A comma-delimited list of parameters to print for each source.  Can be "*" to list all parameters.  If no parameters are specified, a list of the selected sources is printed.')
o.add_option('--divstr', dest='divstr', default=' ',
    help='Divider string to use between source names when printing.  Default is " ".')
opts,args = o.parse_args(sys.argv[1:])

srclist,cutoff = a.scripting.parse_srcs(opts.src)
if opts.cal != None:
    cat = a.cal.get_catalog(opts.cal, srcs=srclist, cutoff=cutoff)
else:
    cat = a.src.get_catalog(srcs=srclist, cutoff=cutoff)
srcnames = cat.keys()
srcnames.sort()
if opts.prms == None:
    print opts.divstr.join(srcnames)
else:
    prms = opts.prms.split(',')
    for s in srcnames:
        p = cat.get_params({s:prms})
        a.fit.print_params(p)
