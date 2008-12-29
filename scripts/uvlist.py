#! /usr/bin/env python
"""
Print headers/variables in a Miriad UV file.  If no keywords are provided,
will print a list of available keywords.

Author: Aaron Parsons
"""

import aipy as a, sys, optparse

o = optparse.OptionParser()
o.set_usage('uvlist.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-k', '--key', dest='key',
    help='Comma-delimited list of keywords to print for each UV file.')
opts,args = o.parse_args(sys.argv[1:])

for uvfile in args:
    print uvfile
    uv = a.miriad.UV(uvfile)
    if opts.key is None: print '    ', uv.items() + uv.vars()
    else:
        for key in opts.key.split(','):
            print '    ', key
            print '        ', uv[key]
    print '-----------------------------------------------------------'
