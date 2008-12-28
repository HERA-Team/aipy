#! /usr/bin/env python
"""Tarball and compress (using bz2) Miriad UV files.
Author: Aaron Parsons
Date: 8/14/07"""

import sys, os
from optparse import OptionParser

p = OptionParser()
p.set_usage('compress_uv.py [options] *.uv')
p.set_description(__doc__)
p.add_option('-d', '--delete', dest='delete', action='store_true',
    help='Delete a uv file after compressing it')
p.add_option('-x', '--expand', dest='expand', action='store_true',
    help='Inflate tar.bz2 files')

opts, args = p.parse_args(sys.argv[1:])

for i in args:
    print i
    if opts.expand:
        rv = os.system('tar xjf %s' % i)
        if rv != 0: break
        continue
    cmp_name = i + '.tar.bz2'
    if os.path.exists(cmp_name):
        print cmp_name, 'exists; skipping...'
        continue
    rv = os.system('tar jcf %s %s' % (cmp_name, i))
    # If tar fails, delete malformed tarball and then exit
    if rv != 0:
        os.system('rm -rf %s' % (cmp_name))
        break
    if opts.delete: os.system('rm -rf %s' % (i))
