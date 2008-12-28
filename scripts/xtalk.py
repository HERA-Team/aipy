#! /usr/bin/env python
"""
A script for removing stable crosstalk from UV files.

Author: Aaron Parsons
Date: 6/03/07
Revisions: None
"""

import aipy.ants, aipy.miriad, aipy.constants, numpy

if __name__ == '__main__':
    import os, sys
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('xtalk.py [options] *.uv')
    p.set_description(__doc__)

    opts, args = p.parse_args(sys.argv[1:])

    for uvfile in args:
        phs_off = {}
        cnt = {}
        print 'Reading', uvfile
        uvi = aipy.miriad.UV(uvfile)

        # Gather all data
        while True:
            p, d = uvi.read_data()
            if d.size == 0: break
            bl = int(p[-1])
            i, j = aipy.miriad.bl2ij(bl)
            if i == j: continue
            try:
                phs_off[bl] += d.filled(0)
                cnt[bl] += numpy.where(d.mask, 0, 1)
            except(KeyError):
                phs_off[bl] = d.filled(0)
                cnt[bl] = numpy.where(d.mask, 0, 1)

        # Average data by channel over entire file
        for bl in phs_off: phs_off[bl] /= numpy.where(cnt[bl] == 0, 1, cnt[bl])
        del(uvi)

        # Generate a pipe for removing average phase bias from data
        def phs_corr_mfunc(uv, preamble, data):
            i, j = aipy.miriad.bl2ij(preamble[-1])
            if i == j: return preamble, data
            return preamble, data - phs_off[preamble[-1]]

        # Apply the pipe to the data
        print 'Working on', uvfile
        uvofile = uvfile+'x'
        if os.path.exists(uvofile):
            print uvofile, 'exists, skipping.'
            continue
        uvi = aipy.miriad.UV(uvfile)
        uvo = aipy.miriad.UV(uvofile, status='new')
        aipy.miriad.pipe_uv(uvi, uvo, mfunc=phs_corr_mfunc)
        del(uvo)
        del(uvi)
