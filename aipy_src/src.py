# Python3 compatibility
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info < (3,):
    FileNotFoundError = IOError

"""
This module provides a front-end interface for accessing sources in 
all catalogs in the _src module of AIPY.
"""

from . import fit, _src
import warnings

def get_catalog(srcs=None, cutoff=None, catalogs=['helm','misc']):
    """Return a source catalog created out of the sources listed in 'srcs',
    or with fluxes above the specified (jy_cutoff,freq_ghz) in 'cutoff'.
    Searches catalogs listed in 'catalogs'."""
    srclist = []
    for c in catalogs:
        try:
            c = getattr(_src, c)
        except(AttributeError):
            continue
        try:
            srclist += c.get_srcs(srcs=srcs, cutoff=cutoff)
        except(FileNotFoundError) as e:
            warnings.warn("Catalog - %s" % str(e), RuntimeWarning)
            continue
    # Add in sources that are already made
    if srcs != None: srclist += [s for s in srcs if type(s) != str]
    return fit.SrcCatalog(srclist)
