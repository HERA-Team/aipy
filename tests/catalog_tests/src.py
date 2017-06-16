'''This module provides a front-end interface for accessing sources in 
all catalogs in the _src module of AIPY.'''
import fit, _src

def get_catalog(srcs=None, cutoff=None, catalogs=['helm','misc']):
    """Return a source catalog created out of the sources listed in 'srcs',
    or with fluxes above the specified (jy_cutoff,freq_ghz) in 'cutoff'.
    Searches catalogs listed in 'catalogs'."""
    srclist = []
    for c in catalogs:
        try: c = getattr(_src, c)
        except(AttributeError): continue
        srclist += c.get_srcs(srcs=srcs, cutoff=cutoff)
    # Add in sources that are already made
    if srcs != None: srclist += [s for s in srcs if type(s) != str]
    return fit.SrcCatalog(srclist)
