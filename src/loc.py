"""Load calibration information contained in specified location-specific 
modules."""

import numpy as n, src, os, sys

def get_freqs(sdf, sfreq, nchan):
    return n.arange(nchan, dtype=n.float) * sdf + sfreq

def get_aa(*args):
    '''Return the AntennaArray specified by loc_key, which should be the
    name of a module somewhere in Python's path that implements a 
    "get_aa(freqs)" function.  That function should return an AntennaArray 
    initialized with appropriate calibration parameters.  This function simply 
    attempts to import that function, and then performs a simple conversion 
    between sdf (channel width in GHz), sfreq (starting frequency), nchan (#
    of channels) --> freqs (an array of channel centers).'''
    if len(args) == 2: loc_key,freqs = args
    else:
        loc_key,sdf,sfreq,nchan = args
        freqs = get_freqs(sdf, sfreq, nchan)
    sys.path.append(os.getcwd())
    exec('from %s import get_aa as _get_aa' % loc_key)
    return _get_aa(freqs)

def get_catalog(loc_key=None, srcs=None, cutoff=None):
    '''Return the source catalog specified by loc_key, which should be the
    name of a module somewhere in Python's path that implements a
    "get_catalog(srcs, cutoff)" function.  That function should return a
    SrcCatalog initialized with appropriate calibration parameters.  This
    function simply attempts to import that function and then pass it the
    specified srcs, cutoff.  If no such function is found, src.get_catalog()
    is called instead.'''
    sys.path.append(os.getcwd())
    try: exec('from %s import get_catalog as _get_catalog' % loc_key)
    except(ImportError): _get_catalog = src.get_catalog
    return _get_catalog(srcs=srcs, cutoff=cutoff)
