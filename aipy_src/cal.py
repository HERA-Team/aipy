"""Load calibration information contained in specified calibration
modules."""

import numpy as np
import os
import sys
from . import src

def get_freqs(sdf, sfreq, nchan):
    return np.arange(nchan, dtype=np.float) * sdf + sfreq

def get_aa(*args):
    '''Return the AntennaArray specified by cal_key, which should be the
    name of a module somewhere in Python's path that implements a 
    "get_aa(freqs)" function.  That function should return an AntennaArray 
    initialized with appropriate calibration parameters.  This function simply 
    attempts to import that function, and then performs a simple conversion 
    between sdf (channel width in GHz), sfreq (starting frequency), nchan (#
    of channels) --> freqs (an array of channel centers).'''
    if len(args) == 2: cal_key,freqs = args
    else:
        cal_key,sdf,sfreq,nchan = args
        freqs = get_freqs(sdf, sfreq, nchan)
    sys.path.append(os.getcwd())
    exec('from %s import get_aa as _get_aa' % cal_key)
    return _get_aa(freqs)

def get_catalog(cal_key=None, srcs=None, cutoff=None, catalogs=['helm','misc']):
    '''Return the source catalog specified by cal_key, which should be the
    name of a module somewhere in Python's path that implements a
    "get_catalog(srcs, cutoff, catalogs)" function.  That function should return a
    SrcCatalog initialized with appropriate calibration parameters.  This
    function simply attempts to import that function and then pass it the
    specified srcs, cutoff.  If no such function is found, src.get_catalog()
    is called instead.'''
    sys.path.append(os.getcwd())
    try: exec('from %s import get_catalog as _get_catalog' % cal_key)
    except(ImportError): _get_catalog = src.get_catalog
    return _get_catalog(srcs=srcs, cutoff=cutoff, catalogs=catalogs)
