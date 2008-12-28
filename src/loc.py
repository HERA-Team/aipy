"""Load an AntennaArray using the calibration information contained in the
specified location-specific module."""

import numpy as n, ant, sim, fit, os, sys, coord

def get_freqs(sdf, sfreq, nchan):
    return n.arange(nchan, dtype=n.float) * sdf + sfreq

def get_aa(loc_key, sdf, sfreq, nchan):
    '''Return the AntennaArray specified by loc_key, which should be the
    name of a module somewhere in Python's path that implements a 
    "get_aa(freqs)" function.  That function should return an AntennaArray 
    initialized with appropriate calibration parameters.  This function simply 
    attempts to import that function, and then performs a simple conversion 
    between sdf (channel width in GHz), sfreq (starting frequency), nchan (#
    of channels) --> freqs (an array of channel centers).'''
    freqs = get_freqs(sdf, sfreq, nchan)
    sys.path.append(os.getcwd())
    exec('from %s import get_aa as _get_aa' % loc_key)
    return _get_aa(freqs)

def get_src_prms(loc_key):
    sys.path.append(os.getcwd())
    try: exec('from %s import src_prms' % loc_key)
    except(ImportError): src_prms = {}
    return src_prms
