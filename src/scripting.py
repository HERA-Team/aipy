"""
Module containing utilities (like parsing of certain command-line arguments) 
for writing scripts.

Author: Aaron Parsons
"""

import miriad, ant, sim, fit, src, numpy as n

def add_standard_options(optparser, ant=False, pol=False, chan=False, 
        loc=False, src=False, dec=False):
    """Add standard command-line options to an optparse.OptionParser() on an 
    opt in basis (i.e. specify =True for each option to be added)."""
    if ant: optparser.add_option ('-a', '--ant', dest='ant', default='cross',
         help='Select antennas/baselines to include.  Options are "all", "auto", "cross", "<ant1 #>_<ant2 #>,..." (a list of baselines), or "<ant1 #>,..." (a list of antennas).  Default is "cross".')
    if pol: optparser.add_option('-p', '--pol', dest='pol', 
        help='Choose polarization (xx, yy, xy, yx) to include.')
    if chan: optparser.add_option('-c', '--chan', dest='chan', default='all',
        help='Select channels (taken after any delay/fringe transforms) to include.  Options are "all", "<chan1 #>,..." (a list of channels), or "<chan1 #>_<chan2 #>" (a range of channels).  Default is "all".')
    if loc: optparser.add_option('-l', '--loc', dest='loc',
        help='Use specified <loc>.py for location-specific calibration.')
    if src: optparser.add_option('-s', '--src', dest='src',
        help='Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".')
    if dec: optparser.add_option('-x', '--decimate', dest='decimate', 
        default=1, type='int',
        help='Use only every Nth integration.  Default is 1.')

def _strarg_to_range(strarg):
    """Split command-line lists/ranges into a list of numbers."""
    strarg = strarg.split(',')
    return [map(float, option.split('_')) for option in strarg]

def uv_selector(uv, ant_str, pol_str):
    """Call uv.select with appropriate options based on string argument for
    antennas (can be 'all', 'auto', 'cross', '0,1,2', or '0_1,0_2') and
    string for polarization ('xx','yy','xy','yx')."""
    if ant_str.startswith('all'): pass
    elif ant_str.startswith('auto'): uv.select('auto', 0, 0)
    elif ant_str.startswith('cross'): uv.select('auto', 0, 0, include=0)
    else:
        antopt = _strarg_to_range(ant_str)
        for opt in antopt:
            try: a1,a2 = opt
            except(ValueError): a1,a2 = opt + [-1]
            uv.select('antennae', a1, a2)
    try: polopt = miriad.str2pol[pol_str]
    except(KeyError): raise ValueError('--pol argument invalid or absent')
    uv.select('polarization', polopt, 0)

def parse_chans(chan_str, nchan):
    """Return array of active channels based on number of channels and
    string argument for chans (can be 'all', '20_30', or '55,56,57')."""
    if chan_str.startswith('all'): chans = n.arange(nchan)
    else:
        chanopt = _strarg_to_range(chan_str)
        if len(chanopt[0]) != 1:
            chanopt = [n.arange(x,y, dtype=n.int) for x,y in chanopt]
        chans = n.concatenate(chanopt)
    return chans.astype(n.int)

def parse_srcs(src_str, force_cat=False):
    """Return src/catalog based on string argument for src (can be "all",
    "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>").  Can force
    a single source to be returned as catalog with force_cat."""
    if src_str.startswith('all'): return src.get_catalog()
    src_opt = src_str.split(',')
    if len(src_opt) == 1:
        src_opt = src_opt[0].split('_')
        if len(src_opt) == 1: s = src.get_src(src_opt[0])
        else:
            ra,dec = src_opt
            s = fit.RadioFixedBody(ra, dec, 0, name='src')
        if force_cat: return src.get_catalog([s])
        else: return s
    else:
        return src.get_catalog(src_opt)

