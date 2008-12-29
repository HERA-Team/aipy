"""
Module containing utilities (like parsing of certain command-line arguments) 
for writing scripts.
"""

import miriad, ant, sim, fit, src, numpy as n

def add_standard_options(optparser, ant=False, pol=False, chan=False, 
        loc=False, loc2=False, src=False, dec=False):
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
    elif loc2: optparser.add_option('-l', '--loc', dest='loc', action='append',
        help='Use specified <loc>.py for location-specific calibration.  Can intersperse locations with UV files to use different calibrations for different files.')
    if src: optparser.add_option('-s', '--src', dest='src',
        help='Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".')
    if dec:
        optparser.add_option('-x', '--decimate', dest='decimate', 
            default=1, type='int',
            help='Use only every Nth integration.  Default is 1.')
        optparser.add_option('--dphs', dest='decphs', 
            default=0, type='int',
            help='Offset to use when decimating (i.e. start counting integrations at this number for the purpose of decimation).  Default is 0.')

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
        antopt1 = [ao[0] for ao in antopt if len(ao) == 1]
        antopt2 = [ao for ao in antopt if len(ao) == 2]
        if len(antopt1) == 1: uv.select('antennae', antopt1[0], -1)
        else:
            for a1 in antopt1:
                for a2 in antopt1:
                    if a1 != a2: uv.select('antennae', a1, a2)
        for opt in antopt2:
            a1,a2 = opt
            if a1 < 0: a1,include = -a1,False
            else: include = True
            uv.select('antennae', a1, a2, include=include)
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

def parse_srcs(src_str):
    """Return (src_list,flux_cutoff) based on string argument for src.
    Can be "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>"."""
    if src_str.startswith('all'): return None, None
    try:
        cutoff = float(src_str)
        return None, cutoff
    except(ValueError): pass
    src_opt = src_str.split(',')
    if len(src_opt) == 1:
        src_opt = src_opt[0].split('_')
        if len(src_opt) == 1: return src_opt, None
        ra,dec = src_opt
        s = fit.RadioFixedBody(ra, dec, name=src_str)
        return [s], None
    else:
        return src_opt, None

def files_to_locs(locations, uvfiles, sysargv):
    import optparse
    o = optparse.OptionParser(); add_standard_options(o, loc2=True)
    locs = {}
    pos = [sysargv.index(uv) for uv in uvfiles]
    curloc = locations[0]
    for i in range(len(pos)):
        locs[curloc] = locs.get(curloc, []) + [sysargv[pos[i]]]
        if i < len(pos)-1 and pos[i+1] - pos[i] > 1:
            opts,args = o.parse_args(sysargv[pos[i]+1:pos[i+1]])
            if not opts.loc is None:
                if curloc is None:
                    locs[opts.loc[-1]] = locs[None]
                    del(locs[None])
                curloc = opts.loc[-1]
    return locs

