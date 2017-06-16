"""
Module containing utilities (like parsing of certain command-line arguments) 
for writing scripts.
"""

import miriad, fit, src, numpy as np, re,phs

def add_standard_options(optparser, ant=False, pol=False, chan=False, 
        cal=False, src=False, prms=False, dec=False, cmap=False, 
        max=False, drng=False):
    """Add standard command-line options to an optparse.OptionParser() on an 
    opt in basis (i.e. specify =True for each option to be added)."""
    if ant: optparser.add_option ('-a', '--ant', dest='ant', default='cross',
         help='Select ant_pol/baselines to include.  Examples: all (all baselines) auto (of active baselines, only i=j) cross (only i!=j) 0,1,2 (any baseline involving listed ants) 0_2,0_3 (only listed baselines) "(0,1)_(2,3)" (same as 0_2,0_3,1_2,2_3. Quotes help bash deal with parentheses) "(-0,1)_(2,-3)" (exclude 0_2,0_3,1_3 include 1_2).  Default is "cross". Select pol by adding appropriate x or y eg 5x_6y.')
    if pol: optparser.add_option('-p', '--pol', dest='pol', default=-1,
        help='Choose polarization (xx, yy, xy, yx) to include.')
    if chan: optparser.add_option('-c', '--chan', dest='chan', default='all',
        help='Select channels (after any delay/delay-rate transforms) to include.  Examples: all (all channels), 0_10 (channels from 0 to 10, including 0 and 10) 0_10_2 (channels from 0 to 10, counting by 2), 0,10,20_30 (mix of individual channels and ranges).  Default is "all".')
    if cal: optparser.add_option('-C', '--cal', dest='cal', 
        help='Use specified <cal>.py for calibration information.')
    if src:
        optparser.add_option('-s', '--src', dest='src',
            help='Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".')
        optparser.add_option('--cat', dest='cat', default='helm,misc',
            help='A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.')
    if prms: optparser.add_option('-P', '--prms', dest='prms',
        help='Parameters (for fitting, usually), can be specified as can be: "obj=prm", "obj=prm/val", "obj=prm/val/sig", "(obj1/obj2)=prm/(val1/val2)/sig", "obj=(prm1/prm2)/val/(sig1/sig2)", comma separated versions of the above, and so on.')
    if dec:
        optparser.add_option('-x', '--decimate', dest='decimate', 
            default=1, type='int',
            help='Use only every Nth integration.  Default is 1.')
        optparser.add_option('--dphs', dest='decphs', 
            default=0, type='int',
            help='Offset to use when decimating (i.e. start counting integrations at this number for the purpose of decimation).  Default is 0.')
    if cmap: optparser.add_option('--cmap', dest='cmap', default='jet',
        help='Colormap for plotting.  Can be gist_earth, gist_heat, gist_stern, gist_yarg, hot, cool, gray, bone, spectral, copper, jet to name a few.  For a more complete list, see pylab.cm.datad.keys().  Default is jet.')
    if max: optparser.add_option( '--max',dest='max',type='float',default=None,
    help='Manually set the maximum color level, in units matching plotting mode.  Default max(data).')
    if drng:
        optparser.add_option('--drng', dest='drng', type='float', default=None,
    help="Dynamic range in color of image, in units matching plotting mode.  Default max(data)-min(data).")

ant_re = r'(\(((-?\d+[xy]?,?)+)\)|-?\d+[xy]?)'
bl_re = '(^(%s_%s|%s),?)' % (ant_re, ant_re, ant_re)
def parse_ants(ant_str, nants):
    """Generate list of (baseline, include, pol) tuples based on parsing of the
    string associated with the 'ants' command-line option."""
    rv,cnt = [], 0
    while cnt < len(ant_str):
        m = re.search(bl_re, ant_str[cnt:])
        if m is None:
            if ant_str[cnt:].startswith('all'): rv = []
            elif ant_str[cnt:].startswith('auto'): rv.append(('auto',1,-1))
            elif ant_str[cnt:].startswith('cross'): rv.append(('auto',0,-1))
            else: raise ValueError('Unparsible ant argument "%s"' % ant_str)
            c = ant_str[cnt:].find(',')
            if c >= 0: cnt += c + 1
            else: cnt = len(ant_str)
        else:
            m = m.groups()
            cnt += len(m[0])
            if m[2] is None:
                ais = [m[8]]
                ajs = range(nants)
            else:
                if m[3] is None: ais = [m[2]]
                else: ais = m[3].split(',')
                if m[6] is None: ajs = [m[5]]
                else: ajs = m[6].split(',')
            include = None
            for i in ais:
                for j in ajs:
                    if type(i) == str and i.startswith('-'):
                         i = i[1:] #nibble the - off the string
                         include = 0 
                    if type(j) == str and j.startswith('-'):
                        j = j[1:]
                        include = 0
                    elif include==0:pass
                    else: include = 1
                    pol = None
                    i,j = str(i),str(j)
                    if not i.isdigit():
                        ai = re.search(r'(\d+)([x,y])',i).groups()
                    if not j.isdigit():
                        aj = re.search(r'(\d+)([x,y])',j).groups()
                    if i.isdigit() and not j.isdigit():
                        pol = ['x'+aj[1],'y'+aj[1]]
                        ai = [i,'']
                    elif not i.isdigit() and j.isdigit():
                        pol = [ai[1]+'x',ai[1]+'y']
                        aj = [j,'']
                    elif not i.isdigit() and not j.isdigit():
                        pol = [ai[1]+aj[1]]
                    if not pol is None:
                        bl = miriad.ij2bl(abs(int(ai[0])),abs(int(aj[0])))
                        for p in pol:
                            rv.append((bl,include,p))
                    else: 
                        bl = miriad.ij2bl(abs(int(i)),abs(int(j)))
                        rv.append((bl,include,-1))
    return rv

def uv_selector(uv, ants=-1, pol_str=-1):
    """Call uv.select with appropriate options based on string argument for
    antennas (can be 'all', 'auto', 'cross', '0,1,2', or '0_1,0_2') and
    string for polarization ('xx','yy','xy','yx')."""
    if ants != -1:
        if type(ants) == str: ants = parse_ants(ants, uv['nants'])
        for cnt,(bl,include,pol) in enumerate(ants):
            if cnt > 0:
                if include: uv.select('or',-1,-1)
                else: uv.select('and',-1,-1)
            if pol == -1: pol = pol_str # default to explicit pol parameter
            if bl == 'auto': uv.select('auto', 0, 0, include=include)
            else:
                i,j = miriad.bl2ij(bl)
                uv.select('antennae', i, j, include=include)
            if pol != -1:
                for p in pol.split(','):
                    polopt = miriad.str2pol[p]
                    uv.select('polarization', polopt, 0)
    elif pol_str != -1:
        for p in pol_str.split(','):
            polopt = miriad.str2pol[p]
            uv.select('polarization', polopt, 0)


def parse_chans(chan_str, nchan, concat=True):
    """Return array of active channels based on number of channels and
    string argument for chans (all, 20_30, or 55,56,57, or 20_30,31,32).
    Channel ranges include endpoints (i.e. 20_30 includes both 20 and 30)."""
    if chan_str.startswith('all'): chanopt = [np.arange(nchan)]
    else:
        chanopt = []
        for co in chan_str.split(','):
            co = map(int, co.split('_'))
            assert(len(co) in [1,2,3])
            if len(co) == 1: chanopt.append(np.array(co))
            elif len(co) == 2: chanopt.append(np.arange(co[0],co[1]+1))
            else: chanopt.append(np.arange(co[0],co[1]+1,co[2]))
    if concat: return np.concatenate(chanopt)
    return chanopt

def parse_srcs(src_str, cat_str):
    """Return (src_list,flux_cutoff,catalogs) based on string argument for src and
    cat.  Can be "all", "<src_name1>,...", "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>", or
    "val/freq" (sources with Jy flux density above val at freq in GHz)."""
    cats = cat_str.split(',')
    if src_str.startswith('all'): return None, None, cats
    if src_str.find('/') != -1:
        cutoff = map(float, src_str.split('/'))
        return None, cutoff, cats
    src_opt = src_str.split(',')
    for i, s in enumerate(src_opt):
        radec = s.split('_')
        try:
            assert(len(radec) == 2)
            ra,dec = radec
            # Check that this is really an ra_dec pair
            for piece in ra.split(':'): piece = float(piece)
            for piece in dec.split(':'): piece = float(piece)
            src_opt[i] = fit.RadioFixedBody(ra,dec, name=s)
        except(AssertionError,ValueError): continue
    return src_opt, None, cats

name = r'([^\(/,\)=]+)'
grp = r'(%s|\((%s(/%s)*)\))' % tuple([name]*3)
prm = r'(%s=%s(/(%s)?(/%s)?)?)' % tuple([grp]*4)
prm_rgx = re.compile(prm)
def parse_prms(prm_str):
    """Return a dict of the form: {'obj': {'prm':(val,sig),...}...} where
    val is a starting value and sig is a known error associated with that 
    start value.  Both default to None if a value is not provided.  The
    string to be parsed can be: "obj=prm", "obj=prm/val", 
    "obj=prm/val/sig", "(obj1/obj2)=prm/(val1/val2)/sig", 
    "obj=(prm1/prm2)/val/(sig1/sig2)", comma separated versions of the above,
    and so on."""
    prms = {}
    for prm in prm_str.split(','):
        m = prm_rgx.match(prm)
        g = m.groups()
        if g[2]: obj = [g[2]]
        else: obj = g[3].split('/')
        if g[8]: plist = [g[8]]
        else: plist = g[9].split('/')
        if g[16]: ival = [float(g[16])]
        elif g[17]: ival = map(float, g[17].split('/'))
        else: ival = [None]
        if g[23]: sval = [float(g[23])]
        elif g[24]: sval = map(float, g[24].split('/'))
        else: sval = [None]
        if len(obj) != 1:
            if len(plist) != 1:
                assert(len(ival) == 1 and len(sval) == 1)
                ival = [ival * len(plist)] * len(obj)
                sval = [sval * len(plist)] * len(obj)
            else:
                if len(ival) == 1: ival = ival * len(obj)
                if len(sval) == 1: sval = sval * len(obj)
                assert(len(ival) == len(obj) and len(sval) == len(obj))
                ival = [[i] for i in ival]
                sval = [[i] for i in sval]
        else:
            if len(plist) != 1:
                if len(ival) == 1: ival = ival * len(plist)
                if len(sval) == 1: sval = sval * len(plist)
                assert(len(ival) == len(plist) and len(sval) == len(plist))
            else:
                assert(len(ival) == 1 and len(sval) == 1)
            ival = [ival]
            sval = [sval]
        for o,il,sl in zip(obj,ival,sval):
            if not prms.has_key(o): prms[o] = {}
            for p,i,s in zip(plist,il,sl):
                prms[o][p] = (i,s)
    return prms

def get_null_aa():
   return phs.AntennaArray([0,0],
       [phs.Antenna(0,0,0,phs.Beam(np.array([0.15])))])
