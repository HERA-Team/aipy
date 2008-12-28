import ant, sim, fit

specials = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter',
    'Saturn', 'Uranus', 'Neptune']
src_data = {
    #             RA          DEC         FLUX  FREQ, INDEX
    #'Sun'  : (None        , None,        18800, .150,  0.10 ), #  y   y
    #'Moon'  : (None        , None,           0, .150,  0.10 ), #  y   y
    #'Mercury': (None        , None,          0, .150,  0.10 ), #  y   y
    #'Venus':   (None        , None,          0, .150,  0.10 ), #  y   y
    #'Mars':    (None        , None,          0, .150,  0.10 ), #  y   y
    #'Jupiter': (None        , None,          0, .150,  0.10 ), #  y   y
    #'Saturn' : (None        , None,          0, .150,  0.10 ), #  y   y
    #'Uranus' : (None        , None,          0, .150,  0.10 ), #  y   y
    #'Neptune': (None        , None,          0, .150,  0.10 ), #  y   y
    'b0320': ('03:20:00.0', '-37:18:00',   259, .408, -0.8  ), #  y
    '3c123': ('04:33:55.4', '+29:35:15',   247, .160, -0.73 ), #  n
    'pic'  : ('05:18:20.2', '-45:49:31',   452, .160, -0.8  ), #  y
    #'crab' : ('05:34:32.0', '+22:00:52',  1500, .159, -0.8  ), #  y   y
    'hyd'  : ('09:15:42.4', '-11:53:13',   243, .160, -1.12 ), #  ?
    '3c273': ('12:26:32.1', '+02:19:14',   102, .160, -0.61 ), #  y?
    'vir'  : ('12:30:49.4', '+12:23:28',   566, .160, -1.54 ), #  y   y
    'cen'  : ('13:22:34.3', '-42:44:15',  1104, .160, -0.52 ), #  y
    'her'  : ('16:48:39.3', '+05:04:17',   378, .160, -1.15 ), #  y
    '3c353': ('17:17:54.1', '-00:55:40',   276, .160, -0.49 ), #  ?
    # Fluxes from Baars 1972:
    'cas'  : ('23:23:25.4', '+58:48:38', 12800, .152, -0.787),
    'cyg'  : ('19:59:28.3', '+40:44:02', 10500, .152, -0.8  ), # Index est graph
    'crab' : ('05:34:32.0', '+22:00:52',  1430, .152, -0.263  ),
    #'vir'  : ('12:30:49.4', '+12:23:28',  605, .408, -.853),  # Flux est graph
    # Fit values:
    'Sun'  : (None        , None,        18800, .150,  0.71 ), # Cal w/ cyg
}

def get_src(s, type='fit'):
    """Return a source created out of the parameters in the dictionary srcs.
    'type' can be 'fit', 'sim', 'ant' depending on which aipy module you
    want to use (and what functionality you need)."""
    ra, dec, st, mfreq, index = src_data[s]
    mdl = eval(type)
    if s in specials:
        src = mdl.RadioSpecial(s, st, meas_freq=mfreq, spec_index=index)
    else:
        src = mdl.RadioFixedBody(ra, dec, st, meas_freq=mfreq,
            spec_index=index, name=s)
    return src

def get_catalog(srcs=None, type='fit'):
    if srcs is None: srcs = src_data.keys()
    srcs = [get_src(s, type=type) for s in srcs]
    mdl = eval(type)
    if mdl == sim: mdl = fit
    return mdl.SrcCatalog(srcs)
