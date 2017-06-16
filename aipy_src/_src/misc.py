import aipy as a, numpy as np

specials = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter',
    'Saturn', 'Uranus', 'Neptune']
src_data = {
    #             RA          DEC         FLUX  FREQ, INDEX ANGSIZE
    'Moon'  : (None        , None,           0, .150,  0., 0.),
    'Mercury': (None        , None,          0, .150,  0., 0.),
    'Venus':   (None        , None,          0, .150,  0., 0.),
    'Mars':    (None        , None,          0, .150,  0., 0.),
    'Jupiter': (None        , None,          0, .150,  0., 0.),
    'Saturn' : (None        , None,          0, .150,  0., 0.),
    'Uranus' : (None        , None,          0, .150,  0., 0.),
    'Neptune': (None        , None,          0, .150,  0., 0.),
    'Sun'  : (None, None, 57000, .150,  2.00, 4.6e-3),
    'for':   ('03:22:41.7', '-37:12:30',   170, .150, -0.8 , 0.),
    'pic'  : ('05:19:49.7', '-45:46:45',   452, .160, -0.8 , 0.),
    'hyd'  : ('09:18:05.7', '-12:05:44',  1860, .160, -2.30, 0.),
    'cen'  : ('13:25:27.6', '-43:01:09', 5144, .150, -0.65, 0.), # From Israel 1998
    'her' : ('16:51:08.15', '4:59:33.3', 300.0, 0.159, -1, 0.000669),
    'sgr'  : ('17:45:40.0', '-29:00:28',  121, .160, -4.21, 0.),
    # Fluxes from Miriad:
    'crab' : ('05:34:32.0', '+22:00:52',  1838, .150, -0.30, 0.),
    'vir'  : ('12:30:49.4', '+12:23:28',   1446, .150, -0.86, 0.),
    'cyg'  : ('19:59:28.3', '+40:44:02', 10900, .150, -0.69, (0.,0.,0.)),
    'cas'  : ('23:23:27.94', '+58:48:42.4',  9160.0, 0.150, -0.73, 0.000000),
}

_misccat = None

def get_srcs(srcs=None, cutoff=None):
    global _misccat
    if _misccat is None:
        _misccat = a.fit.SrcCatalog()
        srclist = []
        for s in src_data:
            ra, dec, st, mfreq, index, srcshape = src_data[s]
            try: len(srcshape)
            except(TypeError): srcshape = (srcshape, srcshape, 0.)
            if s in specials:
                srclist.append(a.fit.RadioSpecial(s, st, mfreq=mfreq, 
                    index=index, srcshape=srcshape))
            else:
                srclist.append(a.fit.RadioFixedBody(ra, dec, jys=st, mfreq=mfreq, 
                    index=index, name=s, srcshape=srcshape))
        _misccat.add_srcs(srclist)
    if srcs is None:
        if cutoff is None: srcs = _misccat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _misccat.keys(): _misccat[s].update_jys(fq)
            srcs = [s for s in _misccat.keys() if _misccat[s].jys[0] > cut]

    srclist = []
    for s in srcs:
        try: srclist.append(_misccat[s])
        except(KeyError): pass
    return srclist

