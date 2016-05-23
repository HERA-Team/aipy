'''The Mauritius Radio Telescope catalog.  Can be downloaded from
http://www.rri.res.in/surveys/MRT/Catalogues/CatalogueMRTa.txt Copy this file
to "CatalogueMRTa.txt" in the _src directory of your AIPY installation.'''

import aipy as a, numpy as np, os

class MRTCatalog(a.fit.SrcCatalog):
    def fromfile(self, filename):
        f = open(filename)
        addsrcs = []
        for L in f.readlines()[4:]:
            s, ra, dec, jys = L.split()[:4]
            addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=s,
                jys=float(jys), index=0, mfreq=.150))
        self.add_srcs(addsrcs)

MRTFILE = os.path.dirname(__file__) + os.sep + 'CatalogueMRTa.txt'
_mrtcat = None

def get_srcs(srcs=None, cutoff=None):
    # Mechanism for delaying instantiation of catalog until it is accessed
    global _mrtcat
    if _mrtcat is None:
        _mrtcat = MRTCatalog()
        _mrtcat.fromfile(MRTFILE)
    if srcs is None:
        if cutoff is None: srcs = _mrtcat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _mrtcat.keys(): _mrtcat[s].update_jys(fq)
            srcs = [s for s in _mrtcat.keys() if _mrtcat[s].jys[0] > cut]
    srclist = []
    for s in srcs:
        try: srclist.append(_mrtcat[s])
        except(KeyError): pass
    return srclist

