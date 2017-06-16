"""
The PAPER catalog here is derived by experimental methods in either facets or 
healpix maps.
"""

import aipy as a, numpy as np, os
    
class PAPERCatalog(a.fit.SrcCatalog):
    def fromfile(self, filename):
        f = open(filename)
        addsrcs = []
        for L in [L for L in f.readlines() if not L.startswith('#')]:
            text = L.split('\t')
            if len(text) <= 3: continue
            try: int(text[0].strip()[0])
            except(ValueError): continue
            ra = text[1]
            dec = text[2]
            name = text[3].strip()
            jys = float(text[4])
            addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                jys=jys, index=0, mfreq=.150))
        self.add_srcs(addsrcs)

PAPERFILE = os.path.dirname(__file__) + os.sep + 'paper.txt'
_papercat = None

def get_srcs(srcs=None, cutoff=None):
    global _papercat
    if _papercat is None:
        _papercat = PAPERCatalog()
        _papercat.fromfile(PAPERFILE)
    if srcs is None:
        if cutoff is None: srcs = _papercat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _papercat.keys(): _papercat[s].update_jys(fq)
            srcs = [s for s in _papercat.keys() if _papercat[s].jys[0] > cut]

    srclist = []
    for s in srcs:
        try: srclist.append(_papercat[s])
        except(KeyError): pass
    return srclist
