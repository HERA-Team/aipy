'''The NVSS Catalog.
Data files are in tab-separated format from Vizier.
To download in the correct format, open a catalog online in Vizier,
select'Tab-Separated-Values' as the Output layout in the drop-down box, set
the maximum entries to 'unlimited', and click 'Sexagesimal' under the box
for 'Target Name or Position'.  Submit the query, and copy the output to a
txt file.  Copy this file to "nvss.txt" in the _src directory of your AIPY
installation.'''

import aipy as a, numpy as np, os
    
class NVSSCatalog(a.fit.SrcCatalog):
    def fromfile(self, filename):
        f = open(filename)
        addsrcs = []
        for L in [L for L in f.readlines() if not L.startswith('#')]:
            text = L.split('\t')
            if len(text) <= 4: continue
            try: int(text[0][0])
            except(ValueError): continue
            ra = text[3].replace(' ',':')
            dec = text[4].replace(' ',':')
            name = text[2].strip()
            jys = float(text[7])/1000.
            addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                jys=jys, index=0, mfreq=1.40))
        self.add_srcs(addsrcs)

NVSSFILE = os.path.dirname(__file__) + os.sep + 'nvss.txt'
_nvsscat = None

def get_srcs(srcs=None, cutoff=None):
    global _nvsscat
    if _nvsscat is None:
        _nvsscat = NVSSCatalog()
        _nvsscat.fromfile(NVSSFILE)
    if srcs is None:
        if cutoff is None: srcs = _nvsscat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _nvsscat.keys(): _nvsscat[s].update_jys(fq)
            srcs = [s for s in _nvsscat.keys() if _nvsscat[s].jys[0] > cut]

    srclist = []
    for s in srcs:
        try: srclist.append(_nvsscat[s])
        except(KeyError): pass
    return srclist
