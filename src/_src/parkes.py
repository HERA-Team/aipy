'''The Parkes Catalog.  
Data files are in tab-separated format from Vizier.
To download in the correct format, open a catalog online in Vizier,
select'Tab-Separated-Values' as the Output layout in the drop-down box, set
the maximum entries to 'unlimited', and click 'Sexagesimal' under the box
for 'Target Name or Position'.  Submit the query, and copy the output to a
txt file.  Copy this file to "parkes.txt" in the _src directory of your AIPY
installation.'''
import aipy as a, numpy as np, os

class ParkesCatalog(a.fit.SrcCatalog):
    def fromfile(self,filename):
        f = open(filename)
        addsrcs = []
        for L in [L for L in f.readlines() if not L.startswith('#')]:
            text = L.split('\t')
            if len(text) <= 4: continue
            try: int(text[0][0])
            except(ValueError): continue
            ra = text[0].replace(' ',':')
            dec = text[1].replace(' ',':')
            name = text[11].strip()
            try: jys = float(text[9])
            except(ValueError): continue
            addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                jys=jys, index=0, mfreq=2.70))
        self.add_srcs(addsrcs)

PARKESFILE = os.path.dirname(__file__) + os.sep + 'parkes.txt'
_parkescat = None

def get_srcs(srcs=None, cutoff=None):
    global _parkescat
    if _parkescat is None:
        _parkescat = ParkesCatalog()
        _parkescat.fromfile(PARKESFILE)
    if srcs is None:
        if cutoff is None: srcs = _parkescat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _parkescat.keys(): _parkescat[s].update_jys(fq)
            srcs = [s for s in _parkescat.keys() if _parkescat[s].jys[0] > cut]
    srclist = []
    for s in srcs:
        try: srclist.append(_parkescat[s])
        except(KeyError): pass
    return srclist
