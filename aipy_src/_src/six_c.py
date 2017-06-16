'''The 6C (Sixth Cambridge) Catalog.
Data files are in tab-separated format from Vizier.
To download in the correct format, open a catalog online in Vizier,
select'Tab-Separated-Values' as the Output layout in the drop-down box, set
the maximum entries to 'unlimited', and click 'Sexagesimal' under the box
for 'Target Name or Position'.  Submit the query, and copy the output to a
txt file.  Copy files to "6c1.txt", "6c2.txt", "6c3.txt", "6c4.txt", "6c5_1.txt",
and "6c5_2.txt" (for the 5 fields) in the _src directory of your AIPY
installation.'''
import aipy as a, numpy as np, os

class SixCCatalog(a.fit.SrcCatalog):
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
            name = '%s_%s' % (ra,dec)
            jys = float(text[4])
            addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                jys=jys, index=0, mfreq=.151))
        self.add_srcs(addsrcs)

SIXCFILES = [os.path.dirname(__file__) + os.sep + cat 
    for cat in ['6c1.txt','6c2.txt','6c3.txt','6c4.txt','6c5_1.txt','6c5_2.txt']]
_sixccat = None

def get_srcs(srcs=None, cutoff=None):
    global _sixccat
    if _sixccat is None:
        _sixccat = SixCCatalog()
        for f in SIXCFILES: _sixccat.fromfile(f)
    if srcs is None:
        if cutoff is None: srcs = _sixccat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _sixccat.keys(): _sixccat[s].update_jys(fq)
            srcs = [s for s in _sixccat.keys() if _sixccat[s].jys[0] > cut]
    srclist = []
    for s in srcs:
        try: srclist.append(_sixccat[s])
        except(KeyError): pass
    return srclist
