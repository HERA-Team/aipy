'''The 7C (Seventh Cambridge) Catalog.
Data files are in tab-separated format from Vizier.
To download in the correct format, open a catalog online in Vizier,
select'Tab-Separated-Values' as the Output layout in the drop-down box, set
the maximum entries to 'unlimited', and click 'Sexagesimal' under the box
for 'Target Name or Position'.  Submit the query, and copy the output to a
txt file.  Copy this file to "7c.txt" in the _src directory of your AIPY
installation.'''
import aipy as a, numpy as np, os

class SevenCCatalog(a.fit.SrcCatalog):
    def fromfile(self, filename):
        f = open(filename)
        addsrcs = []
        for L in [L.strip('\n') for L in f.readlines() if not L.startswith('#')]:
            text = L.split('\t')
            if len(text) <= 4: continue
            try: int(text[0][0])
            except(ValueError): continue
            ra = text[0].replace(' ',':')
            dec = text[1].replace(' ',':')
            name = '%s_%s' % (ra, dec)
            jys = float(text[5])
            addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                jys=jys, index=0, mfreq=0.151))
        self.add_srcs(addsrcs)

SEVENCFILE = os.path.dirname(__file__) + os.sep + '7c.txt'
_sevenccat = None

def get_srcs(srcs=None, cutoff=None):
    global _sevenccat
    if _sevenccat is None:
        _sevenccat = SevenCCatalog()
        _sevenccat.fromfile(SEVENCFILE)
    if srcs is None:
        if cutoff is None: srcs = _sevenccat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _sevenccat.keys(): _sevenccat[s].update_jys(fq)
            srcs = [s for s in _sevenccat.keys() if _sevenccat[s].jys[0] > cut]
    srclist = []
    for s in srcs:
        try: srclist.append(_sevenccat[s])
        except(KeyError): pass
    return srclist

