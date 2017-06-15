'''The Culgoora Catalog.
Data files are in tab-separated format from Vizier.
To download in the correct format, open a catalog online in Vizier,
select'Tab-Separated-Values' as the Output layout in the drop-down box, set
the maximum entries to 'unlimited', and click 'Sexagesimal' under the box
for 'Target Name or Position'.  Submit the query, and copy the output to a
txt file.  Copy this file to "culgoora.txt" in the _src directory of your AIPY
installation.'''
import aipy as a, numpy as np, os

class CulgooraCatalog(a.fit.SrcCatalog):
    metadata = {}
    def get_metadata(self):
        '''Return the actual 80 MHz and 160 MHz measurements used to compute
        the spectral index in the Culgoora catalog. Returns dictionary with
        source names linked to (S80, S160, index) triplets.  In cases where 
        measurements are unavailable, values are entered as None.'''
        return self.metadata
    def fromfile(self,filename):
        '''Parses a tab-delimited text file for the Culgoora array.  When both
        80 MHz and 160 MHz measurements are available, catalog entries use
        mfreq=.160 and include a computed spectral index.  When only one
        measurement is available, mfreq is set for that measurement, and index
        is set to 0.'''
        f = open(filename)
        addsrcs = []
        for L in [L for L in f.readlines() if not L.startswith('#')]:
            text = L.split('\t')
            if len(text) <= 4: continue
            try: int(text[0][0])
            except(ValueError): continue
            ra = text[0].replace(' ',':')
            dec = text[1].replace(' ',':')
            name = text[2].strip()
            try: jys160 = float(text[7])
            except(ValueError): jys160 = None
            try: jys080 = float(text[9])
            except(ValueError): jys080 = None
            try: index = float(text[10])
            except(ValueError): index = None
            # Catch flagged flux-densities
            if index is None and (not jys080 is None and not jys160 is None):
                if text[6].startswith('<'): jys160 = None
                if text[8].startswith('<'): jys080 = None
            self.metadata[name] = (jys080, jys160, index)
            if not index is None:
                addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                    jys=jys160, index=index, mfreq=0.160))
            elif not jys160 is None:
                addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                    jys=jys160, index=0, mfreq=0.160))
            else:
                addsrcs.append(a.fit.RadioFixedBody(ra, dec, name=name,
                    jys=jys080, index=0, mfreq=0.080))
        self.add_srcs(addsrcs)

CULGOORAFILE = os.path.dirname(__file__) + os.sep + 'culgoora.txt'
_culgooracat = None

def get_srcs(srcs=None, cutoff=None):
    global _culgooracat
    if _culgooracat is None:
        _culgooracat = CulgooraCatalog()
        _culgooracat.fromfile(CULGOORAFILE)
    if srcs is None:
        if cutoff is None: srcs = _culgooracat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _culgooracat.keys(): _culgooracat[s].update_jys(fq)
            srcs = [s for s in _culgooracat.keys() if _culgooracat[s].jys[0] > cut]
    srclist = []
    for s in srcs:
        try: srclist.append(_culgooracat[s])
        except(KeyError): pass
    return srclist
