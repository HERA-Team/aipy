'''This module interfaces to the Helmboldt catalog (http://arxiv.org/abs/0707.3418)'''

import aipy as a, numpy as np, os

class HelmboldtFixedBody(a.fit.RadioFixedBody):
    def compute(self, observer):
        a.phs.RadioFixedBody.compute(self, observer)
        self.update_jys(observer.get_afreqs())
    def update_jys(self, afreqs):
        A = np.log10(self._jys)
        try: B,C,D = (list(self.index) + [0,0,0])[:3]
        except(TypeError): B,C,D = (self.index,0,0)
        X = np.log10(afreqs / self.mfreq)
        self.jys = 10**(A + B*X + C*np.exp(D*X))
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'jys':      float(self._jys),
            'ra':       float(self._ra),
            'dec':      float(self._dec),
            'a1':       float(self.srcshape[0]),
            'a2':       float(self.srcshape[1]),
            'th':       float(self.srcshape[2]),
            'dra':      float(self.ionref[0]),
            'ddec':     float(self.ionref[1]),
        }
        try: aprms['index'] = list(self.index)
        except(TypeError): aprms['index'] = float(self.index)
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms

class HelmboldtCatalog(a.fit.SrcCatalog):
    metadata = {}
    def get_metadata(self):
        '''Return info from posfile about constituent measurements from which
        spectral fits were obtained.  Returns dictionary with source names
        linked to lists of (freq,flux,error) triplets. Must first run get_catalog
        to retrieve data.'''
        return self.metadata
    rms = {}
    def get_rms(self):
        """Return a dictionary giving the Helmboldt fit rms for each source.
        rms is only computed where enought data below 325MHz is available to 
        compute a spectral index or Kuehr fit. Null value is indicated with
        a rms=-1. Must first run get_catalog to retrieve data."""
        return self.rms
    ncomp = {}
    def get_ncomp(self):
        """Return a dictionary giving the Helmboldt fit rms for each source.
        rms is only computed where enought data below 325MHz is available to 
        compute a spectral index or Kuehr fit. Null value is indicated with
        a rms=-1. Must first run get_catalog to retrieve data."""
        return self.ncomp
    def fromfile(self, posfile, fitfile):
        srcs = {}
        rms = {}
        ncomp = {}
        # Read RA/DEC
        srclines = [L for L in open(posfile).readlines() if L.startswith('J')]
        for line in srclines:
            srcname = line[:9]
            srcs[srcname] = line[35:57]
            if not self.ncomp.has_key(srcname): self.ncomp[srcname] = int(line[33])
            if not self.metadata.has_key(srcname): self.metadata[srcname] = []
            md = (float(line[58:64])/1e3,float(line[65:73]),float(line[74:81]))
            self.metadata[srcname].append(md)
        for s in srcs:
            ra = srcs[s][:10].strip().replace(' ',':')
            dec = srcs[s][11:].strip().replace(' ',':')
            srcs[s] = [ra, dec]
        # Read spectral data
        srclines = [L for L in open(fitfile).readlines() if L.startswith('J')]
        for line in srclines: 
            srcs[line[:9]].append(map(float, line[13:62].split()))
            try: self.rms[line[:9]] = float(line[63:70])
            except(ValueError): self.rms[line[:9]]=-1

        addsrcs = []
        for s in srcs:
            ra,dec,spec = srcs[s]
            # If there's not a good fit on data, use VLSS value and default index
            if len(spec) < 5 or spec[3] == -99.:
                srctype = a.fit.RadioFixedBody
                jys = spec[0]
                # If there is no index put bogus value of -99
                try: index = spec[1]
                except(IndexError): index = -99
            else:
                srctype = HelmboldtFixedBody
                ABCD = spec[3:]
                jys,index = 10**ABCD[0], ABCD[1:]
            addsrcs.append(srctype(ra, dec, name=s, 
                jys=jys, index=index, mfreq=.074))
        self.add_srcs(addsrcs)

FITFILE = os.path.dirname(__file__) + os.sep + 'helm_fit.txt'
POSFILE = os.path.dirname(__file__) + os.sep + 'helm_pos.txt'

_helmcat = None

def get_srcs(srcs=None, cutoff=None):
    # Mechanism for delaying instantiation of catalog until it is accessed
    global _helmcat
    if _helmcat is None:
        _helmcat = HelmboldtCatalog()
        _helmcat.fromfile(POSFILE, FITFILE)
    if srcs is None:
        if cutoff is None: srcs = _helmcat.keys()
        else:
            cut, fq = cutoff
            fq = np.array([fq])
            for s in _helmcat.keys(): _helmcat[s].update_jys(fq)
            srcs = [s for s in _helmcat.keys() if _helmcat[s].jys[0] > cut]
    srclist = []
    for s in srcs:
        try: srclist.append(_helmcat[s])
        except(KeyError): pass
    return srclist
