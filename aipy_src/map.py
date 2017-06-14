"""
Module for mapping and modeling the entire sky.
"""

import numpy as np, healpix, coord, random,img
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits

# Set a fixed random seed to make scrambling deterministic
random.seed(1)

deg2rad = np.pi / 180.
rad2deg = 180. / np.pi

def pack_sphere(N):
    """Clever formula for putting N points nearly equally spaced on 
    the sphere.  Return xyz coordinates of each point in order from S to N."""
    dz = 2. / N
    z = np.arange(-1+dz/2,1, dz)
    r = np.sqrt(1-z**2)
    dL = np.pi * (3 - np.sqrt(5))
    long = np.arange(0, dL * N, dL)
    return np.array([r*np.cos(long), r*np.sin(long), z])

def _bit_reverse(N, nbits=None):
    """Numerically bit-reverse the number N assuming the specified number
    of bits.  In not specified, will infer the least number of bits that
    can represent N."""
    if nbits is None: nbits = int(np.floor(np.log2(N))) + 1
    ans = 0
    for bit in range(nbits):
        ans += np.bitwise_and(N, 2**bit) * 2**(nbits-2*bit-1)
    return ans

def _bit_reverse_order(N):
    """Generate a list of indices from 0 to N-1 in bit reversed order (or
    the nearest approximation if N is not a power of 2)."""
    nbits = int(np.floor(np.log2(N))) + 1
    indices = _bit_reverse(np.arange(2**nbits), nbits=nbits)
    return indices.compress(indices < N)

def _local_shuffle(L, width=2):
    """Shuffle elements of L without moving them too much by choosing chunks
    of the specified width and only scrambling elements within those chunks."""
    for i in range(int(np.ceil(len(L) / float(width)))):
        chunk = L[width*i:width*(i+1)]
        random.shuffle(chunk)
        L[width*i:width*(i+1)] = chunk

def facet_centers(N, ncrd=2):
    """Return the coordinates of N points equally spaced around the sphere.
    Will return xyz or ra,dec depending on ncrd.  Shuffles the order of the
    pointing centers wrt pack_sphere so that widely spaced points are done 
    first, and then the gaps between them, and then the gaps between those.."""
    assert(ncrd == 2 or ncrd == 3)
    ind1 = np.arange(N); _local_shuffle(ind1)
    ind2 = _bit_reverse_order(N)
    ind = ind1.take(ind2)
    pnts = pack_sphere(N)
    pnts = pnts.take(ind, axis=1)
    if ncrd == 3: return pnts
    else: return coord.eq2radec(pnts)

class Map(object):
    def __init__(self, *args, **kwargs):
        fromfits = kwargs.pop('fromfits', None)
        nindices = kwargs.pop('nindices', None)
        self.args,self.kwargs = args,kwargs
        self.map = healpix.HealpixMap(*args, **kwargs)
        self.wgt = healpix.HealpixMap(*args, **kwargs)
        self.ind = []
        if not fromfits is None: self.from_fits(fromfits)
        self.set_nindices(nindices)
    def set_interpol(self, onoff):
        self.map.set_interpol(onoff)
        self.wgt.set_interpol(onoff)
        for i in self.ind: i.set_interpol(onoff)
    def set_nindices(self, nindices):
        if nindices is None: return
        self.ind = self.ind[:nindices]
        for i in range(nindices - len(self.ind)):
            self.ind.append(healpix.HealpixMap(*self.args, **self.kwargs))
    def get(self, crds):
        fluxes = self.map[crds]
        wgts = self.wgt[crds]
        inds = [i[crds] for i in self.ind]
        return (wgts, fluxes, inds)
    def __getitem__(self, crds):
        """Return the average map/index values at the specified coordinates."""
        w = self.wgt[crds]
        w = np.where(w > 0, w, 1)
        fluxes = self.map[crds] / w
        ind = [i[crds] / w for i in self.ind]
        if len(ind) == 0: return fluxes
        return (fluxes, ind)
    def add(self, crds, wgts, fluxes, inds=[]):
        self.wgt[crds] += wgts
        self.map[crds] += fluxes * wgts
        for n,i in enumerate(inds): self.ind[n][crds] += i * wgts
    def put(self, crds, wgts, fluxes, inds=[]):
        self.wgt[crds] = wgts
        self.map[crds] = fluxes * wgts
        for n,i in enumerate(inds): self.ind[n][crds] = i * wgts
    def reset_wgt(self, wgt=1):
        w = np.where(self.wgt.map > 0, self.wgt.map, 1)
        self.map.map /= w
        for i in self.ind: i.map /= w
        self.wgt.map = np.where(self.wgt.map > 0, wgt, 0)
    def from_map(self, map):
        """Initialize this Map with data from another."""
        self.map.from_hpm(map.map)
        self.wgt.from_hpm(map.wgt)
        for i in range(min(len(self.ind), len(map.ind))): 
            self.ind[i].from_hpm(map.ind[i])
    def __getattr__(self, attr):
        try: object.__getatr__(self, attr)
        except(AttributeError): return self.map.__getattribute__(attr)
    def from_fits(self, filename, hdunum=1):
        self.map.from_fits(filename, hdunum=hdunum, colnum=0)
        self.args = ()
        self.kwargs = {'nside':self.nside(), 'scheme':self.scheme()}
        try:
            self.wgt.from_fits(filename, hdunum=hdunum, colnum=1)
        except(IndexError):
            self.wgt = healpix.HealpixMap(*self.args, **self.kwargs)
            self.wgt.set_map(np.ones_like(self.wgt.map))
        # Figure out how many cols there are (is there a better way?)
        hdu = pyfits.open(filename)[hdunum]
        nind = 0
        while True:
            try: i = hdu.data.field(2 + nind)
            except(IndexError): break
            nind += 1
        del(hdu)
        # Make a spectral index HPM for each additional col in fits file
        for i in range(nind):
            h = healpix.HealpixMap(*self.args, **self.kwargs)
            h.from_fits(filename, hdunum=hdunum, colnum=2+i)
            self.ind.append(h)
    def to_fits(self, filename, format=None, clobber=False,history=''):
        if format is None:
            format = healpix.default_fits_format_codes[self.get_dtype().type]
        hdu0 = pyfits.PrimaryHDU()
        col0 = pyfits.Column(name='signal', format=format, array=self.map.map)
        col1 = pyfits.Column(name='weights', format=format, array=self.wgt.map)
        col_inds = [pyfits.Column(name='sp_index%d'%n,format=format,array=i.map)
            for n,i in enumerate(self.ind)]
        cols = pyfits.ColDefs([col0, col1] + col_inds)
        tbhdu = pyfits.new_table(cols)
        self.map._set_fits_header(tbhdu.header)
        hdulist = pyfits.HDUList([hdu0, tbhdu])
        if history!='':
            history = [h.strip() for h in history.split("\n")]
            for line in history:
                if len(line)>1:
                    if line.startswith('#'):
                        for subline in img.word_wrap(line,80,0,0,'').split("\n"):
                            hdulist[0].header.add_history(subline)
                    else:
                        for subline in img.word_wrap(line,70,5,10,'#').split("\n"):
                            hdulist[0].header.add_history(subline)       
        hdulist.writeto(filename, clobber=clobber)
