"""
A module for mapping and modeling the entire sky.

Author: Aaron Parsons
Date: 11/29/06
Revisions:
"""

import numpy as n, healpix, pyfits

deg2rad = n.pi / 180.
rad2deg = 180. / n.pi

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
        w = n.where(w > 0, w, 1)
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
        w = n.where(self.wgt.map > 0, self.wgt.map, 1)
        self.map.map /= w
        for i in self.ind: i.map /= w
        self.wgt.map = n.where(self.wgt.map > 0, wgt, 0)
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
            self.wgt.set_map(n.ones_like(self.wgt.map))
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
    def to_fits(self, filename, format=None, clobber=False):
        if format is None:
            format = healpix.default_fits_format_codes[self.map.map.dtype.type]
        hdu0 = pyfits.PrimaryHDU()
        col0 = pyfits.Column(name='signal', format=format, array=self.map.map)
        col1 = pyfits.Column(name='weights', format=format, array=self.wgt.map)
        col_inds = [pyfits.Column(name='sp_index%d'%n,format=format,array=i.map)
            for n,i in enumerate(self.ind)]
        cols = pyfits.ColDefs([col0, col1] + col_inds)
        tbhdu = pyfits.new_table(cols)
        self.map._set_fits_header(tbhdu.header)
        hdulist = pyfits.HDUList([hdu0, tbhdu])
        hdulist.writeto(filename, clobber=clobber)
