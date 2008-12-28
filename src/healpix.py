"""
Provides interfaces to Healpix_cxx, which was developed at the 
Max-Planck-Institut fuer Astrophysik and financially supported by the 
Deutsches Zentrum fuer Luft- und Raumfahrt (DLR).
Adds data to the HealpixBase class using numpy arrays, and interfaces to
FITS files using pyfits.

Author: Aaron Parsons
Date: 12/01/2007
Revisions:
"""

import numpy as n, utils, pyfits
from _healpix import HealpixBase

default_fits_format_codes = {
    n.bool_:'L', n.uint8:'B', n.int16:'I', n.int32:'J', n.int64:'K', 
    n.float32:'E', n.float64:'D', n.complex64:'C', n.complex128:'M'
}

class HealpixMap(HealpixBase):
    """Adds a data map to the HealpixBase infrastructure."""
    _use_interpol = True
    def set_interpolation(self, onoff):
        """Choose whether __getitem___ (=HealpixMap[crd]) returns an
        interpolated value or just nearest pixel.  Default upon creating a
        HealpixMap is to interpolate."""
        self._use_interpol = onoff
    def SetData(self, data, ordering="RING"):
        """Assign data to HealpixMap.map.  Infers Nside from # of pixels via
        Npix = 12 * Nside**2."""
        try:
            assert(data.ndim == 1)
            nside = self.npix2nside(data.shape[0])
        except(AssertionError,ValueError):
            raise ValueError("Data must be a 1 dim array with 12*N**2 entries.")
        self.SetNside(nside, ordering)
        self.map = data
    def set_ordering(self, ordering):
        """Reorder the pixels in map to be "RING" or "NEST" ordering."""
        assert(ordering in ["RING", "NEST"])
        if ordering == self.Scheme(): return
        newmap = n.zeros_like(self.map)
        i = n.arange(self.map.shape[0])
        i = self.nest_ring_conv(i, ordering)
        newmap[i] = self.map
        self.SetNside(self.Nside(), ordering)
        self.map = newmap
    def interpolated_get_val(self, crd):
        """Interpolate value of map at provided coordinates from neighboring
        pixels."""
        px, wgts = self.get_interpol(crd)
        return n.sum(self.map[px] * wgts, axis=-1)
    def interpolated_add_val(self, crd, data):
        """Add data into map at provided coordinates, dividing the data among
        adjacent pixels according to their weights."""
        px, wgts = self.get_interpol(crd)
        data.shape += (1,)
        wgts *= data
        data.shape = (data.size,)
        px.shape = (px.size, 1); wgts = wgts.flatten()
        utils.add2array(self.map, px, wgts)
    def __getitem__(self, crd):
        if type(crd) != n.ndarray or crd.ndim == 1: return self.map[crd]
        elif self._use_interpol: return self.interpolated_get_val(crd)
        else: return self.map[self.crd2px(crd)]
    def __setitem__(self, crd, val):
        if type(crd) != n.ndarray or crd.ndim == 1: self.map[crd] = val
        elif self._use_interpol: return self.interpolated_add_val(crd, val)
        else: self.map[self.crd2px(crd)] = val
    def from_fits(self, filename, hdunum=1, colnum=0):
        """Read a HealpixMap from the specified location in a fits file."""
        hdu = pyfits.open(filename)[hdunum]
        data = hdu.data.field(colnum)
        ordering = hdu.header['ORDERING'][:4]
        self.SetData(data, ordering=ordering)
    def _set_fits_header(self, hdr):
        hdr.update('PIXTYPE', 'HEALPIX', 'HEALPIX pixelisation')
        ordering = self.Scheme()
        if ordering == 'NEST': ordering == 'NESTED'
        hdr.update('ORDERING', ordering,
            'Pixel ordering scheme, either RING or NESTED')
        hdr.update('NSIDE', self.Nside(), 'Resolution parameter for HEALPIX')
        hdr.update('FIRSTPIX', 0, "First pixel # (0 based)")
        hdr.update('LASTPIX', self.Npix()-1, "Last pixel # (0 based)")
        hdr.update('INDXSCHM', 'IMPLICIT', "Indexing: IMPLICIT or EXPLICIT")
    def to_fits(self, filename, format=None):
        """Write a HealpixMap to a fits file in the fits format specified by
        'format'.  Default uses mapping of numpy types to fits format types
        stored in default_fits_format_codes."""
        if format is None:
            format = default_fits_format_codes[self.map.dtype.type]
        hdu0 = pyfits.PrimaryHDU()
        col0 = pyfits.Column(name='signal', format=format, array=self.map)
        cols = pyfits.ColDefs([col0])
        tbhdu = pyfits.new_table(cols)
        self._set_fits_header(tbhdu.header)
        hdulist = pyfits.HDUList([hdu0, tbhdu])
        hdulist.writeto(filename)
