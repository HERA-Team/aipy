"""
Provides interfaces to Healpix_cxx, which was developed at the 
Max-Planck-Institut fuer Astrophysik and financially supported by the 
Deutsches Zentrum fuer Luft- und Raumfahrt (DLR).
Adds data to the HealpixBase class using numpy arrays, and interfaces to
FITS files using astropy/pyfits.
"""

import numpy as np, utils
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits
from _healpix import HealpixBase
from _alm import Alm

default_fits_format_codes = {
    np.bool_:'L', np.uint8:'B', np.int16:'I', np.int32:'J', np.int64:'K', 
    np.float32:'E', np.float64:'D', np.complex64:'C', np.complex128:'M'
}

def mk_arr(val, dtype=np.double):
    if type(val) is np.ndarray: return val.astype(dtype)
    return np.array(val, dtype=dtype).flatten()

class HealpixMap(HealpixBase):
    """Collection of utilities for mapping data on a sphere.  Adds a data map 
    to the infrastructure in _healpix.HealpixBase."""
    def __init__(self, *args, **kwargs):
        dtype = kwargs.pop('dtype', np.double)
        interp = kwargs.pop('interp', False)
        fromfits = kwargs.pop('fromfits', None)
        HealpixBase.__init__(self, *args, **kwargs)
        self._use_interpol = interp
        if fromfits is None:
            m = np.zeros((self.npix(),), dtype=dtype)
            self.set_map(m, scheme=self.scheme())
        else: self.from_fits(fromfits)
    def set_interpol(self, onoff):
        """Choose whether __getitem___ (i.e. HealpixMap[crd]) returns an
        interpolated value or just nearest pixel.  Default upon creating a
        HealpixMap is to not interpolate."""
        self._use_interpol = onoff
    def set_map(self, data, scheme="RING"):
        """Assign data to HealpixMap.map.  Infers Nside from # of pixels via
        Npix = 12 * Nside**2."""
        try:
            assert(data.ndim == 1)
            nside = self.npix2nside(data.shape[0])
        except(AssertionError,ValueError):
            raise ValueError("Data must be a 1 dim array with 12*N**2 entries.")
        self.set_nside_scheme(nside, scheme)
        self.map = data
    def get_map(self):
        """Return Healpix data as a 1 dimensional numpy array."""
        return self.map
    def change_scheme(self, scheme):
        """Reorder the pixels in map to be "RING" or "NEST" ordering."""
        assert(scheme in ["RING", "NEST"])
        if scheme == self.scheme(): return
        i = self.nest_ring_conv(np.arange(self.npix()), scheme)
        self[i] = self.map
        self.set_nside_scheme(self.nside(), scheme)
    def __getitem__(self, crd):
        """Access data on a sphere via hpm[crd].
        crd = either 1d array of pixel indices, (th,phi), or (x,y,z), where
        th,phi,x,y,z are numpy arrays of coordinates."""
        if type(crd) is tuple:
            crd = [mk_arr(c, dtype=np.double) for c in crd]
            if self._use_interpol:
                px,wgts = self.crd2px(*crd, **{'interpolate':1})
                return np.sum(self.map[px] * wgts, axis=-1)
            else: px = self.crd2px(*crd)
        else: px = mk_arr(crd, dtype=np.long)
        return self.map[px]
    def __setitem__(self, crd, val):
        """Assign data to a sphere via hpm[crd] = val.  Functionality slightly
        complicated to make repeat coordinates assign sum of values (i.e.
        crd = ([1,1], [2,2]), val = [3,3] will assign 6 to location (1,2).
        crd = either 1d array of pixel indices, (th,phi), or (x,y,z), where
        th,phi,x,y,z are numpy arrays of coordinates."""
        if type(crd) is tuple:
            crd = [mk_arr(c, dtype=np.double) for c in crd]
            px = self.crd2px(*crd)
        else:
            if type(crd) is np.ndarray: assert(len(crd.shape) == 1)
            px = mk_arr(crd, dtype=np.int)
        if px.size == 1:
            if type(val) is np.ndarray: val = mk_arr(val, dtype=self.map.dtype)
            self.map[px] = val
        else:
            m = np.zeros_like(self.map)
            px = px.reshape(px.size,1)
            cnt = np.zeros(self.map.shape, dtype=np.bool)
            val = mk_arr(val, dtype=m.dtype)
            utils.add2array(m, px, val)
            utils.add2array(cnt, px, np.ones(val.shape, dtype=np.bool))
            self.map = np.where(cnt, m, self.map)
    def from_hpm(self, hpm):
        """Initialize this HealpixMap with data from another.  Takes care
        of upgrading or downgrading the resolution, and swaps ordering
        scheme if necessary."""
        if hpm.nside() < self.nside():
            interpol = hpm._use_interpol
            hpm.set_interpol(True)
            px = np.arange(self.npix())
            th,phi = self.px2crd(px, ncrd=2)
            self[px] = hpm[th,phi].astype(self.get_dtype())
        elif hpm.nside() > self.nside():
            px = np.arange(hpm.npix())
            th,phi = hpm.px2crd(px, ncrd=2)
            self[th,phi] = hpm[px].astype(self.get_dtype())
        else:
            if hpm.scheme() == self.scheme():
                self.map = hpm.map.astype(self.get_dtype())
            else:
                i = self.nest_ring_conv(np.arange(self.npix()), hpm.scheme())
                self.map = hpm.map[i].astype(self.get_dtype())
    def from_alm(self, alm):
        """Set data to the map generated by the spherical harmonic
        coefficients contained in alm."""
        self.set_map(alm.to_map(self.nside(), self.scheme()))
    def to_alm(self, lmax, mmax, iter=1):
        """Return an Alm object containing the spherical harmonic components
        of a map (in RING mode) up to the specified lmax,mmax.  Greater
        accuracy can be achieved by increasing iter."""
        assert(self.scheme() == 'RING')
        alm = Alm(lmax,mmax)
        alm.from_map(self.map, iter)
        return alm
    def from_fits(self, filename, hdunum=1, colnum=0):
        """Read a HealpixMap from the specified location in a fits file."""
        hdu = pyfits.open(filename)[hdunum]
        data = hdu.data.field(colnum)
        if not data.dtype.isnative:
            data.dtype = data.dtype.newbyteorder()
            data.byteswap(True)
        scheme= hdu.header['ORDERING'][:4]
        self.set_map(data, scheme=scheme)
    def _set_fits_header(self, hdr):
        hdr.update('PIXTYPE', 'HEALPIX', 'HEALPIX pixelisation')
        scheme = self.scheme()
        if scheme == 'NEST': scheme == 'NESTED'
        hdr.update('ORDERING', scheme,
            'Pixel ordering scheme, either RING or NESTED')
        hdr.update('NSIDE', self.nside(), 'Resolution parameter for HEALPIX')
        hdr.update('FIRSTPIX', 0, "First pixel # (0 based)")
        hdr.update('LASTPIX', self.npix()-1, "Last pixel # (0 based)")
        hdr.update('INDXSCHM', 'IMPLICIT', "Indexing: IMPLICIT or EXPLICIT")
    def get_dtype(self):
        return self.map.dtype
    def to_fits(self, filename, format=None, clobber=True):
        """Write a HealpixMap to a fits file in the fits format specified by
        'format'.  Default uses mapping of numpy types to fits format types
        stored in default_fits_format_codes."""
        if format is None:
            format = default_fits_format_codes[self.get_dtype().type]
        hdu0 = pyfits.PrimaryHDU()
        col0 = pyfits.Column(name='signal', format=format, array=self.map)
        cols = pyfits.ColDefs([col0])
        tbhdu = pyfits.new_table(cols)
        self._set_fits_header(tbhdu.header)
        hdulist = pyfits.HDUList([hdu0, tbhdu])
        hdulist.writeto(filename,clobber=clobber)

