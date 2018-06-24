"""
Provides interfaces to Healpix_cxx, which was developed at the
Max-Planck-Institut fuer Astrophysik and financially supported by the
Deutsches Zentrum fuer Luft- und Raumfahrt (DLR).
Adds data to the HealpixBase class using numpy arrays, and interfaces to
FITS files using astropy/pyfits.
"""

from __future__ import print_function, division, absolute_import

import numpy as np
from . import utils
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits

import healpy

default_fits_format_codes = {
    np.bool_:'L', np.uint8:'B', np.int16:'I', np.int32:'J', np.int64:'K',
    np.float32:'E', np.float64:'D', np.complex64:'C', np.complex128:'M'
}

def mk_arr(val, dtype=np.double):
    if type(val) is np.ndarray: return val.astype(dtype)
    return np.array(val, dtype=dtype).flatten()

HEALPIX_MODES = ('RING','NEST')

class HealpixBase(object):
    """Functionality related to the HEALPix pixelisation."""
    def __init__(self, nside=1, scheme='RING'):
        self._nside = nside
        self._scheme = scheme
    def npix2nside(self, npix):
        """Convert number of pixels to number of sides."""
        return healpy.npix2nside(npix)
    def nest_ring_conv(self, px, scheme):
        """Translate an array of pixel numbers to index data in the scheme specified in 'scheme' ('NEST' or 'RING').  Returns px."""
        mode = {'RING':healpy.nest2ring, 'NEST':healpy.ring2nest}
        if scheme != self._scheme: px = mode[scheme](self._nside, px)
        self._scheme = scheme
        return px
    def set_nside_scheme(self, nside=None, scheme=None):
        """Adjust Nside and Scheme ('RING' or 'NEST')."""
        if nside is not None:
            pow2 = np.log2(nside)
            assert(pow2 == np.around(pow2))
            self._nside = nside
        if scheme is not None:
            assert(scheme in HEALPIX_MODES)
            self._scheme = scheme
    def crd2px(self, c1, c2, c3=None, interpolate=False):
        """Convert 1 dimensional arrays of input coordinates to pixel indices. If only c1,c2 provided, then read them as th,phi.  If c1,c2,c3 provided, read them as x,y,z. If interpolate is False, return a single pixel coordinate.  If interpolate is True, return px,wgts where each entry in px contains the 4 pixels adjacent to the specified location, and wgt contains the 4 corresponding weights of those pixels."""
        is_nest = (self._scheme == 'NEST')
        if not interpolate:
            if c3 is None: # th/phi angle mode
                px = healpy.ang2pix(self._nside, c1, c2, nest=is_nest)
            else: # x,y,z mode
                px = healpy.vec2pix(self._nside, c1, c2, c3, nest=is_nest)
            return px
        else:
            if c3 is not None: # need to translate xyz to th/phi
                c1,c2 = healpy.vec2ang(np.array([c1,c2,c3]).T)
            px,wgts = healpy.get_interp_weights(self._nside, c1, c2, nest=is_nest)
            return px.T, wgts.T
    def px2crd(self, px, ncrd=3):
        """Convert a 1 dimensional input array of pixel numbers to the type of coordinates specified by ncrd.  If ncrd=3 (default), the returned array will have (x,y,z) for each pixel.  Otherwise if ncrd=2, the returned array will have (theta,phi) for each pixel."""
        is_nest = (self._scheme == 'NEST')
        assert(ncrd in (2,3))
        if ncrd == 2: # th/phi mode
            th,phi = healpy.pix2ang(self._nside, px, nest=is_nest)
            return th, phi
        else: # ncrd == 3 -> vec mode
            x,y,z = healpy.pix2vect(self._nside, px, nest=is_nest)
    def order(self):
        """Return the order parameter."""
        return healpy.nside2order(self._nside)
    def nside(self):
        """Return the Nside parameter."""
        return self._nside
    def npix(self):
        """Return the number of pixels in the map."""
        return healpy.nside2npix(self._nside)
    def scheme(self):
        """Return the scheme of the map ('NEST' or 'RING')."""
        return self._scheme

class Alm(object):
    '''Holds coefficients for spherical harmonics up to the specified order, and generates a real-valued HealpixMap from them.  Individual (l,m) coefficients can be accessed by alm[l,m].'''
    def __init__(self, lmax, mmax, dtype=np.complex128):
        assert(lmax >= mmax)
        self._alm = healpy.Alm()
        self._lmax = lmax
        self._mmax = mmax
        self.dtype = dtype
        self.set_to_zero()
    def size(self):
        '''Return the number of Alm coefficients.'''
        return self._alm.getsize(self._lmax,self._mmax)
    def set_to_zero(self):
        """Clear all coefficients."""
        data = np.zeros(self.size(), dtype=self.dtype)
        self.set_data(data)
    def lmax(self):
        """Return the maximum L."""
        return self._lmax
    def mmax(self):
        """Return the maximum M."""
        return self._mmax
    def __getitem__(self, lm):
        l,m = lm
        i = self._alm.getidx(self._lmax, l, m)
        return self.data[i]
    def __setitem__(self, lm, val):
        l,m = lm
        i = self._alm.getidx(self._lmax, l, m)
        self.data[i] = val
    def to_map(self, nside, pixwin=False, fwhm=0.0, sigma=None, 
            pol=True, verbose=False):
        '''Return data for the Healpix map in 'RING' mode with the specified 
        nside (power of 2) generated by these Alm coefficients.

        Parameters
        ----------
        nside : int, scalar
          The nside of the output map.
        pixwin : bool, optional
          Smooth the alm using the pixel window functions. Default: False.
        fwhm : float, scalar, optional
          The fwhm of the Gaussian used to smooth the map (applied on alm)
          [in radians]
        sigma : float, scalar, optional
          The sigma of the Gaussian used to smooth the map (applied on alm)
          [in radians]
        pol : bool, optional
          If True, assumes input alms are TEB. Output will be TQU maps.
          (input must be 1 or 3 alms)
          If False, apply spin 0 harmonic transform to each alm.
          (input can be any number of alms)
          If there is only one input alm, it has no effect. Default: True.

        Returns
        -------
        maps : array or list of arrays
          A Healpix map in RING scheme at nside or a list of T,Q,U maps (if
          polarized input)'''
        return healpy.alm2map(self.get_data(), nside, 
            lmax=self._lmax, mmax=self._mmax,
            pixwin=pixwin, fwhm=fwhm, sigma=sigma, pol=pol, verbose=verbose)
    def from_map(self, data, iter=3, pol=True, use_weights=False, gal_cut=0):
        '''Set the coefficients of this Alm object (with its specified 
        lmax and mmax) from the data of a HealpixMap in 'RING' mode.

        Parameters
        ----------
        data : array-like, shape (Npix,) 
          The input map.
        iter : int, scalar, optional
          Number of iteration (default: 3)
        pol : bool, optional
          If True, assumes input maps are TQU. Output will be TEB alm's.
          (input must be 1 or 3 maps)
          If False, apply spin 0 harmonic transform to each map.
          (input can be any number of maps)
          If there is only one input map, it has no effect. Default: True.
        use_weights: bool, scalar, optional
          If True, use the ring weighting. Default: False.
        gal_cut : float [degrees]
          pixels at latitude in [-gal_cut;+gal_cut] are not taken into account'''
        data = healpy.map2alm(data, lmax=self._lmax, mmax=self._mmax, iter=iter,
            pol=pol, use_weights=use_weights, gal_cut=gal_cut)
        self.set_data(data)
    def lm_indices(self):
        """Return the L and M indices of the coefficients contained in Alm, in the order that they are returned by get_data()."""
        return np.array([self._alm.getlm(self._lmax, i) for i in range(self.size)])
    def get_data(self):
        """Return all of the coefficients contained in Alm, in the order indexed by lm_indices()."""
        return self.data
    def set_data(self, data):
        """Set all of coefficients contained in Alm, in the order indexed by lm_indices()."""
        assert(data.size == self.size())
        self.data = data.astype(self.dtype)

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
