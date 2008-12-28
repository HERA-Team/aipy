"""
A module for gridding UVW data (including W projection), forming images,
and combining (mosaicing) images into a spherical map.

Author: Aaron Parsons
Date: 11/29/06
Revisions:
    10/18/07    arp     Changed calling syntax inside Img.put for an impressive
                        10x speedup.  Also changed Basemap call to not include
                        earth map, for another hefty speed-up.
    11/01/07    arp     Corrected bug for Img.put whereby data in same put call
                        clobbered other data (instead of adding together).
                        Moved put loop into C++ in utils module for speed +
                        correction.
    12/19/07    arp     Removed dependence on Basemap for ra/dec coordinates
                        in image.  Fixed W proj kernel to use correct l,m
                        coordinates.
"""

import numpy as n, utils, healpix, pyfits, coord

deg2rad = n.pi / 180.
rad2deg = 180. / n.pi

def recenter(a, c):
    """Slide the (0,0) point of matrix a to a new location tuple c.  This is
    useful for making an image centered on your screen after performing an
    inverse fft of uv data."""
    s = a.shape
    c = (c[0] % s[0], c[1] % s[1])
    if n.ma.isMA(a):
        a1 = n.ma.concatenate([a[c[0]:], a[:c[0]]], axis=0)
        a2 = n.ma.concatenate([a1[:,c[1]:], a1[:,:c[1]]], axis=1)
    else:
        a1 = n.concatenate([a[c[0]:], a[:c[0]]], axis=0)
        a2 = n.concatenate([a1[:,c[1]:], a1[:,:c[1]]], axis=1)
    return a2

def convolve2d(a, b):
    """Convolve a and b by multiplying in Fourier domain.  Must be same size."""
    return n.fft.ifft2(n.fft.fft2(a) * n.fft.fft2(b))

def gaussian_beam(sigma, shape=0, amp=1., center=(0,0)):
    """Return a 2D gaussian.  Normalized to area under curve = 'amp'.  
    Down by 1/e at distance 'sigma' from 'center'."""
    if type(shape) == type(0): shape = array([2, 2]) * sigma
    def gaussian(x, y):
        nx = n.where(x > shape[0] / 2, x - shape[0], x)
        ny = n.where(y > shape[1] / 2, y - shape[1], y)
        return n.exp(-(nx**2 + ny**2) / sigma**2)
    g = n.fromfunction(gaussian, shape)
    g /= g.sum() / amp
    return recenter(g, center)

class Img:
    """A class for gridding uv data, recording the synthesized beam profile,
    and performing the transforms into image domain."""
    def __init__(self, size=100, res=1):
        """size: The number of wavelengths which the UV matrix spans (this 
            determines the image resolution).
        res: The resolution of the UV matrix (determines image FOV)."""
        self.res = float(res)
        self.size = float(size)
        dim = n.round(self.size / self.res)
        self.uv = n.zeros(shape=(dim,dim), dtype=n.complex64)
        self.bm = n.zeros(shape=(dim,dim), dtype=n.complex64)
    def get_LM(self, center=(0,0)):
        """Get the (l,m) image coordinates for an inverted UV matrix."""
        dim = self.uv.shape[0]
        M,L = n.indices(self.uv.shape)
        L,M = n.where(L > dim/2, L-dim, L), n.where(M > dim/2, M-dim, M)
        L,M = L.astype(n.float)/dim/self.res, M.astype(n.float)/dim/self.res
        #L,M = n.fliplr(L), n.fliplr(M)
        mask = n.where(L**2 + M**2 >= 1, 1, 0)
        L,M = n.ma.array(L, mask=mask), n.ma.array(M, mask=mask)
        return recenter(L, center), recenter(M, center)
    def put(self, uvw, data, wgts=None, apply=True):
        """Grid uv data (w is ignored) onto a UV plane.  Data should already
        have the phase due to w removed.  Assumes the Hermitian conjugate
        data is in uvw already (i.e. the conjugate points are not placed for
        you).  If wgts are not supplied, default is 1 (normal weighting).
        If apply is false, returns uv and bm data without applying it do
        the internally stored matrices."""
        if apply:
            uv = self.uv
            bm = self.bm
        else:
            uv = n.zeros_like(self.uv)
            bm = n.zeros_like(self.bm)
        if wgts is None: wgts = n.ones_like(data)
        inds = n.round(uvw[:,:2] / self.res).astype(n.int)
        ok = n.logical_and(n.abs(inds[:,0]) < uv.shape[0],
            n.abs(inds[:,1]) < uv.shape[1])
        data = data.compress(ok)
        wgts = wgts.compress(ok)
        inds = inds.compress(ok, axis=0)
        utils.add2array(uv, inds, data.astype(uv.dtype))
        utils.add2array(bm, inds, wgts.astype(bm.dtype))
        if not apply: return uv, bm
    def append_hermitian(self, uvw, data, wgts=None):
        """Append to (uvw, data, [wgts]) the points (-uvw, conj(data), [wgts]).
        This is standard practice to get a real-valued image."""
        uvw = n.concatenate([uvw, -uvw], axis=0)
        data = n.concatenate([data, n.conj(data)], axis=0)
        if wgts is None: return uvw, data
        wgts = n.concatenate([wgts, wgts], axis=0)
        return uvw, data, wgts
    def image(self, center=(0,0)):
        """Return the inverse FFT of the UV matrix, with the 0,0 point moved
        to 'center'."""
        return recenter(n.abs(n.fft.ifft2(self.uv)), center)
    def bm_image(self, center=(0,0)):
        """Return the inverse FFT of the sample weightings, with the 0,0 point
        moved to 'center'."""
        return recenter(n.abs(n.fft.ifft2(self.bm)), center)
    def get_coords(self, ra=0, dec=0, center=(0,0)):
        """Return the ra, dec coordinates of each pixel in the image, assuming
        the image is centered on the provided ra, dec (in radians)."""
        y,z = self.get_LM(center)
        shape, mask = y.shape, y.mask
        x = n.sqrt(1 - y**2 - z**2)
        x,y,z = x.flatten(), y.flatten(), -z.flatten()
        vec = n.array([x.filled(0),y.filled(0),z.filled(0)])
        rot_ra = coord.rot_m(ra, n.array([0,0,1]))
        rot_dec = coord.rot_m(dec, n.array([0,1,0]))
        rot = n.dot(rot_ra, rot_dec) # weird notation: rotate ra, then dec
        vec = n.dot(rot, vec)
        th,ra = coord.xyz2thphi(vec)
        dec = n.pi/2 - th
        ra.shape, dec.shape = shape, shape
        return n.ma.array(-ra, mask=mask), n.ma.array(dec, mask=mask)
    def get_top(self, center=(0,0)):
        x,y = self.get_LM(center)
        z = n.sqrt(1 - x**2 - y**2)
        return x,y,z
    def get_eq(self, ha=0, dec=0, center=(0,0)):
        x,y,z = self.get_top(center)
        shape,mask = x.shape, x.mask
        vec = n.array([a.filled().flatten() for a in (x,y,z)])
        m = coord.top2eq_m(ha, dec)
        vec = n.dot(m, vec)
        vec.shape = (3,) + shape
        return n.ma.array(vec, mask=[mask,mask,mask])

class ImgW(Img):
    """A subclass of Img adding W projection functionality (see Cornwell
    et al. 2005 "Widefield Imaging Problems in Radio Astronomy")."""
    def __init__(self, size=100, res=1, wres=.5):
        """wres: the gridding resolution of sqrt(w) when projecting to w=0."""
        Img.__init__(self, size=size, res=res)
        self.wres = wres
    def put(self, uvw, data, wgts=None):
        """Same as Img.put, only now the w component is projected to the w=0
        plane before applying the data to the UV matrix."""
        if len(uvw) == 0: return
        if wgts is None: wgts = n.ones_like(data)
        # Sort uvw in order of w
        order = n.argsort(uvw[:,-1])
        uvw = uvw.take(order, axis=0)
        data = data.take(order)
        w = uvw[:,-1]
        sqrt_w = n.sqrt(n.abs(w)) * n.sign(w)
        i = 0
        while True:
            # Grab a chunk of uvw's that grid w to same point.
            j = sqrt_w.searchsorted(sqrt_w[i]+self.wres)
            avg_w = n.average(w[i:j])
            # Put all uv's down on plane for this gridded w point
            uv, bm = Img.put(self, uvw[i:j,:],data[i:j],wgts[i:j],apply=False)
            # Convolve with the W projection kernel
            ker = n.fromfunction(lambda u,v: self.conv_ker(u,v,avg_w),
                uv.shape)
            self.uv += convolve2d(uv, ker)
            self.bm += convolve2d(bm, ker)
            if j >= len(w): break
            i = j
    def conv_ker(self, u, v, w):
        """Generates the W projection kernel (a function of u, v) for the
        supplied value of w.  See Cornwell et al. 2005 "Widefield Imaging
        Problems in Radio Astronomy" for discussion.  This implementation
        uses a numerically evaluated Fresnel kernel, rather than the
        small-angle approximated one given in the literature."""
        L,M = self.get_LM()
        # This is the exactly evaluated kernel (works better)
        sqrt = n.sqrt(1 - L**2 - M**2)
        G = n.exp(-2*n.pi*1j*w*(sqrt - 1))
        # This is the kernel described by Cornwell using the small angle approx.
        #G = n.exp(n.pi*1j*w*(l**2 + m**2))
        G_hat = n.fft.fft2(G.filled(0))
        return G_hat

class SkyMap:
    """A class for combining data from multiple pointings into a map of the
    sky in cylindrical coordinates."""
    def __init__(self, res=.01, fromfile=None):
        """res: The map resolution, in radians
        fromfile: Initialize from an existing SkyMap file."""
        if not fromfile is None: self.fromfile(fromfile)
        else:
            self.res = float(res)
            d = int(n.pi / self.res)
            self.map = n.zeros((2*d, d), dtype=n.float)
            self.wgt = n.zeros((2*d, d), dtype=n.float)
    def add_data(self, ras, decs, data, weights=None):
        """Add data to the map at the specified ra, dec location.  If
        specified, weights can be given to the data.  Otherwise, equal (=1)
        weighting is used.  Coordinates are in radians."""
        if weights is None: weights = n.ones_like(data)
        ras = ras % (2*n.pi)
        decs = (decs + n.pi/2) % n.pi
        ras = n.round(ras / self.res)
        decs = n.round(decs / self.res)
        data *= weights
        for r,d,dat,wgt in zip(ras.flatten(), decs.flatten(), 
                data.flatten(), weights.flatten()):
            try:
                self.map[r,d] += dat
                self.wgt[r,d] += wgt
            except(IndexError): pass
    def fromfile(self, filename):
        """Read a SkyMap from a file."""
        f = open(filename)
        data = n.fromfile(f)
        self.map = data[:data.size/2]
        self.wgt = data[data.size/2:]
        f.close()
        d = int(n.round(n.sqrt(self.map.size/2)))
        self.map.shape = (2*d, d)
        self.wgt.shape = (2*d, d)
        self.res = n.pi / d
    def tofile(self, filename):
        """Write this SkyMap to a file."""
        f = open(filename, 'w')
        self.map.tofile(f)
        self.wgt.tofile(f)
        f.close()
    def get_map(self):
        """Average the data from all the different pointings together and
        return a single image in cylindrical coordinates.  Axis 0 is
        right ascension [0,2pi); axis 1 will be declination [-pi,pi)."""
        w = n.where(self.wgt > 0, self.wgt, 1)
        return n.abs(self.map / w)

class SkyHMap:
    def __init__(self, nside=128, ordering='RING', fromfits=None):
        self.map = healpix.HealpixMap(nside, ordering=ordering)
        self.wgt = healpix.HealpixMap(nside, ordering=ordering)
        if not fromfits is None: self.from_fits(fromfits)
        else:
            m = n.zeros(self.map.Npix(), dtype=n.float)
            w = n.zeros(self.wgt.Npix(), dtype=n.float)
            self.map.SetData(m, ordering=ordering)
            self.wgt.SetData(w, ordering=ordering)
    def __setitem__(self, crds, value):
        data, weights = value
        self.map[crds] += data * weights
        self.wgt[crds] += weights
    def __getitem__(self, crds):
        m = self.map[crds]
        w = self.wgt[crds]
        return m / n.where(w > 0, w, 1)
    def from_fits(self, filename, hdunum=1):
        self.map.from_fits(filename, hdunum=hdunum, colnum=0)
        self.wgt.from_fits(filename, hdunum=hdunum, colnum=1)
    def to_fits(self, filename, format=None, clobber=False):
        if format is None:
            format = healpix.default_fits_format_codes[self.map.map.dtype.type]
        hdu0 = pyfits.PrimaryHDU()
        col0 = pyfits.Column(name='signal', format=format, array=self.map.map)
        col1 = pyfits.Column(name='weights', format=format, array=self.wgt.map)
        cols = pyfits.ColDefs([col0, col1])
        tbhdu = pyfits.new_table(cols)
        self.map._set_fits_header(tbhdu.header)
        hdulist = pyfits.HDUList([hdu0, tbhdu])
        hdulist.writeto(filename, clobber=clobber)
