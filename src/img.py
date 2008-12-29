"""
Module for gridding UVW data (including W projection), forming images,
and combining (mosaicing) images into spherical maps.

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
    02/25/08    arp     Flipped returned image/beam to correspond to returned
                        coordinates.
"""

import numpy as n, utils, healpix, coord

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
    """Class for gridding uv data, recording the synthesized beam profile,
    and performing transforms into image domain."""
    def __init__(self, size=100, res=1, mf_order=0):
        """size = number of wavelengths which the UV matrix spans (this 
        determines the image resolution).
        res = resolution of the UV matrix (determines image field of view)."""
        self.res = float(res)
        self.size = float(size)
        dim = n.round(self.size / self.res)
        self.shape = (dim,dim)
        self.uv = n.zeros(shape=self.shape, dtype=n.complex64)
        self.bm = []
        for i in range(mf_order+1):
            self.bm.append(n.zeros(shape=self.shape, dtype=n.complex64))
    def get_LM(self, center=(0,0)):
        """Get the (l,m) image coordinates for an inverted UV matrix."""
        dim = self.shape[0]
        M,L = n.indices(self.shape)
        L,M = n.where(L > dim/2, dim-L, -L), n.where(M > dim/2, dim-M, -M)
        L,M = L.astype(n.float)/dim/self.res, M.astype(n.float)/dim/self.res
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
        if len(self.bm) == 1 and len(wgts) != 1: wgts = [wgts]
        assert(len(wgts) == len(self.bm))
        if apply: uv,bm = self.uv,self.bm
        else:
            uv = n.zeros_like(self.uv)
            bm = [n.zeros_like(i) for i in self.bm]
        if wgts is None:
            wgts = []
            for i in range(len(self.bm)):
                if i == 0: wgts.append(n.ones_like(data))
                else: wgts.append(n.zeros_like(data))
        inds = n.round(uvw[:,:2] / self.res).astype(n.int)
        ok = n.logical_and(n.abs(inds[:,0]) < self.shape[0],
            n.abs(inds[:,1]) < self.shape[1])
        data = data.compress(ok)
        inds = inds.compress(ok, axis=0)
        utils.add2array(uv, inds, data.astype(uv.dtype))
        for i,w in enumerate(wgts):
            w = w.compress(ok)
            utils.add2array(bm[i], inds, w.astype(bm[0].dtype))
        if not apply: return uv, bm
    def append_hermitian(self, uvw, data, wgts=None):
        """Append to (uvw, data, [wgts]) the points (-uvw, conj(data), [wgts]).
        This is standard practice to get a real-valued image."""
        uvw = n.concatenate([uvw, -uvw], axis=0)
        data = n.concatenate([data, n.conj(data)], axis=0)
        if wgts is None: return uvw, data
        if len(self.bm) == 1 and len(wgts) != 1: wgts = [wgts]
        assert(len(wgts) == len(self.bm))
        for i,w in enumerate(wgts): wgts[i] = n.concatenate([w, w], axis=0)
        return uvw, data, wgts
    def _gen_img(self, data, center=(0,0)):
        """Return the inverse FFT of the provided data, with the 0,0 point 
        moved to 'center'.  Tranposes to put up=North, right=East."""
        return recenter(n.abs(n.fft.ifft2(data.transpose())), center)
    def image(self, center=(0,0)):
        """Return the inverse FFT of the UV matrix, with the 0,0 point moved
        to 'center'.  Tranposes to put up=North, right=East."""
        return self._gen_img(self.uv, center=center)
    def bm_image(self, center=(0,0), term=None):
        """Return the inverse FFT of the sample weightings (for all mf_order
        terms, or the specified term if supplied), with the 0,0 point
        moved to 'center'.  Tranposes to put up=North, right=East."""
        if not term is None: return self._gen_img(self.bm[term], center=center)
        else: return [self._gen_img(b, center=center) for b in self.bm]
    def get_top(self, center=(0,0)):
        """Return the topocentric coordinates of each pixel in the image."""
        x,y = self.get_LM(center)
        z = n.sqrt(1 - x**2 - y**2)
        return x,y,z
    def get_eq(self, ra=0, dec=0, center=(0,0)):
        """Return the equatorial coordinates of each pixel in the image, 
        assuming the image is centered on the provided ra, dec (in radians)."""
        x,y,z = self.get_top(center)
        shape,mask = x.shape, x.mask
        vec = n.array([a.filled().flatten() for a in (x,y,z)])
        m = coord.top2eq_m(-ra, dec)
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
    def put(self, uvw, data, wgts=None, invker2=None):
        """Same as Img.put, only now the w component is projected to the w=0
        plane before applying the data to the UV matrix."""
        if len(uvw) == 0: return
        if len(self.bm) == 1 and len(wgts) != 1: wgts = [wgts]
        assert(len(wgts) == len(self.bm))
        if wgts is None:
            wgts = []
            for i in range(len(self.bm)):
                if i == 0: wgts.append(n.ones_like(data))
                else: wgts.append(n.zeros_like(data))
        # Sort uvw in order of w
        order = n.argsort(uvw[:,-1])
        uvw = uvw.take(order, axis=0)
        data = data.take(order)
        wgts = [wgt.take(order) for wgt in wgts]
        w = uvw[:,-1]
        sqrt_w = n.sqrt(n.abs(w)) * n.sign(w)
        i = 0
        while True:
            # Grab a chunk of uvw's that grid w to same point.
            j = sqrt_w.searchsorted(sqrt_w[i]+self.wres)
            avg_w = n.average(w[i:j])
            # Put all uv's down on plane for this gridded w point
            wgtsij = [wgt[i:j] for wgt in wgts]
            uv, bm = Img.put(self, uvw[i:j,:],data[i:j],wgtsij,apply=False)
            # Convolve with the W projection kernel
            invker = n.fromfunction(lambda u,v: self.conv_invker(u,v,avg_w),
                uv.shape)
            if not invker2 is None: invker *= invker2
            self.uv += n.fft.ifft2(n.fft.fft2(uv) * invker)
            for i in range(len(self.bm)):
                self.bm[i] += n.fft.ifft2(n.fft.fft2(bm[i]) * invker)
            if j >= len(w): break
            i = j
    def conv_invker(self, u, v, w):
        """Generates the W projection kernel (a function of u,v) for the
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
        G = G.filled(0)
        # Unscramble difference between fft(fft(G)) and G
        G[1:] = n.flipud(G[1:]).copy()
        G[:,1:] = n.fliplr(G[:,1:]).copy()
        return G / G.size

'''
class ImgBW(ImgW):
    """A subclass of ImgW adding better weighting for the primary beam."""
    def __init__(self, size=100, res=1, wres=.5, tres=.2):
        """wres: the gridding resolution of sqrt(w) when projecting to w=0."""
        ImgW.__init__(self, size=size, res=res, wres=wres)
        self.tres = tres
    def put(self, uvw, data, wgts, has, odec, sdec, bm_response):
        """Same as Img.put, only now the w component is projected to the w=0
        plane before applying the data to the UV matrix."""
        if len(uvw) == 0: return
        order = n.argsort(has)
        # Sort uvw in order of time
        uvw = uvw.take(order, axis=0)
        data = data.take(order)
        wgts = wgts.take(order)
        i = 0
        while True:
            # Grab a chunk of uvw's that have same time
            j = has.searchsorted(has[i]+self.tres)
            ha = n.average(has[i:j])
            eqs = self.get_eq(ra=-ha, dec=sdec)
            mask = eqs[0].mask; eqs = eqs.filled(0)
            eqs.shape = (3, eqs.size / 3)
            m = coord.eq2top_m(0, odec)
            top = n.dot(m, eqs)
            z = top[-1]; z.shape = self.uv.shape
            mask = n.logical_or(mask, z < 0)
            resp = bm_response(top); resp.shape = self.uv.shape
            resp = n.where(mask, 0, resp)
            if resp[0,0] <= 0:
                if j >= len(has): break
                i = j
                continue
            resp = resp[0,0] / resp.clip(min(resp[0,0], resp.max()/4),n.Inf)
            #resp = n.where(mask, 0, resp)
            # Transpose here to match transposition in image() and bm_image()
            resp = resp.transpose()
            # Unscramble difference between fft(fft(resp)) and resp
            resp[1:] = n.flipud(resp[1:]).copy()
            resp[:,1:] = n.fliplr(resp[:,1:]).copy()
            resp /= resp.size
            # Put all uv's down on plane for this time
            ImgW.put(self, uvw[i:j,:],data[i:j],wgts[i:j],invker2=resp)
            if j >= len(has): break
            i = j
'''
