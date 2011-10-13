"""
Module for gridding UVW data (including W projection), forming images,
and combining (mosaicing) images into spherical maps.
"""

import numpy as n, utils, coord, pyfits, time
USEDSP = True
if USEDSP: import _dsp

deg2rad = n.pi / 180.
rad2deg = 180. / n.pi

def word_wrap(string, width=80,ind1=0,ind2=0,prefix=''):
    """ word wrapping function.
        string: the string to wrap
        width: the column number to wrap at
        prefix: prefix each line with this string (goes before any indentation)
        ind1: number of characters to indent the first line
        ind2: number of characters to indent the rest of the lines
    """
    awidth = min(width-2-len(prefix+ind1*' '),width-2-len(prefix+ind2*' '))
    words = string.split(' ')
    okwords = []
    chunk = lambda v, l: [v[i*l:(i+1)*l] for i in range(int(n.ceil(len(v)/float(l))))]
    for word in words:
        for okword in chunk(word,awidth):
            okwords.append(okword)
    lines = []
    l = prefix+ind1*' '
    for i,w in enumerate(okwords):
        #print w,len(l+' '+w),width
        if len(l+' ' + w)<width:
            l += ' '+w
        else:
            lines.append(l)
            l = prefix + ind2*' '+w
    lines.append(l)
    return '\n'.join(lines)+'\n'
    
    
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
    if type(shape) == type(0): shape = n.array([2, 2]) * sigma
    def gaussian(x, y):
        nx = n.where(x > shape[0] / 2, x - shape[0], x)
        ny = n.where(y > shape[1] / 2, y - shape[1], y)
        return n.exp(-(nx**2 + ny**2) / sigma**2)
    g = n.fromfunction(gaussian, shape)
    g *= amp
    return recenter(g, center)

def beam_gain(bm):
    #return n.sqrt((n.abs(bm)**2).sum())
    return n.abs(bm).max()

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
        L,M = n.where(L > dim/2, dim-L, -L), n.where(M > dim/2, M-dim, M)
        L,M = L.astype(n.float32)/dim/self.res, M.astype(n.float32)/dim/self.res
        mask = n.where(L**2 + M**2 >= 1, 1, 0)
        L,M = n.ma.array(L, mask=mask), n.ma.array(M, mask=mask)
        return recenter(L, center), recenter(M, center)
    def get_indices(self, u, v):
        """Get the pixel indices corresponding to the provided uv coordinates."""
        if not USEDSP:
            u = n.round(u / self.res).astype(n.int)
            v = n.round(v / self.res).astype(n.int)
            return n.array([-v,u],).transpose()
        else:
            return (-v / self.res).astype(n.float32), (u / self.res).astype(n.float32)
    def get_uv(self):
        """Return the u,v indices of the pixels in the uv matrix."""
        u,v = n.indices(self.shape)
        u = n.where(u < self.shape[0]/2, u, u - self.shape[0])
        v = n.where(v < self.shape[1]/2, v, v - self.shape[1])
        return u*self.res, v*self.res
    def put(self, (u,v,w), data, wgts=None, apply=True):
        """Grid uv data (w is ignored) onto a UV plane.  Data should already
        have the phase due to w removed.  Assumes the Hermitian conjugate
        data is in uvw already (i.e. the conjugate points are not placed for
        you).  If wgts are not supplied, default is 1 (normal weighting).
        If apply is false, returns uv and bm data without applying it do
        the internally stored matrices."""
        if wgts is None:
            wgts = []
            for i in range(len(self.bm)):
                if i == 0: wgts.append(n.ones_like(data))
                else: wgts.append(n.zeros_like(data))
        if len(self.bm) == 1 and len(wgts) != 1: wgts = [wgts]
        assert(len(wgts) == len(self.bm))
        if apply: uv,bm = self.uv,self.bm
        else:
            uv = n.zeros_like(self.uv)
            bm = [n.zeros_like(i) for i in self.bm]
        if not USEDSP:
            inds = self.get_indices(u,v)
            
            ok = n.logical_and(n.abs(inds[:,0]) < self.shape[0],
                n.abs(inds[:,1]) < self.shape[1])
            data = data.compress(ok)
            inds = inds.compress(ok, axis=0)
            utils.add2array(uv, inds, data.astype(uv.dtype))
        else:
            u,v = self.get_indices(u,v)
            _dsp.grid2D_c(uv, u, v, data.astype(uv.dtype))
        
        for i,wgt in enumerate(wgts):
            if not USEDSP:
                wgt = wgt.compress(ok)
                utils.add2array(bm[i], inds, wgt.astype(bm[0].dtype))
            else:
                _dsp.grid2D_c(bm[i], u, v, wgt.astype(bm[0].dtype))
        if not apply: return uv, bm
    def get(self, (u,v,w), uv=None, bm=None):
        """Generate data as would be observed at the provided (u,v,w) based on
        this Img's current uv data.  Phase due to 'w' will be applied to data
        before returning."""
        u,v = u.flatten(), v.flatten()
        if uv is None: uv,bm = self.uv, self.bm[0]
        if not USEDSP:
            # Currently: no interpolation
            inds = self.get_indices(u,v)
            data = self.uv[inds[:,0],inds[:,1]] / self.bm[0][inds[:,0],inds[:,1]]
            data = data.squeeze()
        else:
            u,v = self.get_indices(-u,v)
            u,v = -v,u # XXX necessary, but probably because of axis ordering in FITS files...
            uvdat = n.zeros(u.shape, dtype=n.complex64)
            bmdat = n.zeros(u.shape, dtype=n.complex64)
            _dsp.degrid2D_c(uv, u, v, uvdat)
            _dsp.degrid2D_c(bm, u, v, bmdat)
            #data = uvdat.sum() / bmdat.sum()
            data = uvdat / bmdat
        return data
    def append_hermitian(self, (u,v,w), data, wgts=None):
        """Append to (uvw, data, [wgts]) the points (-uvw, conj(data), [wgts]).
        This is standard practice to get a real-valued image."""
        u = n.concatenate([u, -u], axis=0)
        v = n.concatenate([v, -v], axis=0)
        w = n.concatenate([w, -w], axis=0)
        data = n.concatenate([data, n.conj(data)], axis=0)
        if wgts is None: return n.array((u,v,w)), data
        if len(self.bm) == 1 and len(wgts) != 1: wgts = [wgts]
        assert(len(wgts) == len(self.bm))
        for i,wgt in enumerate(wgts): wgts[i] = n.concatenate([wgt,wgt],axis=0)
        return (u,v,w), data, wgts
    def _gen_img(self, data, center=(0,0)):
        """Return the inverse FFT of the provided data, with the 0,0 point 
        moved to 'center'.  Up=North, Right=East."""
        return recenter(n.fft.ifft2(data).real.astype(n.float32), center)
    def image(self, center=(0,0)):
        """Return the inverse FFT of the UV matrix, with the 0,0 point moved
        to 'center'.  Tranposes to put up=North, right=East."""
        return self._gen_img(self.uv, center=center)
    def bm_image(self, center=(0,0), term=None):
        """Return the inverse FFT of the sample weightings (for all mf_order
        terms, or the specified term if supplied), with the 0,0 point
        moved to 'center'.  Tranposes to put up=North, right=East."""
        if not term is None:
            return self._gen_img(self.bm[term], center=center)
        else:
            return [self._gen_img(b, center=center) for b in self.bm]
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
        if len(mask.shape) == 0: mask = n.zeros(x.shape)
        vec = n.array([a.filled().flatten() for a in (x,y,z)])
        m = coord.top2eq_m(-ra, dec)
        vec = n.dot(m, vec)
        vec.shape = (3,) + shape
        return n.ma.array(vec, mask=[mask,mask,mask])

class ImgW(Img):
    """A subclass of Img adding W projection functionality (see Cornwell
    et al. 2005 "Widefield Imaging Problems in Radio Astronomy")."""
    def __init__(self, size=100, res=1, wres=.5, mf_order=0):
        """wres: the gridding resolution of sqrt(w) when projecting to w=0."""
        Img.__init__(self, size=size, res=res, mf_order=mf_order)
        self.wres = wres
        self.wcache = {}
    def put(self, (u,v,w), data, wgts=None, invker2=None):
        """Same as Img.put, only now the w component is projected to the w=0
        plane before applying the data to the UV matrix."""
        if len(u) == 0: return
        if wgts is None:
            wgts = []
            for i in range(len(self.bm)):
                if i == 0: wgts.append(n.ones_like(data))
                else: wgts.append(n.zeros_like(data))
        if len(self.bm) == 1 and len(wgts) != 1: wgts = [wgts]
        assert(len(wgts) == len(self.bm))
        # Sort uvw in order of w
        order = n.argsort(w)
        u = u.take(order)
        v = v.take(order)
        w = w.take(order)
        data = data.take(order)
        wgts = [wgt.take(order) for wgt in wgts]
        sqrt_w = n.sqrt(n.abs(w)) * n.sign(w)
        i = 0
        while True:
            # Grab a chunk of uvw's that grid w to same point.
            j = sqrt_w.searchsorted(sqrt_w[i]+self.wres)
            print j, len(w)
            avg_w = n.average(w[i:j])
            # Put all uv's down on plane for this gridded w point
            wgtsij = [wgt[i:j] for wgt in wgts]
            uv,bm = Img.put(self, (u[i:j],v[i:j],w[i:j]),
                data[i:j], wgtsij, apply=False)
            # Convolve with the W projection kernel
            invker = n.fromfunction(lambda u,v: self.conv_invker(u,v,avg_w),
                uv.shape)
            if not invker2 is None: invker *= invker2
            self.uv += n.fft.ifft2(n.fft.fft2(uv) * invker)
            for b in range(len(self.bm)):
                self.bm[b] += n.fft.ifft2(n.fft.fft2(bm[b]) * invker)
            if j >= len(w): break
            i = j
    def get(self, (u,v,w)):
        order = n.argsort(w.flat)
        u_,v_,w_ = u.take(order).squeeze(), v.take(order).squeeze(), w.take(order).squeeze()
        sqrt_w = n.sqrt(n.abs(w_)) * n.sign(w_)
        i, d_ = 0, []
        while True:
            # Grab a chunk of uvw's that grid w to same point.
            j = sqrt_w.searchsorted(sqrt_w[i]+self.wres)
            #print j, len(sqrt_w)
            id = n.round(n.average(sqrt_w[i:j]) / self.wres) * self.wres
            if not self.wcache.has_key(id):
                avg_w = n.average(w_[i:j])
                print 'Caching W plane ID=', id
                projker = n.fromfunction(lambda us,vs: self.conv_invker(us,vs,-avg_w), 
                    self.uv.shape).astype(n.complex64)
                uv_wproj = n.fft.ifft2(n.fft.fft2(self.uv) * projker).astype(n.complex64)
                bm_wproj = n.fft.ifft2(n.fft.fft2(self.bm[0]) * projker).astype(n.complex64) # is this right to convolve?
                self.wcache[id] = (uv_wproj, bm_wproj)
                print '%d W planes cached' % (len(self.wcache))
            # Put all uv's down on plane for this gridded w point
            uv_wproj, bm_wproj = self.wcache[id]
            # Could think about improving this by interpolating between w planes.
            d_.append(Img.get(self, (u_[i:j],v_[i:j],w_[i:j]), uv_wproj, bm_wproj))
            if j >= len(sqrt_w): break
            i = j
        d_ = n.concatenate(d_)
        # Put back into original order
        deorder = n.argsort(order)
        return d_.take(deorder)
    def conv_invker(self, u, v, w):
        """Generates the W projection kernel (a function of u,v) for the
        supplied value of w.  See Cornwell et al. 2005 "Widefield Imaging
        Problems in Radio Astronomy" for discussion.  This implementation
        uses a numerically evaluated Fresnel kernel, rather than the
        small-angle approximated one given in the literature."""
        L,M = self.get_LM()
        # This is the exactly evaluated kernel (works better)
        sqrt = n.sqrt(1 - L**2 - M**2).astype(n.complex64)
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
default_fits_format_codes = {
    n.bool_:'L', n.uint8:'B', n.int16:'I', n.int32:'J', n.int64:'K',
    n.float32:'E', n.float64:'D', n.complex64:'C', n.complex128:'M'
}

def to_fits(filename, data, clobber=False,
        axes=('ra---sin','dec---sin','freq','stokes'),
        object='', telescope='', instrument='', observer='', origin='AIPY',
        obs_date=time.strftime('%D'), cur_date=time.strftime('%D'), 
        ra=0, dec=0, d_ra=0, d_dec=0, epoch=2000., 
        freq=0, d_freq=0, bscale=0, bzero=0,history=''):
    """Write image data to a FITS file.  Follows convention of VLA image
    headers.  "axes" describes dimensions of "data" provided.  (ra,dec) are
    the degree coordinates of image center in the specified "epoch". 
    (d_ra,d_dec) are approximate pixel-deltas for ra,dec (approximate because 
    if sine projections of these coordinates are specified--e.g. 
    "ra---sin"--then the deltas change away from the image center).  If a 
    "freq" axis is specified, then "freq" is the frequency of the first entry 
    (in Hz), and "d_freq" is the width of the channel.  The rest are pretty 
    self-explanitory/can be used however you want."""
    data.shape = data.shape + (1,) * (len(axes) - len(data.shape))
    phdu = pyfits.PrimaryHDU(data)
    phdu.data = data.transpose()
    phdu.update_header()
    phdu.header.update('OBJECT', object, comment='SOURCE NAME')
    phdu.header.update('TELESCOP', telescope)
    phdu.header.update('INSTRUME', instrument)
    phdu.header.update('OBSERVER', observer)
    phdu.header.update('DATE-OBS', obs_date, 
        comment='OBSERVATION START DATE DD/MM/YY')
    #phdu.header.update('DATE-MAP', '', 
    #    comment='DATE OF LAST PROCESSING DD/MM/YY')
    phdu.header.update('BSCALE ', bscale,
        comment='REAL = FITS_VALUE * BSCALE + BZERO')
    phdu.header.update('BZERO  ', bzero)
    phdu.header.update('BUNIT  ', 'JY/BEAM ', comment='UNITS OF FLUX')
    phdu.header.update('EQUINOX', epoch, comment='EQUINOX OF RA DEC')
    phdu.header.update('DATAMAX', data.max(), comment='MAX PIXEL VALUE')
    phdu.header.update('DATAMIN', data.min(), comment='MIN PIXEL VALUE')
    for i,ax in enumerate(axes):
        if ax.lower().startswith('ra'): val,delta = (ra, d_ra)
        elif ax.lower().startswith('dec'): val,delta = (dec, d_dec)
        elif ax.lower().startswith('freq'): val,delta = (freq, d_freq)
        elif ax.lower().startswith('stokes'): val,delta = (1, 1)
        else: val,delta = (0,0)
        phdu.header.update('CTYPE%d' % (i+1), ax.upper())
        if ax.lower().startswith('ra') or ax.lower().startswith('dec'):
            phdu.header.update('CRPIX%d' % (i+1), 
                    round(phdu.data.shape[-(i+1)]/2.))
        else:
            phdu.header.update('CRPIX%d' % (i+1), phdu.data.shape[-(i+1)])
        phdu.header.update('CRVAL%d' % (i+1), val)
        phdu.header.update('CDELT%d' % (i+1), delta)
        phdu.header.update('CROTA%d' % (i+1), 0)
    if history!='':
        history = history.split("\n")
        for line in history:
            if len(line)>1:
                for subline in word_wrap(line,70,5,10,'#').split("\n"):
                    phdu.header.add_history(subline)
    phdu.header.update('ORIGIN', origin)
    phdu.header.update('DATE', cur_date, comment='FILE WRITTEN ON DD/MM/YY')
    pyfits.writeto(filename, phdu.data, phdu.header, clobber=True)

def from_fits(filename):
    """Read (data,kwds) from a FITS file.  Matches to_fits() above.  Attempts
    to deduce each keyword listed in to_fits() from the FITS header, but is
    accepting of differences.  Returns values in "kwds" dictionary."""
    phdu = pyfits.open(filename)[0]
    data = phdu.data.transpose()
    kwds = {}
    hitems = (('OBJECT','object'), ('TELESCOP','telescope'),
        ('INSTRUME','instrument'), ('OBSERVER','observer'),
        ('DATE-OBS','obs_date'), ('BSCALE','bscale'), ('BZERO','bzero'),
        ('EQUINOX','epoch'), ('ORIGIN','origin'), ('DATE','cur_date'))
    for fitsname,name in hitems:
        try: kwds[name] = phdu.header[fitsname]
        except(KeyError): pass
    axes = []
    for i in range(phdu.header['NAXIS']):
        try:
            ax = phdu.header['CTYPE%d' % (i+1)].lower()
            axes.append(ax)
            val = phdu.header['CRVAL%d' % (i+1)]
            delta = phdu.header['CDELT%d' % (i+1)]
            if ax.startswith('ra'): kwds['ra'],kwds['d_ra'] = (val,delta)
            elif ax.startswith('dec'): kwds['dec'],kwds['d_dec'] = (val,delta)
            elif ax.startswith('freq'): kwds['freq'],kwds['d_freq']=(val,delta)
            elif ax.startswith('stokes'): pass
            else: pass
        except(KeyError): pass
    kwds['axes'] = axes
    return data, kwds
    
def find_axis(phdu,name):
    #find the axis number for RA
    for k in phdu.header.items():
        if k[0].lower().startswith('ctype'):
            if k[1].lower().startswith(name):            
                return int(k[1][5])


def from_fits_to_fits(infile,outfile,data,kwds,history=None):
    """
    Create a fits file in outfile with data using header from infile using
    kwds to override values in the header.
    See img.to_fits for typical header variables.
    """

    phdu = pyfits.open(infile)[0]
    axes = []
    for i in range(1,phdu.header.get('naxis')):
        type = "CTYPE"+str(i)
        axes.append(phdu.header.get(type))
    print axes
    data.shape = data.shape + (1,) * (len(axes) - len(data.shape))
    phdu.data = data.transpose()
    for i,ax in enumerate(axes):
        if ax.lower().startswith('ra'):
             if kwds.has_key('ra'): val=kwds['ra']
             else: val=None
             if kwds.has_key('d_ra'):delta = kwds['d_ra']
             else: delta=None
        elif ax.lower().startswith('dec'):
             if kwds.has_key('dec'): val=kwds['dec']
             else: val=None
             if kwds.has_key('d_dec'):delta = kwds['d_dec']
             else: delta=None
        elif ax.lower().startswith('freq'):
             if kwds.has_key('freq'): val=kwds['freq']
             else: val=None
             if kwds.has_key('d_freq'):delta = kwds['d_freq']
             else: delta=None
        else: val,delta = None,None
#        elif ax.lower().startswith('dec') and kwds.has_key: val,delta = (dec, d_dec)
#        elif ax.lower().startswith('freq'): val,delta = (freq, d_freq)
#        elif ax.lower().startswith('stokes'): val,delta = (1, 1)
#        else: val,delta = (0,0)
        phdu.header.update('CTYPE%d' % (i+1), ax.upper())
        if ax.lower().startswith('ra') or ax.lower().startswith('dec'):
            phdu.header.update('CRPIX%d' % (i+1), 
                    round(phdu.data.shape[-(i+1)]/2.))
        else:
            phdu.header.update('CRPIX%d' % (i+1), phdu.data.shape[-(i+1)])
        if not val is None: phdu.header.update('CRVAL%d' % (i+1), val);
        if not delta is None: phdu.header.update('CDELT%d' % (i+1), delta)
        phdu.header.update('CROTA%d' % (i+1), 0)
    for k,v in kwds.iteritems():
        try:
            phdu.header.update(k,v)
            print "updated %s with %s"%(k,str(v))
        except(ValueError): 
            print "error on %s "%(k,)
            continue
    if history is None:history = "from_fits_to_fits: from %s to %s"%(infile,
                outfile)
    history = history.split("\n")
    for line in history:
        if len(line)>1:
            for subline in word_wrap(line,70,5,10,'#').split("\n"):
                phdu.header.add_history(subline)
    print phdu.header
    pyfits.writeto(outfile, phdu.data, phdu.header, clobber=True)

    
