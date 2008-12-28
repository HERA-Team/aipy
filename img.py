import numpy, random
from matplotlib.toolkits.basemap import Basemap

#def ang_size(uv_res): return numpy.arcsin(1/uv_res)
def radians(deg): return deg * numpy.pi / 180.
def degrees(rad): return rad * 180. / numpy.pi

#def get_img_coords(uv_size, uv_res, c=(0,0)):
#    uv_res, uv_size = float(uv_res), float(uv_size)
#    v = numpy.arange(-.5/uv_res, .5/uv_res, 1/uv_size)
#    v = numpy.arcsin(v)
#    return c[0] + v, c[1] + v
    
def recenter(a, c):
    """Slide the (0,0) point of matrix a to a new location tuple c.  This is
    useful for making an image centered on your screen after performing an
    inverse fft of uv data."""
    s = a.shape
    c = (c[0] % s[0], c[1] % s[1])
    na = numpy.zeros_like(a)
    na[:s[0] - c[0]], na[s[0] - c[0]:] = a[c[0]:], a[:c[0]]
    nna = numpy.zeros_like(a)
    nna[:,:s[1] - c[1]], nna[:,s[1] - c[1]:] = na[:,c[1]:], na[:,:c[1]]
    return nna

def choose_coords(a, num=1):
    p = numpy.cumsum(a.flat / a.sum())
    r = numpy.random.random((num,))
    s = numpy.searchsorted(p, r)
    x, y = a.shape
    return zip(s / y, s % y)

def convolve2d(a, b):
    return numpy.fft.ifft2(numpy.fft.fft2(a) * numpy.fft.fft2(b))

def gaussian_beam(sigma, shape=0, amp=1., center=(0,0)):
    """Return a 2D gaussian.  Normalized to area under curve = 'amp'.  
    Down by 1/e at distance 'sigma' from 'center'."""
    if type(shape) == type(0): shape = array([2, 2]) * sigma
    def gaussian(x, y):
        nx = numpy.where(x > shape[0] / 2, x - shape[0], x)
        ny = numpy.where(y > shape[1] / 2, y - shape[1], y)
        return numpy.exp(-(nx**2 + ny**2) / sigma**2)
    g = numpy.fromfunction(gaussian, shape)
    g /= g.sum() / amp
    return recenter(g, center)

def deconvolve(im, ker, mdl=None, iter=400, verbose=False):
    dim = im.shape[0]
    ker_pwr = abs(ker).sum()
    if mdl is None: mdl = numpy.zeros_like(im)
    diff = im - abs(convolve2d(mdl, ker))
    score = numpy.average(diff**2)
    for i in range(iter):
        num = numpy.where(abs(diff).sum() > ker_pwr, 1, 0).sum()
        num = random.randint(0, num)
        xy = choose_coords(abs(diff), num=num)
        n_mdl = mdl.copy()
        for a, b in xy:
            n_mdl[a,b] += diff[a,b] * .1 / ker_pwr
        n_mdl = numpy.where(n_mdl < 0, 0, n_mdl)
        n_diff = im - abs(convolve2d(n_mdl, ker))
        n_sqr = n_diff**2
        n_score = numpy.average(n_sqr)
        if verbose: print i, score, n_score
        if n_score < score:
            mdl, diff, sqr, score  = n_mdl, n_diff, n_sqr, n_score
    return mdl


class Img:
    def __init__(self, size=100, res=1):
        self.res = float(res)
        self.size = float(size)
        dim = numpy.round(self.size / self.res)
        self.uv = numpy.zeros(shape=(dim,dim), dtype=numpy.complex64)
        self.bm = numpy.zeros(shape=(dim,dim), dtype=numpy.float32)
    def put(self, uvw, data, wgts=None, apply=True):
        """Assumes hermitian data is in there already."""
        if apply:
            uv = self.uv
            bm = self.bm
        else:
            uv = numpy.zeros_like(self.uv)
            bm = numpy.zeros_like(self.uv)
        if wgts is None: wgts = numpy.ones_like(data)
        inds = numpy.round(uvw[:,:2] / self.res).astype(numpy.int)
        for i, d, w in zip(inds, data, wgts):
            i = tuple(i)
            try:
                uv[i] += d
                #bm[i] += w
                bm[i] += abs(d)
            except(IndexError): pass
        if not apply: return uv, bm
    def append_hermitian(self, uvw, data, wgts=None):
        uvw = numpy.concatenate([uvw, -uvw], axis=0)
        data = numpy.concatenate([data, numpy.conj(data)], axis=0)
        if wgts is None: return uv, data
        wgts = numpy.concatenate([wgts, wgts], axis=0)
        return uvw, data, wgts
    def image(self, center=(0,0)):
        im_img = recenter(numpy.abs(numpy.fft.ifft2(self.uv)), center)
        return im_img
    def bm_image(self, center=(0,0)):
        bm_img = recenter(numpy.abs(numpy.fft.ifft2(self.bm)), center)
        return bm_img
    def __str__(self):
        return str(self.uv)
    def get_coords(self, ra=0, dec=0, fmt='rad'):
        # Create a map to convert the orthographic sky projection
        # into lat, long coordinates
        m = Basemap(projection='ortho',
            lon_0=degrees(ra-numpy.pi), lat_0=degrees(dec), 
            rsphere=1, resolution='c')
        # Convert img coordiates to points on the map (0 to 2)
        x, y = numpy.indices(self.uv.shape)
        x, y = x - self.uv.shape[0]/2., y - self.uv.shape[1]/2.
        x, y = x / self.size + 1, y / self.size + 1
        # Use the orthographic coordinates to invert x,y to long, lat
        x, y = m(x, y, inverse=True)
        # Convert long, lat into ra, dec
        x = x + 180.
        if fmt == 'deg': return x, y
        else: return radians(x), radians(y)


class ImgW(Img):
    def __init__(self, size=100, res=1, wres=.5):
        Img.__init__(self, size=size, res=res)
        self.wres = wres
    def put(self, uvw, data, wgts=None):
        """Assumes hermitian data is in there already."""
        if len(uvw) == 0: return
        if wgts is None: wgts = numpy.ones_like(data)
        # Sort uvw in order of w
        order = numpy.argsort(uvw[:,-1])
        uvw = uvw.take(order, axis=0)
        data = data.take(order)
        w = uvw[:,-1]
        sqrt_w = numpy.sqrt(abs(w)) * numpy.sign(w)
        i = 0
        while True:
            j = sqrt_w.searchsorted(sqrt_w[i]+self.wres)
            avg_w = numpy.average(w[i:j])
            uv, bm = Img.put(self, uvw[i:j,:],data[i:j],wgts[i:j],apply=False)
            ker = numpy.fromfunction(lambda u,v: self.conv_ker(u,v,avg_w),
                uv.shape)
            self.uv += convolve2d(uv, ker)
            self.bm += convolve2d(bm, ker)
            if j >= len(w): break
            i = j
    def conv_ker(self, u, v, w):
        c = self.uv.shape[0] / 2
        u = numpy.where(u > c, u-2*c, u)
        v = numpy.where(v > c, v-2*c, v)
        l, m = numpy.sin(u/self.size), numpy.sin(v/self.size)
        # This is the exactly evaluated kernel (works better)
        valid = l**2 + m**2 <= 1
        sqrt = numpy.sqrt(numpy.where(valid, 1 - l**2 - m**2, 0))
        G = numpy.where(valid, numpy.exp(-2*numpy.pi*1j*w*(sqrt - 1)), 0)
        # This is the kernel described by Cornwell using the small angle approx.
        #G = numpy.exp(numpy.pi*1j*w*(l**2 + m**2))
        G_hat = numpy.fft.fft2(G)
        return G_hat

class SkyMap:
    def __init__(self, res=.01, fromfile=None):
        if not fromfile is None: self.fromfile(fromfile)
        else:
            self.res = float(res)
            d = int(numpy.pi / self.res)
            self.map = numpy.zeros((2*d, d), dtype=numpy.float)
            self.wgt = numpy.zeros((2*d, d), dtype=numpy.float)
    def add_data(self, ras, decs, data, weights=None):
        if weights is None: weights = numpy.ones_like(data)
        ras = ras % (2*numpy.pi)
        decs = (decs + numpy.pi/2) % numpy.pi
        ras = numpy.round(ras / self.res)
        decs = numpy.round(decs / self.res)
        data *= weights
        for r,d,dat,wgt in zip(ras.flatten(), decs.flatten(), 
                data.flatten(), weights.flatten()):
            try:
                self.map[r,d] += dat
                self.wgt[r,d] += wgt
            except(IndexError): pass
    def fromfile(self, filename):
        f = open(filename)
        data = numpy.fromfile(f)
        self.map = data[:data.size/2]
        self.wgt = data[data.size/2:]
        f.close()
        d = int(numpy.round(numpy.sqrt(self.map.size/2)))
        self.map.shape = (2*d, d)
        self.wgt.shape = (2*d, d)
        self.res = numpy.pi / d
    def tofile(self, filename):
        f = open(filename, 'w')
        self.map.tofile(f)
        self.wgt.tofile(f)
        f.close()
    def get_map(self):
        w = numpy.where(self.wgt > 0, self.wgt, 1)
        return self.map / w
