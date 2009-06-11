"""
Module adding simulation support to RadioBodys and AntennaArrays.
Mostly, this means adding gain/amplitude information vs. frequency.
"""

import phs, numpy as n, ephem, coord, healpix

#  ____           _ _       ____            _       
# |  _ \ __ _  __| (_) ___ | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \|  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|____/ \___/ \__,_|\__, |
#                                             |___/ 

class RadioBody:
    """Class defining flux and spectral index of a celestial source."""
    def __init__(self, jys, index):
        """jys = source strength in Janskies
        mfreq = frequency (in GHz) where strength was measured
        index = index of power-law spectral model of source emission."""
        self._jys = jys
        self.index = index
    def update_jys(self, afreqs):
        """Update fluxes relative to the provided observer.  Must be
        called at each time step before accessing information."""
        self.jys = self._jys * (afreqs / self.mfreq)**self.index
    def get_jys(self):
        """Return the fluxes vs. freq that should be used for simulation."""
        return self.jys

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(phs.RadioFixedBody, RadioBody):
    """Class representing a source at fixed RA,DEC.  Adds flux information to
    phs.RadioFixedBody."""
    def __init__(self, ra, dec, name='', epoch=ephem.J2000,
            jys=0., index=-1, mfreq=.150,
            ionref=(0.,0.), srcshape=(0.,0.,0.), **kwargs):
        """ra = source's right ascension (epoch=J2000)
        dec = source's declination (epoch=J2000)
        jys = source strength in Janskies at mfreq)
        mfreq = frequency (in GHz) where source strength was measured
        index = power-law index of source emission vs. freq."""
        phs.RadioFixedBody.__init__(self, ra, dec, mfreq=mfreq, name=name,
            epoch=epoch, ionref=ionref, srcshape=srcshape)
        RadioBody.__init__(self, jys, index)
    def compute(self, observer):
        phs.RadioFixedBody.compute(self, observer)
        RadioBody.update_jys(self, observer.get_afreqs())

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(phs.RadioSpecial, RadioBody):
    """Class representing moving sources (Sun,Moon,planets).  Adds flux
    information to phs.RadioSpecial."""
    def __init__(self, name, jys=0., index=-1., mfreq=.150,
            ionref=(0.,0.), srcshape=(0.,0.,0.), **kwargs):
        """jys = source strength in Janskies at mfreq)
        mfreq = frequency (in GHz) where source strength was measured
        index = power-law index of source emission vs. freq."""
        phs.RadioSpecial.__init__(self, name, mfreq=mfreq,
            ionref=ionref, srcshape=srcshape)
        RadioBody.__init__(self, jys, index)
    def compute(self, observer):
        phs.RadioSpecial.compute(self, observer)
        RadioBody.update_jys(self, observer.get_afreqs())

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/

class SrcCatalog(phs.SrcCatalog):
    """Class for holding a catalog of celestial sources."""
    def get_jys(self, srcs=None):
        """Return list of fluxes of all src objects in catalog."""
        if srcs is None: srcs = self.keys()
        return n.array([self[s].get_jys() for s in srcs])
    def update_jys(self, afreqs):
        for s in self.keys(): self[s].update_jys(afreqs)

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam(phs.Beam):
    """Representation of a flat (gain=1) antenna beam pattern."""
    def response(self, xyz):
        """Return the (unity) beam response as a function of position."""
        x,y,z = n.array(xyz)
        return n.ones((self.afreqs.size, x.size))

class Beam2DGaussian(phs.Beam):
    """Representation of a 2D Gaussian beam pattern, with default setting for 
    a flat beam."""
    def __init__(self, freqs, xwidth=n.Inf, ywidth=n.Inf):
        """xwidth = angular width (radians) in EW direction
        ywidth = angular width (radians) in NS direction"""
        phs.Beam.__init__(self, freqs)
        self.xwidth, self.ywidth = xwidth, ywidth
    def response(self, xyz):
        """Return beam response across active band for specified topocentric 
        coordinates: (x=E,y=N,z=UP). x,y,z may be arrays of multiple 
        coordinates.  Returns 'x' linear polarization (rotate pi/2 for 'y')."""
        x,y,z = xyz
        x,y = n.arcsin(x)/self.xwidth, n.arcsin(y)/self.ywidth
        resp = n.sqrt(n.exp(-(x**2 + y**2)))
        resp = n.resize(resp, (self.afreqs.size, resp.size))
        return resp

class BeamPolynomial(phs.Beam):
    """Representation of a gaussian beam model whose width varies with azimuth
    angle and with frequency."""
    def __init__(self, freqs, poly_azfreq=n.array([[.5]])):
        """poly_azfreq = a 2D polynomial in cos(2*n*az) for first axis and 
        in freq**n for second axis."""
        self.poly = poly_azfreq
        phs.Beam.__init__(self, freqs)
        self.poly = poly_azfreq
        self._update_sigma()
    def select_chans(self, active_chans):
        """Select only enumerated channels to use for future calculations."""
        phs.Beam.select_chans(self, active_chans)
        self.update()
    def _update_sigma(self):
        f = n.resize(self.afreqs, (self.poly.shape[1], self.afreqs.size))
        f = f**n.array([range(self.poly.shape[1])]).transpose()
        self.sigma = n.dot(self.poly, f)
    def update(self):
        phs.Beam.update(self)
        self._update_sigma()
    def response(self, top):
        """Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y')."""
        az,alt = coord.top2azalt(top)
        zang = n.pi/2 - alt
        if zang.size == 1:
            zang = n.array([zang]); zang.shape = (1,)
            az = n.array([az]); az.shape = (1,)
        a = 2 * n.arange(self.poly.shape[0], dtype=n.float)
        a.shape = (1,) + a.shape; az.shape += (1,); zang.shape += (1,)
        a = n.cos(n.dot(az, a))
        a[:,0] = 0.5
        s = n.dot(a, self.sigma)
        return n.sqrt(n.exp(-(zang/s)**2)).transpose()

class BeamAlm(phs.Beam):
    """Representation of a beam model where each pointing has a response
    defined as a polynomial in frequency, and the spatial distributions of 
    these coefficients decomposed into spherical harmonics."""
    def __init__(self, freqs, lmax=8, mmax=8, deg=7, nside=64, coeffs={}):
        """lmax = maximum spherical harmonic term
        mmax = maximum spherical harmonic term in the z direction
        deg = order of polynomial to used for mapping response of each pointing
        nside = resolution of underlying HealpixMap to use
        coeffs = dictionary of polynomial term (integer) and corresponding Alm 
        coefficients (see healpix.py doc)."""
        phs.Beam.__init__(self, freqs)
        self.alm = [healpix.Alm(lmax,mmax) for i in range(deg+1)]
        self.hmap = [healpix.HealpixMap(nside,scheme='RING',interp=True)
            for a in self.alm]
        for c in coeffs:
            if c < len(self.alm): self.alm[-1-c].set_data(coeffs[c])
        self._update_hmap()
    def _update_hmap(self):
        for c,alm in enumerate(self.alm): self.hmap[c].from_alm(self.alm[c])
    def update(self):
        """Update beam model using new set of coefficients.
        coeffs = dictionary of polynomial term (integer) and corresponding Alm 
        coefficients (see healpix.py doc)."""
        phs.Beam.update(self)
        self._update_hmap()
    def response(self, top):
        """Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y')."""
        top = [healpix.mk_arr(c, dtype=n.double) for c in top]
        px,wgts = self.hmap[0].crd2px(*top, **{'interpolate':1})
        poly = n.array([n.sum(h.map[px] * wgts, axis=-1) for h in self.hmap])
        rv = n.polyval(poly, n.reshape(self.afreqs, (self.afreqs.size, 1)))
        return rv

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna(phs.Antenna):
    """Representation of physical location and beam pattern of individual 
    antenna in array."""
    def __init__(self, x, y, z, beam, phsoff=[0.,0.], bp_r=n.array([1]),
            bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0), **kwargs):
        """x,y z = antenna coordinates in equatorial (ns) coordinates
        beam = Beam object (implements response() function)
        phsoff = polynomial phase vs. frequency.  Phs term that is linear
                 with freq is often called 'delay'.
        bp_r = polynomial (in freq) modeling real component of passband
        bp_i = polynomial (in freq) modeling imaginary component of passband
        amp = overall multiplicative scaling of gain
        pointing = antenna pointing (az,alt).  Default is zenith."""
        phs.Antenna.__init__(self, x,y,z, beam=beam, phsoff=phsoff)
        self.set_pointing(*pointing)
        self.bp_r = bp_r
        self.bp_i = bp_i
        self.amp = amp
        self._update_gain()
    def _update_gain(self):
        bp = n.polyval(self.bp_r, self.beam.afreqs) + \
             1j*n.polyval(self.bp_i, self.beam.afreqs)
        self._gain = self.amp * bp
    def update(self):
        phs.Antenna.update(self)
        self._update_gain()
    def set_pointing(self, az=0, alt=n.pi/2, twist=0):
        """Set the antenna beam to point at (az, alt) with specified
        right-hand twist to polarizations.  Polarization y is assumed
        to be +pi/2 azimuth from pol x."""
        y, z = n.array([0,1,0]), n.array([0,0,1])
        # Twist is negative b/c we apply it to the coords, not the beam
        rot = coord.rot_m(-twist, z)
        rot = n.dot(rot, coord.rot_m(alt-n.pi/2, y))
        rot = n.dot(rot, coord.rot_m(-az, z))
        self.rot_pol_x = rot
        self.rot_pol_y = n.dot(coord.rot_m(-n.pi/2, z), rot)
    def bm_response(self, top, pol='x'):
        """Return response of beam for specified polarization."""
        top = n.array(top)
        top = {'x':n.dot(self.rot_pol_x, top), 
               'y':n.dot(self.rot_pol_y, top)}[pol]
        x,y,z = top
        return self.beam.response((x,y,z))
    def passband(self, conj=False):
        if conj: return n.conjugate(self._gain)
        else: return self._gain

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(phs.AntennaArray):
    """Representation of location and time of observation, and response of
    array of antennas as function of pointing and frequency."""
    def set_jultime(self, t=None):
        """Set current time as a Julian date."""
        phs.AntennaArray.set_jultime(self, t=t)
        self.eq2top_m = coord.eq2top_m(-self.sidereal_time(), self.lat)
        self._cache = None
    def passband(self, i, j):
        """Return the passband response of baseline i,j."""
        return self[i].passband() * self[j].passband(conj=True)
    def bm_response(self, i, j, pol='xx'):
        """Return the beam response towards the cached source positions
        for baseline i,j with the specified polarization."""
        assert(pol in ('xx','yy','xy','yx'))
        p1, p2 = pol
        # Check that we have cached results needed.  If not, cache them.
        for c,p in zip([i,j], [p1,p2]):
            if not self._cache.has_key(c): self._cache[c] = {}
            if not self._cache[c].has_key(p):
                x,y,z = self._cache['s_top']
                resp = self[c].bm_response((x,y,z), pol=p).transpose()
                self._cache[c][p] = resp
        return self._cache[i][p1] * n.conjugate(self._cache[j][p2])
    def sim_cache(self, s_eqs, jys=n.array([1.]), mfreqs=0.150,
            ionrefs=(0.,0.), srcshapes=(0,0,0)):
        """Cache intermediate computations given catalog information to speed
        simulation for multiple baselines.  For efficiency, should only be 
        called once per time setting.  MUST be called before sim().
        s_eqs = array of equatorial vectors for all celestial sources
        jys = array of janskies vs. freq for all celestial sources
        mfreqs = array of frequencies where ionrefs were measured
        ionrefs = (dra,ddec), angular offsets in radians for ra/dec at the
            frequency `mfreq'.
        srcshapes = (a1,a2,th) where a1,a2 are angular sizes along the 
            semimajor, semiminor axes, and th is the angle (in radians) of
            the semimajor axis from E."""
        # Get topocentric coordinates of all srcs
        src_top = n.dot(self.eq2top_m, s_eqs)
        # Throw out everything that is below the horizon
        valid = src_top[2,:] > 0
        if n.all(valid == 0): self._cache = {}
        else:
            jys = jys.compress(valid, axis=0)
            try:
                mfreqs = mfreqs.compress(valid)
                mfreqs.shape = (mfreqs.size, 1)
            except(AttributeError): pass
            a1,a2,th = srcshapes
            try:
                a1 = a1.compress(valid)
                a2 = a2.compress(valid)
                th = th.compress(valid)
            except(AttributeError): pass
            dra,ddec = ionrefs
            try:
                dra = dra.compress(valid)
                ddec = ddec.compress(valid)
            except(AttributeError): pass
            src_top = src_top.compress(valid, axis=1)
            s_eqs = s_eqs.compress(valid, axis=1)
            # Get src fluxes vs. freq
            self._cache = {
                'jys':   jys,
                'mfreq': mfreqs,
                's_eqs': s_eqs,
                's_top': src_top,
                's_shp': (a1,a2,th),
                'i_ref': (dra,ddec),
            }
    def sim(self, i, j, pol='xx'):
        """Simulate visibilites for the specified (i,j) baseline and 
        polarization.  sim_cache() must be called at each time step before 
        this will return valid results."""
        assert(pol in ('xx','yy','xy','yx'))
        if self._cache is None:
            raise RuntimeError('sim_cache() must be called before the first sim() call at each time step.')
        elif self._cache == {}:
            return n.zeros_like(self.passband(i,j))
        s_eqs = self._cache['s_eqs']
        u,v,w = self.gen_uvw(i, j, src=s_eqs)
        I_sf = self._cache['jys']
        Gij_sf = self.passband(i,j)
        Bij_sf = self.bm_response(i,j,pol=pol)
        if len(Bij_sf.shape) == 2: Gij_sf = n.reshape(Gij_sf, (1, Gij_sf.size))
        # Get the phase of each src vs. freq, also does resolution effects
        E_sf = n.conjugate(self.gen_phs(s_eqs, i, j, mfreq=self._cache['mfreq'],
            srcshape=self._cache['s_shp'], ionref=self._cache['i_ref'],
            resolve_src=True))
        try: E_sf.shape = I_sf.shape
        except(AttributeError): pass
        # Combine and sum over sources
        GBIE_sf = Gij_sf * Bij_sf * I_sf * E_sf
        Vij_f = GBIE_sf.sum(axis=0)
        return Vij_f
