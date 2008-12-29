"""
Module adding data simulation support to AntennaArrays.  For most classes,
this means adding gain/amplitude information (as a function of frequency).
"""

import ant, numpy as n, ephem, coord, healpix

#  ____           _ _       ____            _       
# |  _ \ __ _  __| (_) ___ | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \|  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|____/ \___/ \__,_|\__, |
#                                             |___/ 

class RadioBody:
    """Class defining flux and spectral index of a celestial source."""
    def __init__(self, janskies, mfreq, index, angsize):
        """janskies = source strength. Can be polynomial in hour angle (-pi,pi)
        mfreq = frequency (in GHz) where strength was measured
        index = index of power-law spectral model of source emission.  Can be
        polynomial in hour angle."""
        try: len(janskies)
        except: janskies = [janskies]
        try: len(index)
        except: index = [index]
        self._janskies = janskies
        self.mfreq = mfreq
        self._index = index
        self.angsize = angsize
    def compute(self, observer):
        """Update fluxes relative to the provided observer.  Must be
        called at each time step before accessing information."""
        self.janskies = n.clip(n.polyval(self._janskies, 
            (observer.sidereal_time()-self.ra+n.pi) % (2*n.pi) - n.pi), 
            0, n.Inf)
        self.index = n.polyval(self._index, observer.sidereal_time())

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(ant.RadioFixedBody, RadioBody):
    """Class representing a source at fixed RA,DEC.  Adds flux information to
    ant.RadioFixedBody."""
    def __init__(self, ra, dec, janskies, mfreq=.150, 
            index=-1., angsize=0., name='', **kwargs):
        """ra = source's right ascension (epoch=J2000)
        dec = source's declination (epoch=J2000)
        janskies = source strength. Can be polynomial in hour angle (-pi,pi)
        mfreq = frequency (in GHz) where strength was measured
        index = index of power-law spectral model of source emission.  Can be
        polynomial in hour angle."""
        ant.RadioFixedBody.__init__(self, ra, dec, name=name)
        RadioBody.__init__(self, janskies, mfreq, index, angsize)
    def compute(self, observer):
        ant.RadioFixedBody.compute(self, observer)
        RadioBody.compute(self, observer)

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(ant.RadioSpecial, RadioBody):
    """Class representing moving sources (Sun,Moon,planets).  Adds flux
    information to ant.RadioSpecial."""
    def __init__(self,name,janskies,mfreq=.150,index=-1.,angsize=0.,**kwargs):
        """janskies = source strength. Can be polynomial in hour angle (-pi,pi)
        mfreq = frequency (in GHz) where strength was measured
        index = index of power-law spectral model of source emission.  Can be
        polynomial in hour angle."""
        ant.RadioSpecial.__init__(self, name)
        RadioBody.__init__(self, janskies, mfreq, index, angsize)
    def compute(self, observer):
        ant.RadioSpecial.compute(self, observer)
        RadioBody.compute(self, observer)

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/ 

class SrcCatalog(ant.SrcCatalog):
    """Class for holding a catalog of celestial sources."""
    def get_fluxes(self, srcs=None):
        """Return list of fluxes of all src objects in catalog."""
        if srcs is None: srcs = self.keys()
        return n.array([self[s].janskies for s in srcs])
    def get_indices(self, srcs=None):
        """Return list of spectral indices of all src objects in catalog."""
        if srcs is None: srcs = self.keys()
        return n.array([self[s].index for s in srcs])
    def get_mfreqs(self, srcs=None):
        """Return list of frequencies where strength is measured for all src 
        objects in catalog."""
        if srcs is None: srcs = self.keys()
        return n.array([self[s].mfreq for s in srcs])
    def get_angsizes(self, srcs=None):
        """Return list of angular sizes of all src objects in catalog."""
        if srcs is None: srcs = self.keys()
        return n.array([self[s].angsize for s in srcs])

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class BeamFlat(ant.Beam):
    """Representation of a flat (gain=1) antenna beam pattern."""
    def response(self, xyz):
        x,y,z = n.array(xyz)
        return n.ones((self.afreqs.size, x.size))

class Beam2DGaussian(ant.Beam):
    """Representation of a 2D Gaussian beam pattern, with default setting for 
    a flat beam."""
    def __init__(self, freqs, xwidth=n.Inf, ywidth=n.Inf):
        """xwidth = angular width (radians) in EW direction
        ywidth = angular width (radians) in NS direction"""
        ant.Beam.__init__(self, freqs)
        self.update(xwidth=xwidth, ywidth=ywidth)
    def update(self, xwidth=None, ywidth=None):
        """Set the width in the x (EW) and y (NS) directions of the gaussian 
        beam."""
        if not xwidth is None: self.xwidth = xwidth
        if not ywidth is None: self.ywidth = ywidth
    def response(self, xyz):
        """Return beam response across active band for specified topocentric 
        coordinates: (x=E,y=N,z=UP). x,y,z may be arrays of multiple 
        coordinates.  Returns 'x' linear polarization (rotate pi/2 for 'y')."""
        x,y,z = xyz
        x,y = n.arcsin(x)/self.xwidth, n.arcsin(y)/self.ywidth
        resp = n.sqrt(n.exp(-(x**2 + y**2)))
        resp = n.resize(resp, (self.afreqs.size, resp.size))
        return resp

class BeamPolynomial(ant.Beam):
    """Representation of a gaussian beam model whose width varies with azimuth
    angle and with frequency."""
    def __init__(self, freqs, poly_azfreq=n.array([[.5]])):
        """poly_azfreq = a 2D polynomial in cos(2*n*az) for first axis and 
        in freq**n for second axis."""
        self.poly = poly_azfreq
        ant.Beam.__init__(self, freqs)
        self.update(poly_azfreq)
    def select_chans(self, active_chans):
        """Select only enumerated channels to use for future calculations."""
        ant.Beam.select_chans(self, active_chans)
        self.update()
    def update(self, poly_azfreq=None):
        """Update beam with new polynomial coefficients.
        poly_azfreq = a 2D polynomial in cos(2*n*az) for first axis and 
        in freq**n for second axis."""
        if poly_azfreq is None: poly_azfreq = self.poly
        elif len(poly_azfreq.shape) == 1: poly_azfreq.shape = self.poly.shape
        self.poly = poly_azfreq
        f = n.resize(self.afreqs, (self.poly.shape[1], self.afreqs.size))
        f = f**n.array([range(self.poly.shape[1])]).transpose()
        self.sigma = n.dot(self.poly, f)
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

class BeamAlm(ant.Beam):
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
        self.alm = [healpix.Alm(lmax,mmax) for i in range(deg+1)]
        self.hmap = [healpix.HealpixMap(nside,scheme='RING',interp=True)
            for a in self.alm]
        ant.Beam.__init__(self, freqs)
        self.update(coeffs)
    def update(self, coeffs={}):
        """Update beam model using new set of coefficients.
        coeffs = dictionary of polynomial term (integer) and corresponding Alm 
        coefficients (see healpix.py doc)."""
        for c in coeffs:
            if c >= len(self.alm): continue
            self.alm[-1-c].set_data(coeffs[c])
            self.hmap[-1-c].from_alm(self.alm[-1-c])
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

class Antenna(ant.Antenna):
    """Representation of physical location and beam pattern of individual 
    antenna in array."""
    def __init__(self, x, y, z, beam, delay=0., offset=0., bp_r=n.array([1]),
            bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0), **kwargs):
        """x,y z = antenna coordinates in equatorial (ns) coordinates
        beam = Beam object (implements response() function)
        delay = signal delay (linear phase vs. frequency)
        offset = signal phase offset (constant phase vs. frequency)
        bp_r = polynomial (in freq) modeling real component of passband
        bp_i = polynomial (in freq) modeling imaginary component of passband
        amp = overall multiplicative scaling of gain
        pointing = antenna pointing (az,alt).  Default is zenith."""
        ant.Antenna.__init__(self, x,y,z, beam=beam, delay=delay, offset=offset)
        self.update_gain(bp_r, bp_i, amp)
        self.update_pointing(*pointing)
    def select_chans(self, active_chans=None):
        """Select only enumerated channels to use for future calculations."""
        ant.Antenna.select_chans(self, active_chans)
        self.update_gain()
    def update_gain(self, bp_r=None, bp_i=None, amp=None):
        """Update passband information using provided parameters.
        bp_r = polynomial (in freq) modeling real component of passband
        bp_i = polynomial (in freq) modeling imaginary component of passband
        amp = overall multiplicative scaling of gain"""
        if not bp_r is None:
            try: len(bp_r)
            except: bp_r = [bp_r]
            self.bp_r = bp_r
        if not bp_i is None:
            try: len(bp_i)
            except: bp_i = [bp_i]
            self.bp_i = bp_i
        if not amp is None: self.amp = amp
        bp = n.polyval(self.bp_r, self.beam.afreqs) + \
             1j*n.polyval(self.bp_i, self.beam.afreqs)
        self._gain = self.amp * bp
    def update_pointing(self, az=0, alt=n.pi/2, twist=0):
        """Set the antenna beam to point at (az, alt) with specified
        right-hand twist to polarizations.  Polarization y is assumed
        to be +pi/2 azimuth from pol x."""
        y, z = n.array([0,1,0]), n.array([0,0,1])
        rot = coord.rot_m(twist, z)
        rot = n.dot(rot, coord.rot_m(alt-n.pi/2, y))
        rot = n.dot(rot, coord.rot_m(-az, z))
        self.rot_pol_x = rot
        self.rot_pol_y = n.dot(coord.rot_m(-n.pi/2, z), rot)
    def bm_response(self, top, pol='x'):
        """Return response of beam for specified polarization."""
        top = n.array(top)
        top = {'x':top, 'y':n.dot(self.rot_pol_y, top)}[pol]
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

class AntennaArray(ant.AntennaArray):
    """Representation of location and time of observation, and response of
    array of antennas as function of pointing and frequency."""
    def set_jultime(self, t=None):
        """Set current time as a Julian date."""
        ant.AntennaArray.set_jultime(self, t=t)
        self.eq2top_m = coord.eq2top_m(-self.sidereal_time(), self.lat)
        self._cache = None
    def passband(self, i, j):
        return self.ants[i].passband() * self.ants[j].passband(conj=True)
    def bm_response(self, i, j, pol='xx'):
        assert(pol in ('xx','yy','xy','yx'))
        p1, p2 = pol
        # Check that we have cached results needed.  If not, cache them.
        for c,p in zip([i,j], [p1,p2]):
            if not self._cache.has_key(c): self._cache[c] = {}
            if not self._cache[c].has_key(p):
                x,y,z = self._cache['s_top']
                resp = self.ants[c].bm_response((x,y,z), pol=p).transpose()
                self._cache[c][p] = resp
        return self._cache[i][p1] * n.conjugate(self._cache[j][p2])
    def sim_cache(self, s_eqs, fluxes=n.array([1.]), indices=n.array([0.]), 
            mfreqs=n.array([.150]), angsizes=None):
        """Cache intermediate computations given catalog information to speed
        simulation for multiple baselines.  For efficiency, should only be 
        called once per time setting.  MUST be called before sim().
        s_eqs = array of equatorial vectors for all celestial sources
        fluxes = array of fluxes for all celestial sources
        indices = array of spectral indices of all celestial sources
        mfreqs = array of frequencies where fluxes were measured
        angsizes = array of angular sizes of all celestial sources."""
        # Get topocentric coordinates of all srcs
        src_top = n.dot(self.eq2top_m, s_eqs)
        # Throw out everything that is below the horizon
        valid = n.logical_and(src_top[2,:] > 0, fluxes > 0)
        if n.all(valid == 0): self._cache = {}
        else:
            fluxes = fluxes.compress(valid)
            indices = indices.compress(valid)
            mfreqs = mfreqs.compress(valid)
            if not angsizes is None: angsizes = angsizes.compress(valid)
            src_top = src_top.compress(valid, axis=1)
            s_eqs = s_eqs.compress(valid, axis=1)
            # Get src fluxes vs. freq
            fluxes.shape = (fluxes.size, 1)
            mfreqs.shape = (mfreqs.size, 1)
            indices.shape = (indices.size, 1)
            freqs = n.resize(self.ants[0].beam.afreqs, 
                (fluxes.size, self.ants[0].beam.afreqs.size))
            self._cache = {
                'I_sf':fluxes * (freqs / mfreqs)**indices,
                's_eqs': s_eqs.transpose(),
                's_top': src_top,
                's_sz': angsizes
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
        I_sf = self._cache['I_sf'] * \
            self.resolve_src(s_eqs, i, j, self._cache['s_sz'])
        Gij_sf = self.passband(i,j)
        Bij_sf = self.bm_response(i,j,pol=pol)
        if len(Bij_sf.shape) == 2: Gij_sf = n.reshape(Gij_sf, (1, Gij_sf.size))
        # Get the phase of each src vs. freq, also does resolution effects
        E_sf = n.conjugate(self.gen_phs(s_eqs, i, j))
        try: E_sf.shape = I_sf.shape
        except(AttributeError): pass
        # Combine and sum over sources
        GBIE_sf = Gij_sf * Bij_sf * I_sf * E_sf
        Vij_f = GBIE_sf.sum(axis=0)
        return Vij_f
