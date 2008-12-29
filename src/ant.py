"""
Module for representing antenna array geometry and for generating
phasing information.
"""
import ephem, math, numpy as n, coord, const, _cephes

class PointingError(Exception):
    """An error to throw if a source is below the horizon."""
    def __init__(self, value): self.parameter = value
    def __str__(self): return repr(self.parameter)

#  _   _ _   _ _ _ _           _____                 _   _                 
# | | | | |_(_) (_) |_ _   _  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | | | | __| | | | __| | | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |_| | |_| | | | |_| |_| | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  \___/ \__|_|_|_|\__|\__, | |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                      |___/        

def juldate2ephem(num):
    """Convert Julian date to ephem date, measured from noon, Jan. 1, 1900."""
    return ephem.date(num - 2415020)

def ephem2juldate(num):
    """Convert ephem date (measured from noon, Jan. 1, 1900) to Julian date."""
    return float(num + 2415020)

#  ____           _ _       ____            _       
# |  _ \ __ _  __| (_) ___ | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \|  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|____/ \___/ \__,_|\__, |
#                                             |___/ 

class RadioBody:
    """The location of a celestial source."""
    def __init__(self, name, mfreq, ionref, srcshape):
        self.src_name = name
        self.mfreq = mfreq
        self.ionref = list(ionref)
        self.srcshape = list(srcshape)
    def compute(self, observer):
        """Update coordinates relative to the provided observer.  Must be
        called at each time step before accessing information."""
        # Generate a map for projecting baselines to uvw coordinates
        self.map = coord.eq2top_m(observer.sidereal_time()-self.ra, self.dec)
    def get_crd(self, crdsys, ncrd=3):
        """Return the coordinates of this location in the desired coordinate
        system ('eq','top') in the current epoch.  If ncrd=2, angular
        coordinates (ra/dec or az/alt) are returned, and if ncrd=3,
        xyz coordinates are returned."""
        assert(crdsys in ('eq','top'))
        assert(ncrd in (2,3))
        if crdsys == 'eq':
            if ncrd == 2: return (self.ra, self.dec)
            return coord.radec2eq((self.ra, self.dec))
        else:
            if ncrd == 2: return (self.az, self.alt)
            return coord.azalt2top((self.az, self.alt))

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(ephem.FixedBody, RadioBody):
    """A source at fixed RA,DEC.  Combines ephem.FixedBody with RadioBody."""
    def __init__(self, ra, dec, mfreq=.150, name='', epoch=ephem.J2000,
            ionref=(0.,0.), srcshape=(0.,0.,0.), **kwargs):
        RadioBody.__init__(self, name, mfreq, ionref, srcshape)
        ephem.FixedBody.__init__(self)
        self._ra, self._dec = ra, dec
        self._epoch = epoch
    def compute(self, observer):
        ephem.FixedBody.compute(self, observer)
        RadioBody.compute(self, observer)

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(RadioBody, object):
    """A moving source (Sun,Moon,planets).  Combines ephem versions of these 
    objects with RadioBody."""
    def __init__(self,name, mfreq=.150,
            ionref=(0.,0.), srcshape=(0.,0.,0.), **kwargs):
        """`name' is used to lookup appropriate ephem celestial object."""
        RadioBody.__init__(self, name, mfreq, ionref, srcshape)
        self.Body = eval('ephem.%s()' % name)
    def __getattr__(self, nm):
        """First try to access attribute from this class, but if that fails, 
        try to get it from the underlying ephem object."""
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Body.__getattribute__(nm)
    def __setattr__(self, nm, val):
        """First try to set attribute for this class, buf if that fails, 
        try to set it for the underlying ephem object."""
        try: object.__setattr__(self, nm, val)
        except(AttributeError): return setattr(self.Body, nm, val)
    def compute(self, observer):
        self.Body.compute(observer)
        RadioBody.compute(self, observer)

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/ 

class SrcCatalog(dict):
    """A catalog of celestial sources."""
    def __init__(self, srcs, **kwargs):
        dict.__init__(self)
        for s in srcs: self.add_src(s)
    def add_src(self, src):
        """Add provided src object (RadioFixedBody,RadioSpecial) to catalog."""
        self[src.src_name] = src
    def get_srcs(self, *args):
        """Return list of all src objects in catalog."""
        return [self[s] for s in args]
    def compute(self, observer):
        """Call compute method of all objects in catalog."""
        for s in self: self[s].compute(observer)
    def get_crds(self, crdsys, ncrd=3, srcs=None):
        """Return coordinates of all objects in catalog."""
        if srcs is None: srcs = self.keys()
        crds = n.array([self[s].get_crd(crdsys, ncrd=ncrd) for s in srcs])
        return crds.transpose()
    def get_mfreqs(self, srcs=None):
        """Return list of frequencies where attributes are measured for all 
        objects in catalog."""
        if srcs is None: srcs = self.keys()
        return n.array([self[s].mfreq for s in srcs])
    def get_ionrefs(self, srcs=None):
        if srcs is None: srcs = self.keys()
        ionref = n.array([self[s].ionref for s in srcs])
        return ionref.transpose()
    def get_srcshapes(self, srcs=None):
        """Return list of angular sizes of all src objects in catalog."""
        if srcs is None: srcs = self.keys()
        a1 = n.array([self[s].srcshape[0] for s in srcs])
        a2 = n.array([self[s].srcshape[1] for s in srcs])
        th = n.array([self[s].srcshape[2] for s in srcs])
        return (a1,a2,th)

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam:
    """Template for representing antenna beam pattern.  Beams also hold
    info about which frequencies are active (i.e. Antennas and AntennaArrays
    access frequencies through Beam)."""
    def __init__(self, freqs, active_chans=None, **kwargs):
        """freqs = frequencies (in GHz) at bin centers across spectrum
        active_chans = indices of channels to use for future calculations."""
        self.freqs = freqs
        self.select_chans(active_chans)
    def select_chans(self, active_chans=None):
        """Select only enumerated channels to use for future calculations."""
        if active_chans is None: active_chans = n.arange(self.freqs.size)
        self.chans = active_chans
        self.afreqs = self.freqs.take(active_chans)

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna:
    """Representation of physical attributes of individual antenna."""
    def __init__(self, x, y, z, beam, delay=0., offset=0., **kwargs):
        """x,y,z = antenna coordinates in equatorial (ns) coordinates
        beam = Beam object
        delay = signal delay (linear phase vs. frequency)
        offset = signal phase offset (constant phase vs. frequency)"""
        self.pos = n.array((x,y,z), n.float64) # must be float64 for mir
        self.beam = beam
        self.delay = delay
        self._offset = offset
        try:
            len(offset)
            self.offset = n.polyval(offset, self.beam.afreqs)
        except(AttributeError,TypeError): self.offset = (offset % 1)
    def select_chans(self, active_chans=None):
        """Select only the specified channels for use in future calculations."""
        self.beam.select_chans(active_chans)
        try:
            len(self._offset)
            self.offset = n.polyval(self._offset, self.beam.afreqs)
        except(AttributeError,TypeError): self.offset = (self._offset % 1)
    def __tuple__(self): return (self.pos[0], self.pos[1], self.pos[2])
    def __list__(self): return [self.pos[0], self.pos[1], self.pos[2]]
    def __add__(self, a): return self.pos + a.pos
    __radd__ = __add__
    def __neg__(self): return -self.pos
    def __sub__(self, a): return self.pos - a.pos
    def __rsub__(self, a): return a.pos - self.pos

#     _                         _                    _   _             
#    / \   _ __ _ __ __ _ _   _| |    ___   ___ __ _| |_(_) ___  _ __  
#   / _ \ | '__| '__/ _` | | | | |   / _ \ / __/ _` | __| |/ _ \| '_ \ 
#  / ___ \| |  | | | (_| | |_| | |__| (_) | (_| (_| | |_| | (_) | | | |
# /_/   \_\_|  |_|  \__,_|\__, |_____\___/ \___\__,_|\__|_|\___/|_| |_|
#                         |___/                                        

class ArrayLocation(ephem.Observer):
    """The location and time of an observation."""
    def __init__(self, location):
        """location = (lat,long,[elev]) of array"""
        ephem.Observer.__init__(self)
        self.pressure = 0
        self.update_location(location)
    def update_location(self, location):
        """Initialize antenna array wth provided location.  May be (lat,long) 
        or (lat,long,elev)."""
        if len(location) == 2: self.lat, self.long = location
        else: self.lat, self.long, self.elev = location
    def get_jultime(self):
        """Get current time as a Julian date."""
        return ephem2juldate(self.date)
    def set_jultime(self, t=None):
        """Set current time as a Julian date."""
        if t is None: t = ephem.julian_date()
        self.set_ephemtime(juldate2ephem(t))
    def set_ephemtime(self, t=None):
        """Set current time as derived from the ephem package."""
        if t is None: t = ephem.now()
        self.date, self.epoch = t, t

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(ArrayLocation):
    """A collection of antennas, their spacings, and location/time of 
    observations."""
    def __init__(self, location, ants, **kwargs):
        """ location = (lat,long,[elev]) of array
        ants = list of Antenna objects."""
        ArrayLocation.__init__(self, location=location)
        self.update_antennas(ants)
        self.select_chans()
    def update_antennas(self, ants):
        """Update antenna array to use the provided list of Antenna objects."""
        self.ants = ants
        self.update()
    def update(self):
        """Update variables derived from antenna parameters/active channels."""
        self._eq2zen = coord.eq2top_m(0., self.lat)
    def set_ephemtime(self, t=None):
        """Set current time as derived from the ephem package.  Recalculates
        matrix for projecting baselines into current positions."""
        ArrayLocation.set_ephemtime(self, t=t)
        self._eq2now = coord.rot_m(-self.sidereal_time(), n.array([0.,0.,1.]))
    def select_chans(self, active_chans=None):
        """Select which channels are used in computations.  Default is all."""
        for a in self.ants: a.select_chans(active_chans)
        self.update()
    def ij2bl(self, i, j):
        """Convert baseline i,j (0 indexed) to Miriad's (i+1) << 8 | (j+1) 
        indexing scheme."""
        return (int(i)+1) << 8 | (int(j)+1)
    def bl2ij(self, bl):
        """Convert Miriad's (i+1) << 8 | (j+1) baseline indexing scheme to 
        i,j (0 indexed)"""
        bl = int(bl)
        return ((bl >> 8) & 255) - 1, (bl & 255) - 1
    def get_baseline(self, i, j, src='z'):
        """Return the baseline corresponding to i,j in various coordinate 
        projections: src='e' for current equatorial, 'z' for zenith 
        topocentric, 'r' for unrotated equatorial, or a RadioBody for
        projection toward that source."""
        bl = self.ants[j] - self.ants[i]
        if type(src) == str:
            if src == 'e': return n.dot(self._eq2now, bl)
            elif src == 'z': return n.dot(self._eq2zen, bl)
            elif src == 'r': return bl
            else: raise ValueError('Unrecognized source:' + src)
        try:
            if src.alt < 0:
                raise PointingError('%s below horizon' % src.src_name)
            m = src.map
        except(AttributeError):
            ra,dec = coord.eq2radec(src)
            m = coord.eq2top_m(self.sidereal_time() - ra, dec)
        return n.dot(m, bl).transpose()
    def get_delay(self, i, j):
        """Return the delay corresponding to baseline i,j."""
        return self.ants[j].delay - self.ants[i].delay
    def get_offset(self, i, j):
        """Return the delay corresponding to baseline i,j."""
        return self.ants[j].offset - self.ants[i].offset
    def gen_uvw(self, i, j, src='z'):
        """Compute uvw coordinates of baseline relative to provided RadioBody, 
        or 'z' for zenith uvw coordinates."""
        x,y,z = self.get_baseline(i,j, src=src)
        afreqs = self.ants[0].beam.afreqs
        if len(x.shape) == 0: return n.array([x*afreqs, y*afreqs, z*afreqs])
        afreqs = n.reshape(afreqs, (1,afreqs.size))
        x.shape += (1,); y.shape += (1,); z.shape += (1,)
        return n.array([n.dot(x,afreqs), n.dot(y,afreqs), n.dot(z,afreqs)])
    def gen_phs(self, src, i, j, doref=False, dores=False,
            mfreq=.150, ionref=(0.,0.), srcshape=(0.,0.,0.)):
        """Return phasing that is multiplied to data to point to src."""
        try:
            if doref: ionref = src.ionref
            if dores: srcshape = srcshape
        except(AttributeError): pass
        u,v,w = self.gen_uvw(i,j,src=src)
        if doref: dw = self.refract(u, v, ionref=ionref)
        else: dw = 0.
        if dores: res = self.resolve_src(u, v, srcshape=srcshape)
        else: res = 1.
        t = self.get_delay(i,j)
        o = self.get_offset(i,j)
        afreqs = self.ants[0].beam.afreqs
        phs = res * n.exp(-1j*2*n.pi*(w + dw + t*afreqs + o))
        return phs.squeeze()
    def resolve_src(self, u, v, srcshape=(0,0,0)):
        """Adjust amplitudes to reflect resolution effects for a uniform 
        elliptical disk characterized by srcshape:
        srcshape = (a1,a2,th) where a1,a2 are angular sizes along the 
            semimajor, semiminor axes, and th is the angle (in radians) of
            the semimajor axis from E."""
        a1,a2,th = srcshape
        try:
            if len(u.shape) > len(a1.shape): 
                a1.shape += (1,); a2.shape += (1,); th.shape += (1,)
        except(AttributeError): pass
        ru = a1 * (u*n.cos(th) - v*n.sin(th))
        rv = a2 * (u*n.sin(th) + v*n.cos(th))
        x = 2 * n.pi * n.sqrt(ru**2 + rv**2)
        # Use first Bessel function of the first kind (J_1)
        return n.where(x == 0, 1, 2 * _cephes.j1(x)/x).squeeze()
    def refract(self, u, v, mfreq=.150, ionref=(0.,0.)):
        """Calibrate a frequency-dependent source offset by scaling measured
        offsets at a given frequency.  Generates dw, a change in the
        projection of a baseline towards that source, which can be used to
        fix the computed phase of that source.
        ionref = (dra, ddec) where dra, ddec are angle offsets (in radians)
            of the source along ra/dec axes at the specified mfreq.
        u,v = u,v components of baseline, which are used to compute the
            change in w given angle offsets and the small angle approx."""
        dra,ddec = ionref
        f2 = (self.ants[0].beam.afreqs / mfreq)**2
        f2.shape = (1, f2.size)
        try:
            if len(dra.shape) < len(f2.shape):
                dra.shape += (1,); ddec.shape += (1,)
        except(AttributeError): pass
        dw = dra / f2 * u + ddec / f2 * v
        return dw
    def phs2src(self, data, src, i, j, doref=False, dores=False,
            mfreq=.150, ionref=(0.,0.), srcshape=(0.,0.,0.)):
        """Apply phasing to zenith-phased data to point to src."""
        return data * self.gen_phs(src, i, j, doref=doref, dores=dores,
            mfreq=mfreq, ionref=ionref, srcshape=srcshape)
    def unphs2src(self, data, src, i, j, doref=False, dores=False,
            mfreq=.150, ionref=(0.,0.), srcshape=(0.,0.,0.)):
        """Remove phasing from src-phased data to point to zenith."""
        return data * n.conjugate(self.gen_phs(src, i, j, doref=doref, 
            dores=dores, mfreq=mfreq, ionref=ionref, srcshape=srcshape))
