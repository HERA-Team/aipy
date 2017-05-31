"""
Module for representing antenna array geometry and for generating
phasing information.
"""
import ephem, math, numpy as np, coord, const, _cephes
from miriad import ij2bl, bl2ij

class PointingError(Exception):
    """An error to throw if a source is below the horizon."""
    def __init__(self, value): self.parameter = value
    def __str__(self): return str(self.parameter)

#  _   _ _   _ _ _ _           _____                 _   _                 
# | | | | |_(_) (_) |_ _   _  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | | | | __| | | | __| | | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |_| | |_| | | | |_| |_| | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  \___/ \__|_|_|_|\__|\__, | |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                      |___/        

def juldate2ephem(num):
    """Convert Julian date to ephem date, measured from noon, Dec. 31, 1899."""
    return ephem.date(num - 2415020.)

def ephem2juldate(num):
    """Convert ephem date (measured from noon, Dec. 31, 1899) to Julian date."""
    return float(num + 2415020.)

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
    def __str__(self):
        return "%s"% (self.src_name,)
    def compute(self, observer):
        """Update coordinates relative to the provided observer.  Must be
        called at each time step before accessing information."""
        # Generate a map for projecting baselines to uvw coordinates
        self.map = coord.eq2top_m(observer.sidereal_time()-self.ra, self.dec)
    def get_crds(self, crdsys, ncrd=3):
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
    def __str__(self):
        if self._dec<0: return RadioBody.__str__(self) + '\t' + str(self._ra) +'\t'+ str(self._dec)
        else: return RadioBody.__str__(self) + '\t' + str(self._ra) +'\t'+'+' + str(self._dec)
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
    """A catalog of celestial sources.  Can be initialized with a list
    of src objects, of as an empty catalog."""
    def __init__(self, *srcs, **kwargs):
        dict.__init__(self)
        self.add_srcs(*srcs)
    def add_srcs(self, *srcs):
        """Add src object(s) (RadioFixedBody,RadioSpecial) to catalog."""
        if len(srcs) == 1 and getattr(srcs[0], 'src_name', None) == None:
            srcs = srcs[0]
        for s in srcs: self[s.src_name] = s
    def get_srcs(self, *srcs):
        """Return list of all src objects in catalog."""
        if len(srcs) == 0: srcs = self.keys()
        elif len(srcs) == 1 and type(srcs[0]) != str:
            return [self[s] for s in srcs[0]]
        else: return [self[s] for s in srcs]
    def compute(self, observer):
        """Call compute method of all objects in catalog."""
        for s in self: self[s].compute(observer)
    def get_crds(self, crdsys, ncrd=3, srcs=None):
        """Return coordinates of all objects in catalog."""
        if srcs is None: srcs = self.keys()
        crds = np.array([self[s].get_crds(crdsys, ncrd=ncrd) for s in srcs])
        return crds.transpose()
    def get(self, attribute, srcs=None):
        """Return the specified source attribute (e.g. "mfreq" for src.mfreq)
        in an array for all src names in 'srcs'.  If not provided, defaults to 
        all srcs in catalog."""
        if srcs is None: srcs = self.keys()
        return np.array([getattr(self[s], attribute) for s in srcs]).transpose()

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam:
    """Template for representing antenna beam pattern.  Beams also hold
    info about which frequencies are active (i.e. Antennas and AntennaArrays
    access frequencies through Beam)."""
    def __init__(self, freqs, **kwargs):
        """freqs = frequencies (in GHz) at bin centers across spectrum."""
        self.freqs = freqs
        self.chans = np.arange(self.freqs.size)
        self._update_afreqs()
    def _update_afreqs(self):
        self.afreqs = self.freqs.take(self.chans)
    def update(self):
        self._update_afreqs()
    def select_chans(self, active_chans=None):
        """Select only enumerated channels to use for future calculations."""
        if active_chans is None: active_chans = np.arange(self.freqs.size)
        self.chans = active_chans
        self.update()

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna:
    """Representation of physical attributes of individual antenna."""
    def __init__(self, x, y, z, beam, phsoff=[0.,0.], **kwargs):
        """x,y,z = antenna coordinates in equatorial (ns) coordinates
        beam = Beam object
        phsoff = polynomial phase vs. frequency.  Phs term that is linear
                 with freq is often called 'delay'."""
        self.pos = np.array((x,y,z), np.float64) # must be float64 for mir
        self.beam = beam
        self._phsoff = phsoff
        self._update_phsoff()
    def select_chans(self, active_chans=None):
        """Select only the specified channels for use in future calculations."""
        self.beam.select_chans(active_chans)
        self.update()
    def _update_phsoff(self):
        self.phsoff = np.polyval(self._phsoff, self.beam.afreqs)
    def update(self):
        self._update_phsoff()
    def __iter__(self): return self.pos.__iter__()
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
        if len(location) == 2: self.lat, self.long = location
        else: self.lat, self.long, self.elev = location
        self._update_eq2zen()
    def _update_eq2zen(self):
        self._eq2zen = coord.eq2top_m(0., self.lat)
    def update(self):
        self._update_eq2zen()
    def get_jultime(self):
        """Get current time as a Julian date."""
        return ephem2juldate(self.date)
    def set_jultime(self, t=None):
        """Set current time as a Julian date."""
        if t is None: t = ephem.julian_date()
        self.set_ephemtime(juldate2ephem(t))
    def set_ephemtime(self, t=None):
        """Set current time as derived from the ephem package.  Recalculates
        matrix for projecting baselines into current positions."""
        if t is None: t = ephem.now()
        self.date, self.epoch = t, t
        self._eq2now = coord.rot_m(-self.sidereal_time(), np.array([0.,0.,1.]))

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
        self.ants = ants
    def __iter__(self): return self.ants.__iter__()
    def __getitem__(self, *args): return self.ants.__getitem__(*args)
    def __setitem__(self, *args): return self.ants.__setitem__(*args)
    def __len__(self): return self.ants.__len__()
    def update(self):
        ArrayLocation.update(self)
        for a in self: a.update()
    def select_chans(self, active_chans=None):
        """Select which channels are used in computations.  Default is all."""
        for a in self: a.select_chans(active_chans)
        self.update()
    def ij2bl(self, i, j):
        """Convert baseline i,j (0 indexed) to Miriad's (i+1) << 8 | (j+1) 
        indexing scheme."""
        return ij2bl(i,j)
    
    def bl2ij(self, bl):
        """Convert Miriad's (i+1) << 8 | (j+1) baseline indexing scheme to 
        i,j (0 indexed)"""
        return bl2ij(bl)

    def bl_indices(self, auto=True, cross=True):
        """Return bl indices for baselines in the array."""
        if auto:
            if cross: return [self.ij2bl(i,j) 
                for i in range(len(self)) 
                for j in range(i,len(self))]
            else: return [self.ij2bl(i,i) for i in range(len(self))]
        else:
            if cross: return [self.ij2bl(i,j) 
                for i in range(len(self)) 
                for j in range(i+1,len(self))]
            else: return []
    def get_afreqs(self):
        """Return array of frequencies that are active for simulation."""
        return self[0].beam.afreqs
    def get_freqs(self):
        """Return array of (all) frequencies."""
        return self[0].beam.freqs
    def get_baseline(self, i, j, src='z'):
        """Return the baseline corresponding to i,j in various coordinate 
        projections: src='e' for current equatorial, 'z' for zenith 
        topocentric, 'r' for unrotated equatorial, or a RadioBody for
        projection toward that source."""
        bl = self[j] - self[i]
        if type(src) == str:
            if src == 'e': return np.dot(self._eq2now, bl)
            elif src == 'z': return np.dot(self._eq2zen, bl)
            elif src == 'r': return bl
            else: raise ValueError('Unrecognized source:' + src)
        try:
            if src.alt < 0:
                raise PointingError('%s below horizon' % src.src_name)
            m = src.map
        except(AttributeError):
            ra,dec = coord.eq2radec(src)
            m = coord.eq2top_m(self.sidereal_time() - ra, dec)
        return np.dot(m, bl).transpose()
    def get_phs_offset(self, i, j):
        """Return the frequency-dependent phase offset of baseline i,j."""
        return self[j].phsoff - self[i].phsoff
    def gen_uvw(self, i, j, src='z', w_only=False):
        """Compute uvw coordinates of baseline relative to provided RadioBody, 
        or 'z' for zenith uvw coordinates.  If w_only is True, only w (instead
        of (u,v,w) will be returned)."""
        x,y,z = self.get_baseline(i,j, src=src)
        afreqs = self.get_afreqs()
        afreqs = np.reshape(afreqs, (1,afreqs.size))
        if len(x.shape) == 0:
            if w_only: return z*afreqs
            else: return np.array([x*afreqs, y*afreqs, z*afreqs])
        #afreqs = np.reshape(afreqs, (1,afreqs.size))
        x.shape += (1,); y.shape += (1,); z.shape += (1,)
        if w_only: return np.dot(z,afreqs)
        else: return np.array([np.dot(x,afreqs), np.dot(y,afreqs), np.dot(z,afreqs)])
    def gen_phs(self, src, i, j, mfreq=.150, ionref=None, srcshape=None, 
            resolve_src=False):
        """Return phasing that is multiplied to data to point to src."""
        if ionref is None:
            try: ionref = src.ionref
            except(AttributeError): pass
        if not ionref is None or resolve_src: u,v,w = self.gen_uvw(i,j,src=src)
        else: w = self.gen_uvw(i,j,src=src, w_only=True)
        if not ionref is None: w += self.refract(u, v, mfreq=mfreq, ionref=ionref)
        o = self.get_phs_offset(i,j)
        phs = np.exp(-1j*2*np.pi*(w + o))
        if resolve_src:
            if srcshape is None:
                try: srcshape = src.srcshape
                except(AttributeError): pass
            if not srcshape is None: phs *= self.resolve_src(u, v, srcshape=srcshape)
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
        ru = a1 * (u*np.cos(th) - v*np.sin(th))
        rv = a2 * (u*np.sin(th) + v*np.cos(th))
        x = 2 * np.pi * np.sqrt(ru**2 + rv**2)
        # Use first Bessel function of the first kind (J_1)
        return np.where(x == 0, 1, 2 * _cephes.j1(x)/x).squeeze()
    def refract(self, u_sf, v_sf, mfreq=.150, ionref=(0.,0.)):
        """Calibrate a frequency-dependent source offset by scaling measured
        offsets at a given frequency.  Generates dw, a change in the
        projection of a baseline towards that source, which can be used to
        fix the computed phase of that source.
        ionref = (dra, ddec) where dra, ddec are angle offsets (in radians)
            of sources along ra/dec axes at the specified mfreq.
        u_sf,v_sf = u,v components of baseline, used to compute the
            change in w given angle offsets and the small angle approx.  Should
            be numpy arrays with sources (s) along the 1st axis and
            freqs (f) along the 2nd."""
        dra,ddec = ionref
        s,f = u_sf.shape
        try: dra.shape = (s,1)
        except(AttributeError): pass
        try: ddec.shape = (s,1)
        except(AttributeError): pass
        try: mfreq.shape = (s,1)
        except(AttributeError): pass
        f2 = self.get_afreqs()**2 ; f2.shape = (1, f)
        return (dra*u_sf + ddec*v_sf) * mfreq**2 / f2
    def phs2src(self, data, src, i, j, mfreq=.150, ionref=None, srcshape=None):
        """Apply phasing to zenith-phased data to point to src."""
        return data * self.gen_phs(src, i, j, 
            mfreq=mfreq, ionref=ionref, srcshape=srcshape, resolve_src=False)
    def unphs2src(self,data,src, i, j, mfreq=.150, ionref=None, srcshape=None):
        """Remove phasing from src-phased data to point to zenith."""
        return data / self.gen_phs(src, i, j,
            mfreq=mfreq, ionref=ionref, srcshape=srcshape, resolve_src=False)
