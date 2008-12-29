"""
Module for representing the geometry of an array of antennas and generating
phasing information.

Author: Aaron Parsons
Date: 11/07/2006
Revisions: 
    12/05/2006  arp Redid subclassing of ephem FixedBody, and added RadioSun.
    01/12/2007  arp Fixed H (hour angle) to be computed in sidereal time, not
                    solar.
    01/28/2007  arp Sync'd to latest miriad.py (0.0.4) version.  Hefty (4-5x)
                    speed up!  Also converted to use n-1.0.1
    05/14/2007  arp Fixed bug where coordinates were precessed to 2000, not
                    the current date.  
    02/06/2008  arp Set ephem pressure to 0 to remove optical distortion
                    correction.  Moved get_uvw_map to coord.eq2top_m.
"""
import ephem, math, numpy as n, coord, const

class PointingError(Exception):
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
    """Class representing coordinates of a celestial source."""
    def __init__(self, name='', **kwargs):
        self.src_name = name
    def compute(self, observer):
        """Update coordinates relative to the provided observer.  Must be
        called at each time step before accessing information."""
        # Generate a map for projecting baselines to uvw coordinates
        self.map = coord.eq2top_m(observer.sidereal_time()-self.ra, self.dec)
        # Generate vector pointing at src in various coordinates
        self.e_vec = coord.radec2eq((self.ra, self.dec))
        self.t_vec = coord.azalt2top((self.az, self.alt))
    def get_crd(self, crdsys, ncrd=3):
        """Return the coordinates of this location in the desired coordinate
        system ('eq','top') in the current epoch.  If ncrd=2, angular
        coordinates (ra/dec,az/alt) are returned, and if ncrd=3,
        xyz coordinates are returned."""
        assert(crdsys in ('eq','top'))
        assert(ncrd in (2,3))
        if crdsys == 'eq':
            if ncrd == 2: return (self.ra, self.dec)
            return self.e_vec
        else:
            if ncrd == 2: return (self.az, self.alt)
            return self.t_vec

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(ephem.FixedBody, RadioBody):
    """Class representing a source at fixed RA,DEC.  Combines ephem.FixedBody 
    with RadioBody."""
    def __init__(self, ra, dec, name='', **kwargs):
        """ra = source's right ascension (epoch=J2000)
        dec = source's declination (epoch=J2000)"""
        RadioBody.__init__(self, name=name)
        ephem.FixedBody.__init__(self)
        self._ra, self._dec = ra, dec
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
    """Class representing moving sources (Sun,Moon,planets).  Combines ephem
    versions of these objects with RadioBody."""
    def __init__(self, name, **kwargs):
        """name is used to lookup appropriate ephem celestial object."""
        RadioBody.__init__(self, name=name)
        self.Body = eval('ephem.%s()' % name)
    def __getattr__(self, nm):
        """When getting an attribute, first try to get it from this class,
        and if that fails, try to get it from the underlying ephem object."""
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Body.__getattribute__(nm)
    def __setattr__(self, nm, val):
        """When setting an attribute, first try to set it for this class,
        and if that fails, try to set it for the underlying ephem object."""
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
    """Class for holding a catalog of celestial sources."""
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
    """Representation of physical attributes of individual antenna in array."""
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
        """Select only enumerated channels to use for future calculations."""
        self.beam.select_chans(active_chans)
        self.offset = n.polyval(self._offset, self.beam.afreqs)
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
    """Representation of location and time of an observation."""
    def __init__(self, location):
        """location = (lat, long, [elev]) of array"""
        ephem.Observer.__init__(self)
        self.pressure = 0
        self.update_location(location)
    def update_location(self, location):
        """Initialize antenna array wth provided location.  May be (lat, long) 
        or (lat, long, elev)."""
        if len(location) == 2: self.lat, self.long = location
        else: self.lat, self.long, self.elev = location
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
    """Representation of collection of antennas, their spacings, and
       location/time of observations."""
    def __init__(self, location, ants, **kwargs):
        """ location = (lat, long, [elev]) of array
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
        bls,dlys,offs = [],[],[]
        self.bl_order = {}
        # Calculate deltas between antenna pairs
        for i, ai in enumerate(self.ants):
            for j, aj in enumerate(self.ants[i:]):
                bl = self.ij2bl(i, j+i)
                bls.append(aj - ai)
                dlys.append(aj.delay - ai.delay)
                offs.append(aj.offset - ai.offset)
                self.bl_order[bl] = len(bls) - 1
        self.bls,self.dlys,self.offs = n.array(bls),n.array(dlys),n.array(offs)
        # Compute (static) zenith baselines
        m = coord.eq2top_m(0., self.lat)
        self.zbls = n.dot(m, self.bls.transpose()).transpose()
    def set_ephemtime(self, t=None):
        """Set current time as derived from the ephem package.  Recalculates
        matrix for projecting baselines into current positions."""
        ArrayLocation.set_ephemtime(self, t=t)
        # Rotate baselines in eq coords to current sidereal time
        m = coord.rot_m(-self.sidereal_time(), n.array([0.,0.,1.]))
        self.ebls = n.dot(m, self.bls.transpose()).transpose()
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
        b = self.bl_order[self.ij2bl(i,j)]
        if type(src) == str:
            if src == 'e': return self.ebls[b]
            elif src == 'z': return self.zbls[b]
            elif src == 'r': return self.bls[b]
            else: raise ValueError('Unrecognized source:' + src)
        if src.alt < 0: raise PointingError('%s below horizon' % src.src_name)
        return n.dot(src.map, self.bls[b])
    def get_delay(self, i, j):
        """Return the delay corresponding to baseline i,j."""
        return self.dlys[self.bl_order[self.ij2bl(i,j)]]
    def get_offset(self, i, j):
        """Return the delay corresponding to baseline i,j."""
        return self.offs[self.bl_order[self.ij2bl(i,j)]]
    def gen_uvw(self, i, j, src='z'):
        """Compute uvw coordinates of baseline relative to provided RadioBody, 
        or 'z' for zenith uvw coordinates."""
        xyz = self.get_baseline(i,j, src=src)
        afreqs = self.ants[0].beam.afreqs
        if len(xyz.shape) == 1: xyz = n.resize(xyz, (afreqs.size, xyz.size))
        return xyz * n.reshape(afreqs, (afreqs.size, 1))
    def gen_phs(self, src, i, j, angsize=None):
        """Return phasing that is multiplied to data to point to src.  If
        angsize is provided, amplitudes will be adjusted to reflect resolution
        effects for a uniform disc of the provided angular radius (radians)."""
        try: src = src.e_vec
        except(AttributeError): pass
        bl = self.get_baseline(i,j,src='e')
        bl_dot_s = n.dot(src, bl)
        if len(bl_dot_s.shape) >= 1: bl_dot_s.shape = (bl_dot_s.size, 1)
        t = self.get_delay(i,j)
        o = self.get_offset(i,j)
        afreqs = self.ants[0].beam.afreqs
        afreqs = n.reshape(afreqs, (1,afreqs.size))
        # Rudimentarily account for resolution effects
        if not angsize is None:
            angsize.shape = bl_dot_s.shape
            amp = angsize * n.sqrt(n.dot(bl,bl) - bl_dot_s**2)
            amp = n.sinc(n.dot(amp, afreqs))
        else: amp = 1.
        phs = amp * n.exp(-1j*2*n.pi*(n.dot(bl_dot_s + t, afreqs) + o))
        return phs.squeeze()
    def phs2src(self, data, src, i, j):
        """Apply phasing to zenith-phased data to point to src."""
        return data * self.gen_phs(src, i, j)
    def unphs2src(self, data, src, i, j):
        """Remove phasing from src-phased data to point to zenith."""
        return data * n.conjugate(self.gen_phs(src, i, j))
    def rmsrc(self, data, srcs, i, j, swath=0, norm_bandpass=None):
        """Remove src flux from data by phasing to each src and removing signal
        within swath bins of 0 in lag domain (Fourier transform of the 
        frequency axis)."""
        if type(srcs) != list: srcs = [srcs]
        for src in srcs:
            try: phs = self.gen_phs(src, i, j)
            except(PointingError): continue
            data *= phs
            if swath == 0:
                if norm_bandpass is None: data -= n.ma.average(data)
                else: data -= n.ma.sum(data) * norm_bandpass
            else:
                img = n.fft.ifft(data.filled(0))
                img[swath:-swath] = 0
                data -= n.fft.fft(img)
            data /= phs
        return data
