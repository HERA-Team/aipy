"""
Module for representing an array of antennas.  Two different coordinate
systems are employed in this module:
Topocentric (measured in cm):
    x axis: +N (-S)
    y axis: +E (-W)
    z axis: elevation
Equatorial (measured in nanoseconds):
    x axis: radial (plane of equator)
    y axis: +E (-W)
    z axis: +N celestial pole (-S)

Author: Aaron Parsons
Date: 11/07/2006
Revisions: 
    12/05/2006  arp Redid subclassing of ephem FixedBody, and added RadioSun.
    01/12/2007  arp Fixed H (hour angle) to be computed in sidereal time, not
                    solar.
    01/28/2007  arp Sync'd to latest miriad.py (0.0.4) version.  Hefty (4-5x)
                    speed up!  Also converted to use numpy-1.0.1
    05/14/2007  arp Fixed bug where coordinates were precessed to 2000, not
                    the current date.  
"""
import pyephem as ephem, math, numpy
import const

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
    """Convert Julian date to an ephem date, which is measured from noon,
    Jan. 1, 1900."""
    return ephem.date(num - 2415020)

def ephem2juldate(num):
    """Convert ephem date (measured from noon, Jan. 1, 1900) to a Julian 
    date."""
    return float(num + 2415020)

def gen_uvw_map(H, d):
    """Return a mapping to project uvw to a given hour angle (H) and dec (d).
    Traded speed of math.sin for numpy.sin to handle freq-dependent H, d (due
    to ionosphere refraction)."""
    sin_H, cos_H = numpy.sin(H), numpy.cos(H)
    sin_d, cos_d = numpy.sin(d), numpy.cos(d)
    zero = numpy.zeros_like(H)
    map =  numpy.array([[    sin_H    ,       cos_H  ,       zero  ],
                        [ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
                        [  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])
    if len(map.shape) == 3: map = map.transpose([2, 0, 1])
    return map

def azalt2top(az, alt):
    """Convert azimuth, altitude into topocentric coordinates."""
    sin_alt = math.sin(alt); cos_alt = math.cos(alt)
    sin_az = math.sin(az); cos_az = math.cos(az)
    return numpy.array([cos_az*cos_alt, sin_az*cos_alt, sin_alt])

#  ____           _ _       ____            _       
# |  _ \ __ _  __| (_) ___ | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \|  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|____/ \___/ \__,_|\__, |
#                                             |___/ 

class RadioBody:
    """A class redefining ephem's sense of brightness for radio astronomy."""
    def __init__(self, name='', **kwargs):
        self.src_name = name
        self._obs_cache = None
    def in_cache(self, observer, refract):
        try: ocache = observer.cache
        except(AttributeError): return False
        if refract: ocache = -ocache
        if self._obs_cache == ocache: return True
    def cache(self, observer, refract):
        try: ocache = observer.cache
        except(AttributeError): return
        if refract: ocache = -ocache
        self._obs_cache = ocache
    def gen_uvw_map(self, observer, refract=False):
        """Generate a uvw map useful for projecting baselines.  Can compensate
        for ionospheric refraction for a given plasma frequency f_c."""
        H = float(observer.sidereal_time() - self.ra)
        d = self.dec
        if refract:
            afreqs = observer.ants[0].beam.afreqs
            H, d = self.add_refraction(H, d, afreqs, observer.lat, self.f_c)
        map = gen_uvw_map(H, d)
        return map
    def add_refraction(self, H, d, f, lat,
            f_c=.012, d_ion=200e5, R_e=6000e5):
        """Modify H (hour angle) and d (declination) for ionospheric refraction,
        given an observing frequency (f), observer latitude (lat), plasma
        frequency (f_c in GHz), ionospheric height (d_ion), and earth radius
        (R_e).  Complicated math to change a simple Z angle refraction into
        H and d perturbations.  For now, only compensates for spherically
        symmetric ionosphere.  See 1982PASAu...4..386S."""
        Z = math.pi/2 - self.alt    # Change alt into zenith angle
        A = self.az                 # Azimuth
        sin_Z, cos_Z, tan_Z = math.sin(Z), math.cos(Z), math.tan(Z)
        sin_lat, cos_lat = math.sin(lat), math.cos(lat)
        sin_A, cos_A = math.sin(A), math.cos(A)
        sin_d, cos_d = math.sin(d), math.cos(d)
        sin_H, cos_H = math.sin(H), math.cos(H)
        # Change in zenith angle due to spherically symmetrical ionosphere
        d_Zs = (f_c/f)**2 * (d_ion/R_e) * tan_Z / cos_Z**2
        # Change in zenith angle due to large scale horizontal density gradients
        d_Zw = 0 # ignoring wedge component of zen_ang for refraction for now
        # Total zenith angle displacement is sum of 2 components
        d_Z = d_Zs + d_Zw
        # Change in az angle due to large scale horizontal density gradients
        d_Aw = 0 # ignoring wedge component of azimuth for refraction for now
        # Total azimuth angle displacement is only from wedge component
        d_A = d_Aw
        # Convert zenith angle displacement into declination component
        d_d = (d_Z*(-sin_Z*sin_lat + cos_Z*cos_lat*cos_A) \
                - d_A*sin_Z*cos_lat*sin_A) / cos_d
        # Convert zenith angle displacement into hour-angle component
        d_H = (d_d*sin_H*sin_d - d_A*cos_A*sin_Z - d_Z*sin_A*cos_Z) \
                / (cos_H*cos_d)
        return H + d_H, d + d_d
        

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(ephem.FixedBody, RadioBody):
    """A class combining ephem's FixedBody with a RadioBody."""
    def __init__(self, ra, dec, f_c=.012, name='', **kwargs):
        """ra:           source's right ascension (epoch=2000)
        dec:          source's declination (epoch=2000)
        f_c:          the plasma frequency (a measure of the electron density)
                      in the direction of the source, which can cause the 
                      source to appear displaced from its quiescent position 
                      as a function of frequency."""
        RadioBody.__init__(self, name=name)
        ephem.FixedBody.__init__(self)
        self._ra, self._dec, self.f_c = ra, dec, f_c
    def compute(self, observer, refract=False):
        if self.in_cache(observer, refract): return
        ephem.FixedBody.compute(self, observer)
        self.map = self.gen_uvw_map(observer, refract=refract)
        self.cache(observer, refract)

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(RadioBody, object):
    """A class combining ephem's Sun, Moon, planets, etc. with a RadioBody."""
    def __init__(self, name, f_c=.012, **kwargs):
        """f_c:          the plasma freq (a measure of the electron density)
                      in the direction of the source, which can cause the 
                      source to appear displaced from its quiescent position 
                      as a function of frequency."""
        RadioBody.__init__(self, name=name)
        self.Body = eval('ephem.%s()' % name)
        self.f_c = f_c
    def __getattr__(self, nm):
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Body.__getattribute__(nm)
    def __setattr__(self, nm, val):
        try: object.__setattr__(self, nm, val)
        except(AttributeError): return setattr(self.Body, nm, val)
    def compute(self, observer, refract=False):
        if self.in_cache(observer, refract): return
        self.Body.compute(observer)
        self.map = self.gen_uvw_map(observer, refract=refract)
        self.cache(observer, refract)

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/ 

class SrcCatalog(dict):
    """A class for holding celestial sources."""
    def __init__(self, srcs, **kwargs):
        dict.__init__(self)
        for s in srcs: self.add_src(s)
        self._obs_cache = None
    def add_src(self, src):
        self[src.src_name] = src
    def in_cache(self, observer, refract):
        try: ocache = observer.cache
        except(AttributeError): return False
        if refract: ocache = -ocache
        if self._obs_cache == ocache: return True
    def cache(self, observer, refract):
        try: ocache = observer.cache
        except(AttributeError): return
        if refract: ocache = -ocache
        self._obs_cache = ocache
    def get_srcs(self, *args):
        return [self[s] for s in args]
    def compute(self, observer, refract=False):
        if self.in_cache(observer, refract): return
        for s in self: self[s].compute(observer, refract=refract)
        self.cache(observer, refract)

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam:
    """Template for an antenna's beam pattern on the sky.  Beams are also
    responsible for holding information about which frequencies are active
    (i.e. Antenna and AntennaArray access frequencies through the Beam)."""
    def __init__(self, freqs, active_chans=None, **kwargs):
        """freqs:        frequencies (in GHz) at bin centers across spectrum
        active_chans: channels to be selected for future freq calculations"""
        self.freqs = freqs
        self.select_chans(active_chans)
    def select_chans(self, active_chans=None):
        """Choose only 'active_chans' for future freq calculations."""
        if active_chans is None: active_chans = numpy.arange(self.freqs.size)
        self.chans = active_chans
        self.afreqs = self.freqs.take(active_chans)

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna:
    """A representation of the physical location an individual antenna in 
    an array, and possibly a systematic delay associated with it."""
    def __init__(self, x, y, z, beam, delay=0., offset=0., **kwargs):
        """x, y, z:    Antenna coordinates in equatorial (ns) coordinates"""
        self.pos = numpy.array((x,y,z), numpy.float64) # must be float64 for mir
        self.beam = beam
        self.delay = delay
        self.offset = offset
    def select_chans(self, active_chans=None):
        self.beam.select_chans(active_chans)
    def __tuple__(self):
        return (self.pos[0], self.pos[1], self.pos[2])
    def __list__(self):
        return [self.pos[0], self.pos[1], self.pos[2]]
    def __add__(self, a):
        return self.pos + a.pos
    __radd__ = __add__
    def __neg__(self):
        return -self.pos
    def __sub__(self, a):
        return self.pos - a.pos
    def __rsub__(self, a):
        return a.pos - self.pos

#     _                         _                    _   _             
#    / \   _ __ _ __ __ _ _   _| |    ___   ___ __ _| |_(_) ___  _ __  
#   / _ \ | '__| '__/ _` | | | | |   / _ \ / __/ _` | __| |/ _ \| '_ \ 
#  / ___ \| |  | | | (_| | |_| | |__| (_) | (_| (_| | |_| | (_) | | | |
# /_/   \_\_|  |_|  \__,_|\__, |_____\___/ \___\__,_|\__|_|\___/|_| |_|
#                         |___/                                        

class ArrayLocation(ephem.Observer):
    """Collected information about where and when an array is."""
    def __init__(self, location=None, uv=None):
        """location:   location of the array in (lat, long, [elev])
        uv:         Miriad UV file"""
        ephem.Observer.__init__(self)
        self.cache = 0
        if not uv is None: self.from_uv(uv)
        else:
            if location is None:
                raise ValueError('Must provide either uv or location.')
            self.update_location(location)
    def update_location(self, location):
        """Initialize the antenna array for the provided location.  Locations
        may be (lat, long) or (lat, long, elev)."""
        if len(location) == 2: self.lat, self.long = location
        else: self.lat, self.long, self.elev = location
        # Establish conversion between equatorial and topocentric
        self._top2eq_map = numpy.array(
            [[-math.sin(self.lat), 0, math.cos(self.lat)],
             [          0        , 1,          0        ],
             [ math.cos(self.lat), 0, math.sin(self.lat)]])
        self.cache += 1
    def set_jultime(self, t=None):
        """Set the current time to a Julian date."""
        if t is None: t = ephem.julian_date()
        self.set_ephemtime(juldate2ephem(t))
    def set_ephemtime(self, t=None):
        """Set the current time to a time derived from the ephem package."""
        if t is None: t = ephem.now()
        self.date = t
        self.epoch = t
        self.cache += 1
    def src_eq_vec(self, body):
        """Return an equatorial vector pointing at body."""
        return self.top2eq(azalt2top(body.az, body.alt), no_norm=True)
    def from_uv(self, uv):
        """Update location from 'latitud' and 'longitu' in Miriad UV file."""
        location = (uv['latitud'], uv['longitu'])
        self.update_location(location)
    def top2eq(self, *args, **kwargs):
        """Convert topocentric antenna coordinates to equatorial coordinates,
        given the current latitude."""
        if len(args) == 3: top_vec = numpy.array(args)
        elif len(args) == 1: top_vec = args[0]
        else: raise ValueError('Wrong number of arguments.')
        rv = numpy.dot(self._top2eq_map, top_vec)
        if kwargs.has_key('no_norm'): return rv
        else: return rv / const.len_ns
    def eq2top(self, *args, **kwargs):
        """Convert equatorial antenna coordinates to topocentric,
        given the current latitude."""
        if len(args) == 3: eq_vec = numpy.array(args)
        elif len(args) == 1: eq_vec = args[0]
        else: raise ValueError('Wrong number of arguments.')
        rv = numpy.dot(numpy.transpose(self._top2eq_map), eq_vec)
        if kwargs.has_key('no_norm'): return rv
        else: return rv * const.len_ns

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(ArrayLocation):
    """A representation of a collection of antennas, their spacings, and
       information about when and where the array is."""
    def __init__(self, ants=None, location=None, uv=None, **kwargs):
        """ants:       a list of Antenna instances
        location:   location of the array in (lat, long, [elev])
        uv:         Miriad UV file"""
        ArrayLocation.__init__(self, location=location, uv=uv)
        if not uv is None: self.from_uv(uv)
        else:
            if ants is None:
                raise ValueError('Must provide either ants or uv.')
            self.update_antennas(ants)
            self.select_chans()
    def from_uv(self, uv):
        """Update antenna positions, array location, and frequencies from 
        Miriad UV file."""
        ArrayLocation.from_uv(self, uv)
        # Generate frequency information
        sfreq = uv['sfreq']
        sdf = uv['sdf']
        freqs = numpy.arange(uv['nchan'], dtype=numpy.float) * sdf + sfreq
        beam = Beam(freqs)
        # Generate antenna positions
        ants = uv['antpos']
        ants.shape = (3, uv['nants'])
        ants = ants.transpose()
        # Should get delay information..., and what about offsets?
        ants = [Antenna(x,y,z, beam=beam, delay=0.) for x,y,z in ants]
        self.update_antennas(ants)
    def update_antennas(self, ants):
        """Update the antenna array to use the provided list of antennas."""
        self.ants = ants
        self.update()
    def update(self):
        """If antenna parameters or active channels have been changed, this
        function updates variables derived from them.  Increments the cache
        value to force recalculation of cached values in source projections."""
        bls = []    # Baseline vectors (in ns equatorial coordinates)
        dlys = []
        offs = []
        self.baseline_order = {}
        for i, ai in enumerate(self.ants):
            for j, aj in enumerate(self.ants[i:]):
                j += i
                bls.append(aj - ai)
                dlys.append(aj.delay - ai.delay)
                offs.append(aj.offset - ai.offset)
                bl = self.ij2bl(i, j)
                self.baseline_order[bl] = len(bls) - 1
        self.baselines = numpy.array(bls)
        self.delays = numpy.array(dlys)
        self.offsets = numpy.array(offs)
        self.cache += 1
    def select_chans(self, active_chans=None):
        """Select which channels are used in computations.  Default is all."""
        for a in self.ants: a.select_chans(active_chans)
        self.update()
    def ij2bl(self, i, j=None):
        """Convert from i,j (counting from 0) baseline notation to
        (i+1) << 8 | (j+1) baseline notation.  If j is not provided, assume
        i is in the desired baseline notation."""
        if j is None: return int(i)
        else: return (int(i)+1) << 8 | (int(j)+1)
    def bl2ij(self, bl):
        """Convert from (i+1) << 8 | (j+1) baseline notation to
        i,j (counting from 0) baseline notation."""
        bl = int(bl)
        return ((bl >> 8) & 255) - 1, (bl & 255) - 1
    def get_baseline(self, i, j):
        """Return the baseline corresponding to i,j."""
        return self.baselines[self.baseline_order[self.ij2bl(i,j)]]
    def get_delay(self, i, j):
        """Return the delay corresponding to i,j."""
        return self.delays[self.baseline_order[self.ij2bl(i,j)]]
    def get_offset(self, i, j):
        """Return the delay corresponding to i,j."""
        return self.offsets[self.baseline_order[self.ij2bl(i,j)]]
    def get_projected_baseline(self, i, j, src=None):
        """Project equatorial baselines toward a source position."""
        # If no body is provided to point to, return zenith baselines
        if src is None: uvw_map = gen_uvw_map(0., self.lat)
        else:
            if src.alt < 0:
                raise PointingError('%s is below horizon' % src.src_name)
            uvw_map = src.map
        # Counting on freqs to not change often (otherwise body's cached map
        # will raise an error in the following dot product).
        proj_bl = numpy.dot(uvw_map, self.get_baseline(i, j))
        return proj_bl
    def gen_phs(self, src, i, j, with_coord=False):
        """Return the phasing to be multiplied to data to point to src."""
        xyz = self.get_projected_baseline(i, j, src=src)
        z = xyz[...,2]
        t = self.get_delay(i, j)
        o = self.get_offset(i, j)
        afreqs = self.ants[0].beam.afreqs
        phs = numpy.exp(-1j*(2*numpy.pi * (z+t) * afreqs + o))
        if with_coord:
            if len(xyz.shape) == 1:
                xyz = numpy.resize(xyz, (len(afreqs), xyz.size))
            return phs, xyz * numpy.reshape(afreqs, (afreqs.size, 1))
        else: return phs
    def phs2src(self, data, src, i, j, with_coord=False):
        """Apply phasing to zenith data to point to src."""
        rv = self.gen_phs(src, i, j, with_coord=with_coord)
        if with_coord: return data * rv[0], rv[1]
        else: return data * rv
    def unphs2src(self, data, src, i, j):
        """Remove phasing from src data to point to zenith."""
        return data * numpy.conjugate(self.gen_phs(src, i, j))
    def rmsrc(self, data, srcs, i, j, swath=0, norm_bandpass=None):
        """Remove src flux from data.  Can aggressively remove srcs (allowing
        for inexact calibration) by removing 'swath' adjacent bins in 
        lag space."""
        if type(srcs) != list: srcs = [srcs]
        for src in srcs:
            try: phs = self.gen_phs(src, i, j)
            except(PointingError): continue
            data *= phs
            if swath == 0:
                if norm_bandpass is None: data -= numpy.ma.average(data)
                else: data -= numpy.ma.sum(data) * norm_bandpass
            else:
                img = numpy.fft.ifft(data.filled(0))
                img[swath:-swath] = 0
                data -= numpy.fft.fft(img)
            data /= phs
        return data
