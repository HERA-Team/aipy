"""
Module for representing an array of antennas.  Two different coordinate
systems hare employed in this module:
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
import ephem, math, numpy
import constants

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
    """Return a mapping to project uvw to a given hour angle (H) and dec (d)."""
    sin_H, cos_H = math.sin(H), math.cos(H)
    sin_d, cos_d = math.sin(d), math.cos(d)
    return numpy.array([[    sin_H    ,       cos_H  ,         0   ],
                        [ -sin_d*cos_H,   sin_d*sin_H,      cos_d  ],
                        [  cos_d*cos_H,  -cos_d*sin_H,      sin_d  ]])

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
    def __init__(self):
        self.prev_sidereal_time = 0     # Used to avoid redundant map calc.
    def gen_uvw_map(self, observer):
        """Generate a uvw map useful for projecting baselines."""
        self.compute(observer)
        if self.alt < 0: raise PointingError('%s is below horizon' % self.name)
        t = observer.sidereal_time()
        if t != self.prev_sidereal_time:
            self.prev_sidereal_time = t
            H = float(t - self.ra)
            d = self.dec
            self.map = gen_uvw_map(H, d)
        return self.map

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(ephem.FixedBody, RadioBody):
    """A class combining ephem's FixedBody with a RadioBody."""
    def __init__(self, ra, dec, name=''):
        """ra:           source's right ascension
        dec:          source's declination"""
        RadioBody.__init__(self)
        ephem.FixedBody.__init__(self)
        self._ra = ra
        self._dec = dec
        self.name = name
    def __repr__(self):
        """Return a string which can be used to hash src."""
        return 'RadioFixedBody(%s,%s,name=%s)' % (self._ra,self._dec,self.name)

#  ____           _ _      ____              
# |  _ \ __ _  __| (_) ___/ ___| _   _ _ __  
# | |_) / _` |/ _` | |/ _ \___ \| | | | '_ \ 
# |  _ < (_| | (_| | | (_) |__) | |_| | | | |
# |_| \_\__,_|\__,_|_|\___/____/ \__,_|_| |_|

class RadioSun(RadioBody, object):
    """A class combining ephem's Sun with a RadioBody."""
    def __init__(self):
        RadioBody.__init__(self)
        self.Sun = ephem.Sun()
    def __getattr__(self, nm):
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Sun.__getattribute__(nm)
    def __setattr__(self, nm, val):
        try: object.__setattr__(self, nm, val)
        except(AttributeError): return setattr(self.Sun, nm, val)
    def __repr__(self):
        return 'RadioSun()'

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna:
    """A representation of the physical location an individual antenna in 
    an array, and possibly a systematic delay associated with it."""
    def __init__(self, x, y, z, delay=0.):
        """x, y, z:    Antenna coordinates in equatorial (ns) coordinates"""
        self.pos = numpy.array((x,y,z), numpy.float64) # must be float64 for mir
        self.delay = delay
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
    def set_jultime(self, t=None):
        """Set the current time to a Julian date."""
        if t is None: t = ephem.julian_date()
        self.set_ephemtime(juldate2ephem(t))
    def set_ephemtime(self, t=None):
        """Set the current time to a time drived from the ephem package."""
        if t is None: t = ephem.now()
        self.date = t
        self.epoch = t
    def src_eq_vec(self, body):
        """Return an equatorial vector pointing at body."""
        body.compute(self)
        return self.top2eq(azalt2top(body.az, body.alt), no_norm=True)
    def from_uv(self, uv):
        """Update location from 'latitud' and 'longitu' in Miriad UV file."""
        location = (uv['latitud'], uv['longitu'])
        self.update_location(location)

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(ArrayLocation):
    """A representation of a collection of antennas, their spacings, and
       information about when and where the array is."""
    def __init__(self, antennas=None, location=None, freqs=None, uv=None):
        """antennas:   a list of Antenna instances
        location:   location of the array in (lat, long, [elev])
        uv:         Miriad UV file"""
        ArrayLocation.__init__(self, location=location, uv=uv)
        if not uv is None: self.from_uv(uv)
        else:
            if antennas is None:
                raise ValueError('Must provide either antennas or uv.')
            self.update_antennas(antennas)
            self.freqs = freqs
            self.set_active_chans()
    def from_uv(self, uv):
        """Update antenna positions, array location, and frequencies from 
        Miriad UV file."""
        ArrayLocation.from_uv(self, uv)
        antennas = uv['antpos']
        antennas.shape = (3, uv['nants'])
        antennas = antennas.transpose()
        self.update_antennas(antennas)
        sfreq = uv['sfreq']
        sdf = uv['sdf']
        self.freqs = numpy.arange(uv['nchan'], dtype=numpy.float) * sdf + sfreq
        self.set_active_chans()
    def set_active_chans(self, chans=None):
        """Select which channels are used in computations.  Default is all."""
        if chans is None: self.active_freqs = self.freqs
        else: self.active_freqs = numpy.take(self.freqs, chans)
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
    def update_antennas(self, antennas):
        """Initialize the antenna array using a list of antennas.  Generates
        zenith baselines."""
        self.antennas = antennas
        self.n_ants = len(antennas)
        bls = []
        dlys = []
        self.baseline_order = {}
        for i in range(self.n_ants):
            for j in range(i, self.n_ants):
                bls.append(antennas[j] - antennas[i])
                try: dlys.append(antennas[j].delay - antennas[i].delay)
                except(AttributeError): dlys.append(0)
                bl = self.ij2bl(i, j)
                self.baseline_order[bl] = len(bls) - 1
        self.baselines = numpy.array(bls)
        self.delays = numpy.array(dlys)
    def get_baseline(self, i, j=None):
        """Return the baseline corresponding to i,j (see ij2bl for details)."""
        bl = self.ij2bl(i, j)
        return self.baselines[self.baseline_order[bl]]
    def top2eq(self, *args, **kwargs):
        """Convert topocentric antenna coordinates to equatorial coordinates,
        given the current latitude."""
        if len(args) == 3: top_vec = numpy.array(args)
        elif len(args) == 1: top_vec = args[0]
        else: raise ValueError('Wrong number of arguments.')
        rv = numpy.dot(self._top2eq_map, top_vec)
        if kwargs.has_key('no_norm'): return rv
        else: return rv / constants.len_ns
    def eq2top(self, *args, **kwargs):
        """Convert equatorial antenna coordinates to topocentric,
        given the current latitude."""
        if len(args) == 3: eq_vec = numpy.array(args)
        elif len(args) == 1: eq_vec = args[0]
        else: raise ValueError('Wrong number of arguments.')
        rv = numpy.dot(numpy.transpose(self._top2eq_map), eq_vec)
        if kwargs.has_key('no_norm'): return rv
        else: return rv * constants.len_ns
    def get_projected_baseline(self, i, j=None, body=None):
        """Rotate zenith (the default pointing) to the actual
        pointing, return the projected baseline i,j."""
        # If no body is provided to point to, return zenith baselines
        if body is None: uvw_map = gen_uvw_map(0., self.lat)
        else: uvw_map = body.gen_uvw_map(self)
        proj_bl = numpy.dot(uvw_map, self.get_baseline(i, j))
        return proj_bl
    def get_delay(self, i, j=None):
        """Return the delay corresponding to i,j (see ij2bl for details)."""
        bl = self.ij2bl(i, j)
        return self.delays[self.baseline_order[bl]]
    def gen_phs(self, src, i, j=None):
        """Return the phasing which must be applied to data to point to src."""
        bl = self.ij2bl(i, j)
        z = self.get_projected_baseline(bl, body=src)[2]
        t = self.get_delay(bl)
        return numpy.exp(-2*numpy.pi*1j * (z+t) * self.active_freqs)
    def phs2src(self, data, src, i, j=None):
        """Apply phasing to zenith data to point to src."""
        return data * self.gen_phs(src, i, j)
    def unphs2src(self, data, src, i, j=None):
        """Remove phasing from src data to point to zenith."""
        return data * numpy.conjugate(self.gen_phs(src, i, j))
    def rmsrc(self, data, src, i, j=None, swath=0):
        """Remove src flux from data.  Can aggressively remove srcs (allowing
        for inexact calibration) by removing 'swath' adjacent bins in 
        lag space."""
        try: phs = self.gen_phs(src, i, j)
        except(PointingError): return data
        d = data * phs
        if swath == 0: d -= numpy.ma.average(d)
        else:
            img = numpy.fft.ifft(d.filled(0))
            img[swath:-swath] = 0
            d -= numpy.fft.fft(img)
        d /= phs
        return d
