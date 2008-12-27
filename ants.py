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
    def __init__(self, strength, freqs, meas_freq, spec_index, 
            active_chans=None, t0=0):
        """strength:     source flux measured at 'meas_freq'
        freqs:        frequencies (in GHz) at bin centers across spectrum
        meas_freq:    frequency (in GHz) where 'strength' was measured
        spec_index:   index of power-law spectral model of source emission
        active_chans: channels to be selected for future freq calculations"""
        self.freqs = freqs
        self._meas_freq = meas_freq
        self._strength = strength
        self._t0 = t0
        self._spec_index = spec_index
        self.select_chans(active_chans)
        self.prev_sidereal_time = 0     # Used to avoid redundant map calc.
    def select_chans(self, active_chans):
        """Choose only 'active_chans' for future freq calculations."""
        if active_chans is None: active_chans = numpy.arange(self.freqs.size)
        self.chans = numpy.array(active_chans)
        self.update(self._strength, self._spec_index)
    def update(self, strength, spec_index):
        """Update parameters for source strength and spectral index."""
        try: len(strength)
        except(TypeError): strength = [strength]
        self._strength = strength
        try: len(spec_index)
        except(TypeError): spec_index = [spec_index]
        self._spec_index = spec_index
        self._emission = (self.freqs.take(self.chans) / self._meas_freq)
    def emission(self, observer):
        # This may be a redundant compute
        # Could optimize to buffer result for same times
        # or if strength, spec_index polynomials are of order 1
        self.compute(observer)
        t = observer.date - self._t0
        cur_strength = numpy.polyval(self._strength, t)
        cur_spec_index = numpy.polyval(self._spec_index, t)
        return cur_strength * self._emission**cur_spec_index
    def gen_uvw_map(self, observer):
        """Generate a uvw map useful for projecting baselines."""
        self.compute(observer)
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
    def __init__(self, ra, dec, freqs, 
            strength=1., meas_freq=.150, spec_index=-1., active_chans=None):
        """ra:           source's right ascension
        dec:          source's declination
        freqs:        frequencies (in GHz) at bin centers across spectrum
        strength:     source flux measured at 'meas_freq'
        meas_freq:    frequency (in GHz) where 'strength' was measured
        spec_index:   index of power-law spectral model of source emission
        active_chans: channels to be selected for future freq calculations"""
        ephem.FixedBody.__init__(self)
        self._ra = ra
        self._dec = dec
        RadioBody.__init__(self, strength, freqs, meas_freq, spec_index, 
            active_chans=active_chans)
    def __repr__(self):
        """This string is used to hash this source to avoid redundant
        computations."""
        return 'RadioFixedBody(%f, %f, strength=%s, spec_index=%s)' % \
            (self._ra, self._dec, self._strength, self._spec_index)

#  ____           _ _      ____              
# |  _ \ __ _  __| (_) ___/ ___| _   _ _ __  
# | |_) / _` |/ _` | |/ _ \___ \| | | | '_ \ 
# |  _ < (_| | (_| | | (_) |__) | |_| | | | |
# |_| \_\__,_|\__,_|_|\___/____/ \__,_|_| |_|

class RadioSun(RadioBody, object):
    """A representation of the Sun in radio.  Uses blackbelt kung-fu to 
    subclass ephem.Sun()."""
    def __init__(self, freqs, strength=1., meas_freq=.150, spec_index=-1.,
            active_chans=None):
        """freqs:        frequencies (in GHz) at bin centers across spectrum
        strength:     source flux measured at 'meas_freq'
        meas_freq:    frequency (in GHz) where 'strength' was measured
        spec_index:   index of power-law spectral model of source emission
        active_chans: channels to be selected for future freq calculations"""
        RadioBody.__init__(self, strength, freqs, meas_freq, spec_index, 
            active_chans=active_chans)
        self.Sun = ephem.Sun()
    def __getattr__(self, nm):
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Sun.__getattribute__(nm)
    def __setattr__(self, nm, val):
        try: object.__setattr__(self, nm, val)
        except(AttributeError): return setattr(self.Sun, nm, val)
    def __repr__(self):
        return 'RadioSun(%f, %f, strength=%s, spec_index=%s)' % \
            (self._ra, self._dec, self._strength, self._spec_index)

#  ____                           _     _     _   
# / ___|  ___  _   _ _ __ ___ ___| |   (_)___| |_ 
# \___ \ / _ \| | | | '__/ __/ _ \ |   | / __| __|
#  ___) | (_) | |_| | | | (_|  __/ |___| \__ \ |_ 
# |____/ \___/ \__,_|_|  \___\___|_____|_|___/\__|

class SourceList:
    """A class for holding celestial sources."""
    def __init__(self, src_dict, active_chans=None):
        self.names = src_dict.keys()
        self.names.sort()
        self.sources = [src_dict[s] for s in self.names]
        self.select_chans(active_chans)
    def select_chans(self, active_chans):
        """Choose only 'active_chans' for future freq calculations."""
        for s in self.sources: s.select_chans(active_chans)

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

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class ArrayLocation(ephem.Observer):
    """Collected information about where and when an array is."""
    def __init__(self, location):
        """location:   location of the array in (lat, long, [elev])"""
        ephem.Observer.__init__(self)
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

class AntennaArray(ArrayLocation):
    """A representation of a collection of antennas, their spacings, and
       information about when and where the array is."""
    def __init__(self, antennas, location):
        """antennas:   a list of Antenna instances
        location:   location of the array in (lat, long, [elev])"""
        ArrayLocation.__init__(self, location)
        self.update_antennas(antennas)
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
        self.baseline_order = {}
        for i in range(self.n_ants):
            for j in range(i, self.n_ants):
                bls.append(antennas[j] - antennas[i])
                bl = self.ij2bl(i, j)
                self.baseline_order[bl] = len(bls) - 1
        self.baselines = numpy.array(bls)
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

#  ___ _         _       _                        _                     
# | _ \ |_  ___ /_\  _ _| |_ ___ _ _  _ _  __ _  /_\  _ _ _ _ __ _ _  _ 
# |  _/ ' \(_-</ _ \| ' \  _/ -_) ' \| ' \/ _` |/ _ \| '_| '_/ _` | || |
# |_| |_||_/__/_/ \_\_||_\__\___|_||_|_||_\__,_/_/ \_\_| |_| \__,_|\_, |
#                                                                  |__/ 

class PhsAntennaArray(AntennaArray):
    """A class which adds phasing functionality to AntennaArray."""
    def __init__(self, antennas, location):
        AntennaArray.__init__(self, antennas, location)
        self.update_antennas(antennas)
    def update_antennas(self, antennas):
        """Initialize the antenna array using a list of antennas.  Generates
        zenith baselines and relative delays."""
        AntennaArray.update_antennas(self, antennas)
        dlys = []
        for i in range(self.n_ants):
            for j in range(i, self.n_ants):
                dlys.append(antennas[j].delay - antennas[i].delay)
        self.delays = numpy.array(dlys)
    def get_delay(self, i, j=None):
        """Return the delay corresponding to i,j (see ij2bl for details)."""
        bl = self.ij2bl(i, j)
        return self.delays[self.baseline_order[bl]]
    def phs2src(self, data, src, i, j=None):
        bl = self.ij2bl(i, j)
        w = self.get_projected_baseline(bl, body=src)[2]
        t = self.get_delay(bl)
        phs = numpy.exp(-2*numpy.pi*1j * (w+t) * src.freqs)
        return data * phs
    def unphs2src(self, data, src, i, j=None):
        bl = self.ij2bl(i, j)
        w = self.get_projected_baseline(bl, body=src)[2]
        t = self.get_delay(bl)
        phs = numpy.exp(-2*numpy.pi*1j * (w+t) * src.freqs)
        return data * numpy.conjugate(phs)
