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
    """Convert Julian date to an ephem date, which is measured from noon,
    Jan. 1, 1900."""
    return ephem.date(num - 2415020)

def ephem2juldate(num):
    """Convert ephem date (measured from noon, Jan. 1, 1900) to a Julian 
    date."""
    return float(num + 2415020)

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
    def compute(self, observer):
        # Generate a map for projecting baselines to uvw coordinates
        self.map = coord.eq2top_m(observer.sidereal_time()-self.ra, self.dec)
        # Generate a vector coordinates pointing at src
        self.e_vec = coord.radec2eq((self.ra, self.dec))
        self.t_vec = coord.azalt2top((self.az, self.alt))
    def get_crd(self, crdsys, ncrd=3):
        """Return the coordinates of this location in the desired coordinate
        system ('eq','top') in the current epoch.  If ncrd=2, angular
        coordinates (ra/dec,lat/long,az/alt) are returned, and if ncrd=3,
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
    """A class combining ephem's FixedBody with a RadioBody."""
    def __init__(self, ra, dec, name='', **kwargs):
        """ra:           source's right ascension (epoch=2000)
        dec:          source's declination (epoch=2000)"""
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
    """A class combining ephem's Sun, Moon, planets, etc. with a RadioBody."""
    def __init__(self, name, **kwargs):
        RadioBody.__init__(self, name=name)
        self.Body = eval('ephem.%s()' % name)
    def __getattr__(self, nm):
        try: return object.__getattr__(self, nm)
        except(AttributeError): return self.Body.__getattribute__(nm)
    def __setattr__(self, nm, val):
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
    """A class for holding celestial sources."""
    def __init__(self, srcs, **kwargs):
        dict.__init__(self)
        for s in srcs: self.add_src(s)
    def add_src(self, src):
        self[src.src_name] = src
    def get_srcs(self, *args):
        return [self[s] for s in args]
    def compute(self, observer):
        for s in self: self[s].compute(observer)
    def get_crds(self, crdsys, ncrd=3):
        crds = n.array([s.get_crd(crdsys, ncrd=ncrd) for s in self.values()])
        return crds.transpose()

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
        if active_chans is None: active_chans = n.arange(self.freqs.size)
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
    def __init__(self, x, y, z, beam, delay=0., **kwargs):
        """x, y, z:    Antenna coordinates in equatorial (ns) coordinates"""
        self.pos = n.array((x,y,z), n.float64) # must be float64 for mir
        self.beam = beam
        self.delay = delay
    def select_chans(self, active_chans=None):
        self.beam.select_chans(active_chans)
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
    """Collected information about where and when an array is."""
    def __init__(self, location):
        """location:   location of the array in (lat, long, [elev])"""
        ephem.Observer.__init__(self)
        self.pressure = 0
        self.update_location(location)
    def update_location(self, location):
        """Initialize the antenna array for the provided location.  Locations
        may be (lat, long) or (lat, long, elev)."""
        if len(location) == 2: self.lat, self.long = location
        else: self.lat, self.long, self.elev = location
    def set_jultime(self, t=None):
        """Set the current time to a Julian date."""
        if t is None: t = ephem.julian_date()
        self.set_ephemtime(juldate2ephem(t))
    def set_ephemtime(self, t=None):
        """Set the current time to a time derived from the ephem package."""
        if t is None: t = ephem.now()
        self.date, self.epoch = t, t

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(ArrayLocation):
    """A representation of a collection of antennas, their spacings, and
       information about when and where the array is."""
    def __init__(self, location, ants, **kwargs):
        """ location:   location of the array in (lat, long, [elev])
        ants:       a list of Antenna instances"""
        ArrayLocation.__init__(self, location=location)
        self.update_antennas(ants)
        self.select_chans()
    def update_antennas(self, ants):
        """Update the antenna array to use the provided list of antennas."""
        self.ants = ants
        self.update()
    def update(self):
        """If antenna parameters or active channels have been changed, this
        function updates variables derived from them.  Increments the cache
        value to force recalculation of cached values in source projections."""
        bls,dlys,offs = [],[],[] 
        self.bl_order = {}
        for i, ai in enumerate(self.ants):
            for j, aj in enumerate(self.ants[i:]):
                bl = self.ij2bl(i, j+i)
                bls.append(aj - ai)
                dlys.append(aj.delay - ai.delay)
                self.bl_order[bl] = len(bls) - 1
        self.bls,self.dlys,self.offs = n.array(bls),n.array(dlys),n.array(offs)
        # Compute (static) zenith baselines
        m = coord.eq2top_m(0., self.lat)
        self.zbls = n.dot(m, self.bls.transpose()).transpose()
    def set_ephemtime(self, t=None):
        ArrayLocation.set_ephemtime(self, t=t)
        # Rotate baselines in eq coords to current sidereal time
        m = coord.rot_m(-self.sidereal_time(), n.array([0.,0.,1.]))
        self.ebls = n.dot(m, self.bls.transpose()).transpose()
    def select_chans(self, active_chans=None):
        """Select which channels are used in computations.  Default is all."""
        for a in self.ants: a.select_chans(active_chans)
        self.update()
    def ij2bl(self, i, j):
        """Convert baseline i,j (0 indexed) to Miriad's (i+1) << 8 | (j+1)"""
        return (int(i)+1) << 8 | (int(j)+1)
    def bl2ij(self, bl):
        """Convert Miriad's (i+1) << 8 | (j+1) to i,j (0 indexed)"""
        bl = int(bl)
        return ((bl >> 8) & 255) - 1, (bl & 255) - 1
    def get_baseline(self, i, j, src='z'):
        """Return the baseline corresponding to i,j in various coordinate 
        projections: src='e' for current equatorial, 'z' for zenith 
        topocentric, 'r' for unrotated equatorial, or a RadioBody for
        a projection toward that source."""
        b = self.bl_order[self.ij2bl(i,j)]
        if type(src) == str:
            if src == 'e': return self.ebls[b]
            elif src == 'z': return self.zbls[b]
            elif src == 'r': return self.bls[b]
            else: raise ValueError('Unrecognized source:' + src)
        if src.alt < 0: raise PointingError('%s below horizon' % src.src_name)
        return n.dot(src.map, self.bls[b])
    def get_delay(self, i, j):
        """Return the delay corresponding to i,j."""
        return self.dlys[self.bl_order[self.ij2bl(i,j)]]
    def gen_uvw(self, i, j, src='z'):
        """Compute uvw coordinates for a provided RadioBody, or 'z' for
        zenith uvw coordinates."""
        xyz = self.get_baseline(i,j, src=src)
        afreqs = self.ants[0].beam.afreqs
        if len(xyz.shape) == 1: xyz = n.resize(xyz, (afreqs.size, xyz.size))
        return xyz * n.reshape(afreqs, (afreqs.size, 1))
    def gen_phs(self, src, i, j):
        """Return the phasing to be multiplied to data to point to src."""
        try: src = src.e_vec
        except(AttributeError): pass
        bl_dot_s = n.dot(src, self.get_baseline(i,j,src='e'))
        if len(bl_dot_s.shape) >= 1: bl_dot_s.shape = (bl_dot_s.size, 1)
        t = self.get_delay(i,j)
        afreqs = self.ants[0].beam.afreqs
        afreqs = n.reshape(afreqs, (1,afreqs.size))
        phs = n.exp(-1j*(2*n.pi*n.dot(bl_dot_s + t, afreqs)))
        return phs.squeeze()
    def phs2src(self, data, src, i, j):
        """Apply phasing to zenith data to point to src."""
        return data * self.gen_phs(src, i, j)
    def unphs2src(self, data, src, i, j):
        """Remove phasing from src data to point to zenith."""
        return data * n.conjugate(self.gen_phs(src, i, j))
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
                if norm_bandpass is None: data -= n.ma.average(data)
                else: data -= n.ma.sum(data) * norm_bandpass
            else:
                img = n.fft.ifft(data.filled(0))
                img[swath:-swath] = 0
                data -= n.fft.fft(img)
            data /= phs
        return data
