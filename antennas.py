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
                    Conjugated sim data (-2pi1j, instead of 2pi1j).
    01/01/2007  arp Added gain information for antennas.  Expanded sim_data
                    to have optional compute of gradient.
    01/12/2007  arp Fixed H (hour angle) to be computed in sidereal time, not
                    solar.
    01/28/2007  arp Sync'd to latest miriad.py (0.0.4) version.  Hefty (4-5x)
                    speed up!  Also converted to use numpy-1.0.1
    03/02/2007  arp Changed sim_data to do 1 baseline at a time.  Substantial
                    restructuring of parameter passing.  More documentation.
"""

# Copyright (C) 2007 Aaron Parsons
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

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
            active_chans=None):
        """strength:     source flux measured at 'meas_freq'
        freqs:        frequencies (in GHz) at bin centers across spectrum
        meas_freq:    frequency (in GHz) where 'strength' was measured
        spec_index:   index of power-law spectral model of source emission
        active_chans: channels to be selected for future freq calculations"""
        self.freqs = freqs
        self._meas_freq = meas_freq
        self._strength = strength
        self._spec_index = spec_index
        self.select_chans(active_chans)
        self.prev_sidereal_time = 0     # Used to avoid redundant map calc.
    def select_chans(self, active_chans):
        """Choose only 'active_chans' for future freq calculations."""
        if active_chans is None:
            try: active_chans = self.chans
            except: active_chans = numpy.arange(self.freqs.size)
        self.chans = active_chans
        self.update(self._strength, self._spec_index)
    def update(self, strength, spec_index):
        """Update parameters for source strength and spectral index."""
        #self._strength, self._spec_index = strength, spec_index
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
        t = observer.date - 39051.3774421
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

#  ____                       
# | __ )  ___  __ _ _ __ ___  
# |  _ \ / _ \/ _` | '_ ` _ \ 
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam:
    """Template for an antenna's beam pattern on the sky."""
    def __init__(self, freqs, active_chans=None):
        """freqs:        frequencies (in GHz) at bin centers across spectrum
        active_chans: channels to be selected for future freq calculations"""
        self.freqs = freqs
        self.select_chans(active_chans)
    def select_chans(self, active_chans):
        """Choose only 'active_chans' for future freq calculations."""
        if active_chans is None: active_chans = numpy.arange(self.freqs.size)
        self.chans = active_chans
        self.active_freqs = self.freqs.take(active_chans)
    def response(self, zang, az, pol=1):
        """Return the beam response across the band for input zenith angle 
        (zang), and azimuth (az).  Rotate beam model 90 degrees if pol == 2."""
        return numpy.ones_like(self.active_freqs)

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

# A default passband fit from PAPER's measured receiver gain 
# (in Receiver_Gain.txt)
DFLT_GAIN_POLY = [-8.25e1, 8.34e1, -3.50e1, 7.79e1, -9.71e-1, 6.41e-2, -1.75e-2]

class Antenna:
    """A representation of the physical location and beam pattern of an
    individual antenna in an array."""
    def __init__(self, x, y, z, beam, delay=0., offset=0., 
            gain_poly=DFLT_GAIN_POLY, active_chans=None, 
            pointing=(0.,numpy.pi/2)):
        """x, y, z:    Antenna coordinates in equatorial (ns) coordinates
        beam:       Object with function 'response(zang, az)'
        delay:      Cable/systematic delay in ns
        offset:     Frequency-independent phase offset
        gain_poly:  Polynomial fit of passband
        pointing:   Antenna pointing=(az, alt).  Default is zenith"""
        self.pos = numpy.array((x,y,z), numpy.float64) # must be float64 for c2m
        self.beam = beam
        self.delay = delay
        self.offset = offset
        self.gain_poly = gain_poly
        self.select_chans(active_chans)
        self.update_pointing(pointing)
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
    def select_chans(self, active_chans):
        self.beam.select_chans(active_chans)
        self.update_gain(self.gain_poly)
    def update_gain(self, gain_poly):
        """Set a passband gain based on the polynomial fit 'gain_poly'.  Select
        only 'active_chans' for future freq calculations."""
        try: len(gain_poly)
        except(TypeError): gain_poly = [gain_poly]
        self.gain_poly = gain_poly
        self.gain = numpy.polyval(gain_poly, self.beam.active_freqs)
    def update_pointing(self, azalt):
        """Set the antenna beam to point at azalt=(az, alt)."""
        self.pointing = azalt
    def response(self, azalt, pol=1):
        """Return the total antenna response to a source at azalt=(az, alt),
        including beam response, per-frequency gain, and a phase offset."""
        zang = ephem.separation(self.pointing, azalt)
        beam_resp = self.beam.response(zang, azalt[0], pol=pol)
        offset = numpy.exp(2*math.pi*1j*self.offset)
        return beam_resp * self.gain * offset

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 
class AntennaArray(ephem.Observer):
    """A representation of a collection of antennas, their spacings, and
       information about when and where the array is."""
    def __init__(self, antennas, location):
        """antennas:   a list of Antenna instances
        location:   location of the array in (lat, long, [elev])"""
        ephem.Observer.__init__(self)
        self.update_location(location)
        self.update_antennas(antennas)
    def gen_bl(self, i, j=None):
        """Convert from i,j (counting from 0) baseline notation to
        (i+1) << 8 | (j+1) baseline notation.  If j is not provided, assume
        i is in the desired baseline notation."""
        if j is None: return int(i)
        else: return (int(i)+1) << 8 | (int(j)+1)
    def gen_ij(self, bl):
        """Convert from (i+1) << 8 | (j+1) baseline notation to
        i,j (counting from 0) baseline notation."""
        bl = int(bl)
        return ((bl >> 8) & 255) - 1, (bl & 255) - 1
    def update_antennas(self, antennas):
        """Initialize the antenna array using a list of antennas.  Generates
        zenith baselines and relative delays."""
        self.antennas = antennas
        self.n_ants = len(antennas)
        bls = []; dlys = []
        self.baseline_order = {}
        for i in range(self.n_ants):
            for j in range(i, self.n_ants):
                bls.append(antennas[j] - antennas[i])
                dlys.append(antennas[j].delay - antennas[i].delay)
                bl = self.gen_bl(i, j)
                self.baseline_order[bl] = len(bls) - 1
        self.baselines = numpy.array(bls)
        self.delays = numpy.array(dlys)
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
    def get_baseline(self, i, j=None):
        """Return the baseline corresponding to i,j (see gen_bl for details)."""
        bl = self.gen_bl(i, j)
        return self.baselines[self.baseline_order[bl]]
    def get_delay(self, i, j=None):
        """Return the delay corresponding to i,j (see gen_bl for details)."""
        bl = self.gen_bl(i, j)
        return self.delays[self.baseline_order[bl]]
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
    def set_jultime(self, t=None):
        """Set the current time to a Julian date."""
        if t is None: t = ephem.julian_date()
        self.set_ephemtime(juldate2ephem(t))
    def set_ephemtime(self, t=None):
        """Set the current time to a time drived from the ephem package."""
        if t is None: t = ephem.now()
        self.date = t
    def get_projected_baseline(self, i, j=None, body=None):
        """Rotate zenith (the default pointing) to the actual
        pointing, return the projected baseline i,j."""
        # If no body is provided to point to, return zenith baselines
        if body is None: uvw_map = gen_uvw_map(0., self.lat)
        else: uvw_map = body.gen_uvw_map(self)
        proj_bl = numpy.dot(uvw_map, self.get_baseline(i, j))
        return proj_bl
    def src_eq_vec(self, body):
        """Return an equatorial vector pointing at body."""
        body.compute(self)
        return self.top2eq(azalt2top(body.az, body.alt), no_norm=True)

class SimAntennaArray(AntennaArray):
    """A class which adds simulation functionality to AntennaArray."""
    def __init__(self, antennas, location, active_chans=None):
        """antennas:     a list of Antenna instances
        location:     location of the array in (lat, long, [elev])
        active_chans: channels to be selected for future freq calculations"""
        AntennaArray.__init__(self, antennas, location)
        self.select_chans(active_chans)
        self.prev_date = 0              # Used to avoid redundant illuminations
        self.src_hash = 0               # Used to avoid redundant illuminations
    def select_chans(self, active_chans):
        for a in self.antennas: a.select_chans(active_chans)
        self.freqs = self.antennas[0].beam.active_freqs
    def illuminate(self, ant, srcs, pol=1):
        """Find the degree to which each source in the list 'srcs' is 
        illuminated by the beam pattern of 'ant'.  Useful for creating 
        simulation data."""  
        if self.prev_date != self.date or hash(str(srcs)) != self.src_hash:
            self.prev_illuminations = {1:{}, 2:{}}
            self.prev_date = self.date
            self.src_hash = hash(str(srcs))
        try: return self.prev_illuminations[pol][ant]
        except(KeyError): pass
        a = self.antennas[ant]
        nchan = self.freqs.size
        # <GAS> -> Gain * Antenna beam * Source flux
        GAS_sf = numpy.zeros((len(srcs), nchan), dtype=numpy.float)
        for n, s in enumerate(srcs):
            s.compute(self)
            # Skip if source is below horizon
            if s.alt < 0: continue
            GAS_sf[n] = a.response((s.az, s.alt), pol=pol) * s.emission(self)
        self.prev_illuminations[pol][ant] = GAS_sf
        return GAS_sf
    def sim_data(self, srcs, ant1, ant2=None, calc_grad=False, stokes=-5):
        r"""Calculates visibilities at a given time for a list of RadioBodys,
        for an array of frequencies according to the Measurement Equation:
            V_{ij}(\nu,t) = \phi_{ij\nu}(t) + \Sigma_n{g_i(\nu) g_j^*(\nu) 
                            A_{i\nu}(\hat S_n(t)) A_{j\nu}(\hat S_n(t))
                            S_n\left(\nu\over\nu_0\right)^{\alpha_n}
                            e^{2\pi\vec b_{ij}\cdot\hat S_n(t) 
                            + 2\pi\nu\tau_{ij}}}"""
        bl = self.gen_bl(ant1, ant2)
        i, j = self.gen_ij(bl)
        if   stokes == -5: pol1, pol2 = 1, 1
        elif stokes == -6: pol1, pol2 = 2, 2
        elif stokes == -7: pol1, pol2 = 1, 2
        elif stokes == -8: pol1, pol2 = 2, 1
        else: 
            raise ValueError('Unsupported stokes value: %d not in  [-4,-8].' \
                % (stokes))
        # <GBS> -> Gain * Baseline beam * Source flux
        # <sf> -> source, freq (matrix axes)
        GBS_sf = numpy.conjugate(self.illuminate(i, srcs, pol=pol1))
        GBS_sf *= self.illuminate(j, srcs, pol=pol2)
        # <W> -> w component of baseline
        W_sf = numpy.zeros(GBS_sf.shape, dtype=numpy.float)
        for n, s in enumerate(srcs):
            # Nanosec coords
            z = self.get_projected_baseline(bl, body=s)[2]
            # Actual wavenumber coords
            W_sf[n] = z * self.freqs
        # <T> -> tau -> baseline delay
        T__f = self.get_delay(bl) * self.freqs
        # <E> -> exp(w component + baseline delay)
        E_sf = numpy.exp(2*numpy.pi*1j * (W_sf + T__f))
        # <GBSE> -> total per-source visibility calculation
        GBSE_sf = GBS_sf * E_sf
        # <V> -> Visibility data -> sum of GBSE over sources
        V_f = GBSE_sf.sum(axis=0)
        if calc_grad: return V_f, GBSE_sf
        else: return V_f

#  _____         _   _                     _     
# |_   _|__  ___| |_| |__   ___ _ __   ___| |__  
#   | |/ _ \/ __| __| '_ \ / _ \ '_ \ / __| '_ \ 
#   | |  __/\__ \ |_| |_) |  __/ | | | (__| | | |
#   |_|\___||___/\__|_.__/ \___|_| |_|\___|_| |_|

if __name__ == '__main__':
    import params, miriad, sys, os
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('antennas.py [options] *.uv')
    p.set_description(__doc__)
    p.add_option('-d', '--diff', dest='diff', action='store_true',
        help='Subtract simulated data from measured data.')
    p.add_option('-f', '--fit_file', dest='fit_file', default=None,
        help='File where latest fit is stored.')
    opts, args = p.parse_args(sys.argv[1:])

    def normalize(a):
        """Normalize each number to magnitude 1, but preserve phase info."""
        aa = abs(a)
        return numpy.ma.where(aa > 0, a/aa, 0)

    aa = params.antenna_array
    if not opts.fit_file is None: aa.fromfile(opts.fit_file)
    chans = numpy.arange(params.NCHAN)
    for a in args:
        if opts.diff: outfile = a + '.dif'
        else: outfile = a + '.sim'
        if os.path.exists(outfile):
            print 'Skipping:', a, '(file exists)'
            continue
        else: print 'Working on:', a
        uvi = miriad.UV(a)
        uvo = miriad.UV(outfile, 'new')
        if opts.diff:
            def mfunc(uv, preamble, data, vars):
                u, v, w, t, bl = preamble
                aa.set_jultime(t)
                sim_data = numpy.ma.array(
                    aa.sim_data(params.sources.values(), bl, stokes=uv['pol']), mask=0)
                return preamble, d - sim_data
            hist_lines = 'DIF: Model subtracted from data.\n'
        else:
            def mfunc(uv, preamble, data, vars):
                u, v, w, t, bl = preamble
                aa.set_jultime(t)
                sim_data = numpy.ma.array(
                    aa.sim_data(params.sources.values(), bl, stokes=uv['pol']), mask=0)
                return preamble, sim_data
            hist_lines = 'SIM: Data overwritten with simulation data.\n'
        miriad.map_uv(uvi, uvo, mfunc=mfunc, append2history=hist_lines, 
            send_time_blks=False)
        del(uvi); del(uvo)
