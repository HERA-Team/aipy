"""
Module for adding data simulation support to AntennaArrays.  For most classes,
this means adding gain/amplitude information (as a function of frequency).

Author: Aaron Parsons
Date: 11/07/06
Revisions:
    12/05/06    arp Conjugated sim_data (-2pi1j, instead of 2pi1j).
    01/01/07    arp Added gain information for antennas.  Expanded sim_data to 
                    have optional compute of gradient.
    03/02/07    arp Changed sim_data to do 1 baseline at a time.  Substantial
                    restructuring of parameter passing.  More documentation.
    05/15/07    arp Split part of Antennas into PhsAntennas to have more
                    streamlined support for phasing antenna data.
    10/10/07    arp Switched bandpass parameterization from polynomial to
                    splines.
"""

import ant, numpy, ephem
from scipy.interpolate import splrep, splev

#  _   _ _   _ _ _ _           _____                 _   _                 
# | | | | |_(_) (_) |_ _   _  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | | | | __| | | | __| | | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |_| | |_| | | | |_| |_| | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  \___/ \__|_|_|_|\__|\__, | |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                      |___/        

def splinefit(x, y, smooth=1e-4, order=3):
    """Return a spline interpolation of y as a function of x."""
    return splrep(x, y, s=smooth, k=order)

def splineval(spline, x):
    """Return the values of spline evaluated at x."""
    return splev(x, spline)

#  ____           _ _       ____            _       
# |  _ \ __ _  __| (_) ___ | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \|  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|____/ \___/ \__,_|\__, |
#                                             |___/ 

class RadioBody:
    """A class redefining ephem's sense of brightness for radio astronomy."""
    def __init__(self, strength, meas_freq=.150, spec_index=-1, 
            ang_size=0.):
        """strength:     source flux measured at 'meas_freq'
        meas_freq:    frequency (in GHz) where 'strength' was measured
        spec_index:   index of power-law spectral model of source emission"""
        self._strength = strength
        self._meas_freq = meas_freq
        self._spec_index = spec_index
        self._ang_size = ang_size
    def gen_emission(self, observer):
        """Extrapolate source emission from measured freq using power law."""
        afreqs = observer.ants[0].beam.afreqs
        emission = (afreqs / self._meas_freq)**self._spec_index
        return emission * self._strength

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(ant.RadioFixedBody, RadioBody):
    """A class adding simulation capability to ant.RadioFixedBody"""
    def __init__(self, ra, dec, strength, f_c=.012, meas_freq=.150, 
            spec_index=-1., ang_size=0., name='', **kwargs):
        """ra:           source's right ascension
        dec:          source's declination
        strength:     source flux measured at 'meas_freq'
        meas_freq:    frequency (in GHz) where 'strength' was measured
        spec_index:   index of power-law spectral model of source emission"""
        ant.RadioFixedBody.__init__(self, ra, dec, f_c=f_c, name=name)
        RadioBody.__init__(self, strength, meas_freq=meas_freq,
            spec_index=spec_index, ang_size=ang_size)
    def compute(self, observer, refract=False):
        if self.in_cache(observer, refract): return
        ant.RadioFixedBody.compute(self, observer, refract=refract)
        self.emission = self.gen_emission(observer)
        self.cache(observer, refract)

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(ant.RadioSpecial, RadioBody):
    """A class adding simulation capability to ant.RadioSun"""
    def __init__(self, name, strength, f_c=.012, meas_freq=.150, spec_index=-1.,
            ang_size=8.7e-3, **kwargs):
        """strength:     source flux measured at 'meas_freq'
        meas_freq:    frequency (in GHz) where 'strength' was measured
        spec_index:   index of power-law spectral model of source emission"""
        ant.RadioSpecial.__init__(self, name, f_c=f_c)
        RadioBody.__init__(self, strength, meas_freq=meas_freq,
            spec_index=spec_index, ang_size=ang_size)
    def compute(self, observer, refract=False):
        if self.in_cache(observer, refract): return
        ant.RadioSpecial.compute(self, observer, refract=refract)
        self.emission = self.gen_emission(observer)
        self.cache(observer, refract)

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam(ant.Beam):
    def response(self, zang, az, pol=1):
        """Return the beam response across the band for input zenith angle
        (zang), and azimuth (az).  Rotate beam model 90 degrees if pol == 2."""
        return numpy.ones_like(self.afreqs)

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna(ant.Antenna):
    """A representation of the physical location and beam pattern of an
    individual antenna in an array."""
    def __init__(self, x, y, z, beam, delay=0., offset=0.,
            spline=None, amp=1, pointing=(0.,numpy.pi/2), **kwargs):
        """x, y, z:    Antenna coordinates in equatorial (ns) coordinates
        beam:       Object with function 'response(zang, az)'
        delay:      Cable/systematic delay in ns
        offset:     Frequency-independent phase offset
        spline:     Spline fit of passband.  Default is flat.
        pointing:   Antenna pointing=(az, alt).  Default is zenith"""
        ant.Antenna.__init__(self, x,y,z, beam=beam, delay=delay, offset=offset)
        # Implement a flat = 1 passband if no spline is provided.
        if spline is None:
            spline = (numpy.array([ 0.,  0.,  0.,  0., 1., 1.,  1.,  1.]),
                numpy.array([ 1.,  1.,  1.,  1.]), 3)
        self.update_gain(spline, amp)
        self.update_pointing(pointing)
    def select_chans(self, active_chans=None):
        ant.Antenna.select_chans(self, active_chans)
        self.update_gain(spline=self.spline)
    def update_gain(self, spline=None, amp=None):
        if not spline is None: self.spline = spline
        if not amp is None: self.amp = amp
        self.gain = self.amp * splineval(self.spline, self.beam.afreqs)
    def update_pointing(self, azalt):
        """Set the antenna beam to point at azalt=(az, alt)."""
        self.pointing = azalt
    def response(self, azalt, pol=1):
        """Return the total antenna response to a source at azalt=(az, alt),
        including beam response, per-frequency gain, and a phase offset."""
        zang = ephem.separation(self.pointing, azalt)
        beam_resp = self.beam.response(zang, azalt[0], pol=pol)
        return beam_resp * self.gain

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(ant.AntennaArray):
    """A class which adds simulation functionality to AntennaArray."""
    def __init__(self, simants, location, active_chans=None, **kwargs):
        """simants:         a list of Antenna instances
        location:     location of the array in (lat, long, [elev])
        active_chans: channels to be selected for future freq calculations"""
        ant.AntennaArray.__init__(self, ants=simants, location=location)
        self.select_chans(active_chans)
    def select_chans(self, active_chans=None):
        for a in self.ants: a.select_chans(active_chans)
    def illuminate(self, ant, srcs, pol=1):
        """Find the degree to which each source in the list 'srcs' is
        illuminated by the beam pattern of 'ant'.  Useful for creating
        simulation data."""
        a = self.ants[ant]
        nchan = a.beam.afreqs.size
        # <GAS> -> Gain * Antenna beam * Source flux
        GAS_sf = numpy.zeros((len(srcs), nchan), dtype=numpy.float)
        for n, s in enumerate(srcs):
            # Skip if source is below horizon
            if s.alt < 0: continue
            GAS_sf[n] = a.response((s.az, s.alt), pol=pol) * s.emission
        return GAS_sf
    def sim_data(self, srcs, i, j=None, stokes=-5):
        r"""Calculates visibilities at a given time for a list of RadioBodys,
        for an array of frequencies according to the Measurement Equation:
            V_{ij}(\nu,t) = \phi_{ij\nu}(t) + \Sigma_n{g_i(\nu) g_j^*(\nu)
                            A_{i\nu}(\hat S_n(t)) A_{j\nu}(\hat S_n(t))
                            S_n\left(\nu\over\nu_0\right)^{\alpha_n}
                            e^{2\pi\vec b_{ij}\cdot\hat S_n(t)
                            + 2\pi\nu\tau_{ij}}}"""
        if j is None: i, j = self.bl2ij(i)
        bl = self.ij2bl(i, j)
        afreqs = self.ants[0].beam.afreqs
        if   stokes == -5: pol1, pol2 = 1, 1
        elif stokes == -6: pol1, pol2 = 2, 2
        elif stokes == -7: pol1, pol2 = 1, 2
        elif stokes == -8: pol1, pol2 = 2, 1
        else:
            raise ValueError('Unsupported stokes value: %d not in  [-5,-8].' \
                % (stokes))
        # <GBS> -> Gain * Baseline beam * Source flux
        # <sf> -> source, freq (matrix axes)
        GBS_sf = numpy.conjugate(self.illuminate(i, srcs, pol=pol1))
        GBS_sf *= self.illuminate(j, srcs, pol=pol2)
        # <P> -> phase -> exp(1j*(2*numpy.pi * (z+t) * freq + offset))
        P_sf = numpy.zeros(GBS_sf.shape, dtype=numpy.complex)
        for n, s in enumerate(srcs):
            # If PointingError, P_sf[s] = 0, so flux from source will be nulled
            try: phs, uvw = self.gen_phs(s, bl, with_coord=True)
            except(ant.PointingError): continue
            u, v, w = uvw[:,0], uvw[:,1], uvw[:,2]
            P_sf[n] = numpy.conjugate(phs)
            # Take into account effects of resolving source
            GBS_sf[n] *= numpy.sinc(s._ang_size * numpy.sqrt(u**2+v**2))
        # <GBSP> -> total per-source visibility calculation
        GBSP_sf = GBS_sf * P_sf
        # <V> -> Visibility data -> sum of GBSP over sources
        V_f = GBSP_sf.sum(axis=0)
        return V_f

