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
    12/12/07    arp Switched bp again, this time to interpolation from a
                    decimated bandpass.
"""

import ant, numpy as n, ephem, coord
from interp import interpolate

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
    def response(self, topo_xyz):
        """Return the beam response across the active band for the specified
        topocentric coordinates (with z pointing towards zenith, x pointing
        north, and y pointing west).  1st axis should be xyz, 2nd axis should
        be multiple coordinates."""
        return n.ones_like(self.afreqs)

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna(ant.Antenna):
    """A representation of the physical location and beam pattern of an
    individual antenna in an array."""
    def __init__(self, x, y, z, beam, delay=0., offset=0.,
            bp=None, amp=1, pointing=(0.,n.pi/2,0), **kwargs):
        """x, y, z:    Antenna coordinates in equatorial (ns) coordinates
        beam:       Object with function 'response(zang, az)'
        delay:      Cable/systematic delay in ns
        offset:     Frequency-independent phase offset
        bp:         Decimated sampling of passband. Default is flat.
        pointing:   Antenna pointing=(az, alt).  Default is zenith"""
        ant.Antenna.__init__(self, x,y,z, beam=beam, delay=delay, offset=offset)
        # Implement a flat = 1 passband if no bp is provided
        if bp is None: bp = n.ones(beam.freqs.size, dtype=n.float)
        self.update_gain(bp, amp)
        self.update_pointing(*pointing)
    def select_chans(self, active_chans=None):
        ant.Antenna.select_chans(self, active_chans)
        self.update_gain()
    def update_gain(self, bp=None, amp=None):
        if not bp is None: self.bp = bp
        if not amp is None: self.amp = amp
        bp = self.bp
        dec_factor = self.beam.freqs.size / self.bp.size
        if dec_factor != 1: bp = interpolate(self.bp, dec_factor)
        gain = self.amp * bp
        self.gain = gain.take(self.beam.chans)
    def update_pointing(self, az=0, alt=n.pi/2, twist=0):
        """Set the antenna beam to point at (az, alt) with the specified
        right-hand twist to the polarizations.  Polarization y is assumed
        to be +pi/2 azimuth from pol x."""
        y, z = n.array([0,1,0]), n.array([0,0,1])
        rot = coord.rotmatrix(twist, z)
        rot = n.dot(rot, coord.rotmatrix(alt-n.pi/2, y))
        rot = n.dot(rot, coord.rotmatrix(-az, z))
        self.rot_pol_x = rot
        self.rot_pol_y = n.dot(coord.rotmatrix(-n.pi/2, z), rot)
    def response(self, xyz_topo, pol='x'):
        """Return the total antenna response to the specified topocentric 
        coordinates (with z pointing towards zenith, x pointing north, and y 
        pointing west).  This includes including beam response, per-frequency 
        gain, and a phase offset.  1st axis should be xyz, 2nd axis should be 
        multiple coordinates."""
        assert(pol in ['x','y'])
        if pol == 'x': rot = self.rot_pol_x
        else: rot = self.rot_pol_y
        beam_resp = self.beam.response(n.dot(rot, xyz_topo))
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
    def illuminate(self, ant, srcs, pol='x'):
        """Find the degree to which each source in the list 'srcs' is
        illuminated by the beam pattern of 'ant'.  Useful for creating
        simulation data."""
        a = self.ants[ant]
        nchan = a.beam.afreqs.size
        # <GAS> -> Gain * Antenna beam * Source flux
        GAS_sf = n.zeros((len(srcs), nchan), dtype=n.float)
        for n, s in enumerate(srcs):
            # Skip if source is below horizon
            if s.alt < 0: continue
            GAS_sf[n] = a.response((s.az, s.alt), pol=pol) * s.emission
        return GAS_sf
    def sim_data(self, srcs, i, j, pol='xx'):
        r"""Calculates visibilities at a given time for a list of RadioBodys,
        for an array of frequencies according to the Measurement Equation:
            V_{ij}(\nu,t) = \phi_{ij\nu}(t) + \Sigma_n{g_i(\nu) g_j^*(\nu)
                            A_{i\nu}(\hat S_n(t)) A_{j\nu}(\hat S_n(t))
                            S_n\left(\nu\over\nu_0\right)^{\alpha_n}
                            e^{2\pi\vec b_{ij}\cdot\hat S_n(t)
                            + 2\pi\nu\tau_{ij}}}"""
        assert(pol in ['xx', 'yy', 'xy', 'yx'])
        bl = self.ij2bl(i, j)
        afreqs = self.ants[0].beam.afreqs
        pol1, pol2 = pol
        # <GBS> -> Gain * Baseline beam * Source flux
        # <sf> -> source, freq (matrix axes)
        GBS_sf = n.conjugate(self.illuminate(i, srcs, pol=pol1))
        GBS_sf *= self.illuminate(j, srcs, pol=pol2)
        # <P> -> phase -> exp(1j*(2*n.pi * (z+t) * freq + offset))
        P_sf = n.zeros(GBS_sf.shape, dtype=n.complex)
        for n, s in enumerate(srcs):
            # If PointingError, P_sf[s] = 0, so flux from source will be nulled
            try: phs, uvw = self.gen_phs(s, i, j, with_coord=True)
            except(ant.PointingError): continue
            u, v, w = uvw[:,0], uvw[:,1], uvw[:,2]
            P_sf[n] = n.conjugate(phs)
            # Take into account effects of resolving source
            GBS_sf[n] *= n.sinc(s._ang_size * n.sqrt(u**2+v**2))
        # <GBSP> -> total per-source visibility calculation
        GBSP_sf = GBS_sf * P_sf
        # <V> -> Visibility data -> sum of GBSP over sources
        V_f = GBSP_sf.sum(axis=0)
        return V_f

