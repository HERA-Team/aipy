"""
Module for adding data simulation support to AntennaArrays.

Author: Aaron Parsons
Date: 11/07/2006
Revisions:
    12/05/2006  arp Conjugated sim_data (-2pi1j, instead of 2pi1j).
    01/01/2007  arp Added gain information for antennas.  Expanded sim_data to 
                    have optional compute of gradient.
    03/02/2007  arp Changed sim_data to do 1 baseline at a time.  Substantial
                    restructuring of parameter passing.  More documentation.
    05/15/2007  arp Split part of SimAntennas into PhsAntennas to have more
                    streamlined support for phasing antenna data.
"""

import ants, numpy

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

#  ____  _              _          _                         
# / ___|(_)_ __ ___    / \   _ __ | |_ ___ _ __  _ __   __ _ 
# \___ \| | '_ ` _ \  / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  ___) | | | | | | |/ ___ \| | | | ||  __/ | | | | | | (_| |
# |____/|_|_| |_| |_/_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|
                                                           
# A default passband fit from PAPER's measured receiver gain
# (in Receiver_Gain.txt)
DFLT_GAIN_POLY = [-8.25e1, 8.34e1, -3.50e1, 7.79e1, -9.71e-1, 6.41e-2, -1.75e-2]

class SimAntenna(ants.Antenna):
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
        ants.Antenna.__init__(self, x, y, z, delay=delay)
        self.beam = beam
        self.offset = offset
        self.gain_poly = gain_poly
        self.select_chans(active_chans)
        self.update_pointing(pointing)
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

#  ___ _         _       _                        _                     
# / __(_)_ __   /_\  _ _| |_ ___ _ _  _ _  __ _  /_\  _ _ _ _ __ _ _  _ 
# \__ \ | '  \ / _ \| ' \  _/ -_) ' \| ' \/ _` |/ _ \| '_| '_/ _` | || |
# |___/_|_|_|_/_/ \_\_||_\__\___|_||_|_||_\__,_/_/ \_\_| |_| \__,_|\_, |
#                                                                  |__/ 

class SimAntennaArray(ants.PhsAntennaArray):
    """A class which adds simulation functionality to AntennaArray."""
    def __init__(self, simantennas, location, active_chans=None):
        """simantennas:     a list of SimAntenna instances
        location:     location of the array in (lat, long, [elev])
        active_chans: channels to be selected for future freq calculations"""
        ants.PhsAntennaArray.__init__(self, simantennas, location)
        self.select_chans(active_chans)
    def select_chans(self, active_chans):
        for a in self.antennas: a.select_chans(active_chans)
        self.freqs = self.antennas[0].beam.active_freqs
    def illuminate(self, ant, srcs, pol=1):
        """Find the degree to which each source in the list 'srcs' is
        illuminated by the beam pattern of 'ant'.  Useful for creating
        simulation data."""
        a = self.antennas[ant]
        nchan = self.freqs.size
        # <GAS> -> Gain * Antenna beam * Source flux
        GAS_sf = numpy.zeros((len(srcs), nchan), dtype=numpy.float)
        for n, s in enumerate(srcs):
            s.compute(self)
            # Skip if source is below horizon
            if s.alt < 0: continue
            GAS_sf[n] = a.response((s.az, s.alt), pol=pol) * s.emission(self)
        return GAS_sf
    def sim_data(self, srcs, ant1, ant2=None, calc_grad=False, stokes=-5):
        r"""Calculates visibilities at a given time for a list of RadioBodys,
        for an array of frequencies according to the Measurement Equation:
            V_{ij}(\nu,t) = \phi_{ij\nu}(t) + \Sigma_n{g_i(\nu) g_j^*(\nu)
                            A_{i\nu}(\hat S_n(t)) A_{j\nu}(\hat S_n(t))
                            S_n\left(\nu\over\nu_0\right)^{\alpha_n}
                            e^{2\pi\vec b_{ij}\cdot\hat S_n(t)
                            + 2\pi\nu\tau_{ij}}}"""
        bl = self.ij2bl(ant1, ant2)
        i, j = self.bl2ij(bl)
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

#  ____  _                 _       _             
# / ___|(_)_ __ ___  _   _| | __ _| |_ ___  _ __ 
# \___ \| | '_ ` _ \| | | | |/ _` | __/ _ \| '__|
#  ___) | | | | | | | |_| | | (_| | || (_) | |   
# |____/|_|_| |_| |_|\__,_|_|\__,_|\__\___/|_|   

class Simulator(SimAntennaArray, ants.SourceList):
    """Contains all information for simulating interferometric data from a
    list of sources."""
    def __init__(self, antennas, location, src_dict, active_chans=None):
        SimAntennaArray.__init__(self, antennas, location,
            active_chans=active_chans)
        ants.SourceList.__init__(self, src_dict, active_chans=active_chans)
        self.set_activity()
        self.chans = self.antennas[0].beam.chans
    def select_chans(self, active_chans):
        """Choose only 'active_chans' for future freq calculations."""
        try: SimAntennaArray.select_chans(self, active_chans)
        except(AttributeError): pass
        try: ants.SourceList.select_chans(self, active_chans)
        except(AttributeError): pass
    def set_activity(self, antennas=[], baselines=[], stokes=[],
            sources=[]):
        """Sets the activity of baselines and stokes parameters.  This
        determines the behavior of 'is_active', which can be used as a tool
        for simulating only the data you are interested in.  You can also
        select which sources are used in 'sim_data'."""
        # Generate all active baselines from provided baselines and antennas
        antennas.sort()
        for i, a1 in enumerate(antennas):
            for a2 in antennas[i:]:
                print a1, a2
                baselines.append(self.ij2bl(a1, a2))
        if len(baselines) == 0: baselines = self.baseline_order.keys()
        self.active_baselines = {}
        for bl in self.baseline_order.keys(): self.active_baselines[bl] = False
        for bl in baselines: self.active_baselines[bl] = True
        # Generate active stokes parameters
        if len(stokes) == 0: stokes = [-5, -6, -7, -8]
        self.active_stokes = {-5:False, -6:False, -7:False, -8:False}
        for s in stokes: self.active_stokes[s] = True
        # Generate active sources used in sim_data
        if len(sources) == 0: sources = self.sources
        self.active_sources = sources
    def is_active(self, bl, stokes):
        """Returns whether the given bl, stokes parameter is active for
        simulation, as dictated by 'set_activity'."""
        return self.active_baselines[bl] and self.active_stokes[stokes]
    def sim_data(self, ant1, ant2=None, calc_grad=False, stokes=-5):
        """Use the active sources defined by 'set_activity' to create sim
        data for the given baseline (=ant1, or (ant1, ant2) if ant2 is
        provided) and stokes parameter."""
        return SimAntennaArray.sim_data(self, self.active_sources, ant1,
            ant2=ant2, calc_grad=calc_grad, stokes=stokes)

