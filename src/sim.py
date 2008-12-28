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
    03/06/08    arp Vectorized simulation w/ 4x increase in speed and
                    support for pixel-based simulation.
"""

import ant, numpy as n, ephem, coord, healpix
from interp import interpolate

#  ____           _ _       ____            _       
# |  _ \ __ _  __| (_) ___ | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \|  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|____/ \___/ \__,_|\__, |
#                                             |___/ 

class RadioBody:
    """A class redefining ephem's sense of brightness for radio astronomy."""
    def __init__(self, janskies, mfreq=.150, index=-1):
        """janskies:     source strength
        mfreq:    frequency (in GHz) where strength was measured
        index:   index of power-law spectral model of source emission"""
        self.janskies = janskies
        self.mfreq = mfreq
        self.index = index

#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(ant.RadioFixedBody, RadioBody):
    """A class adding simulation capability to ant.RadioFixedBody"""
    def __init__(self, ra, dec, janskies, mfreq=.150, 
            index=-1., name='', **kwargs):
        ant.RadioFixedBody.__init__(self, ra, dec, name=name)
        RadioBody.__init__(self, janskies, mfreq=mfreq, index=index)

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(ant.RadioSpecial, RadioBody):
    """A class adding simulation capability to ant.RadioSun"""
    def __init__(self, name, janskies, mfreq=.150, index=-1., **kwargs):
        ant.RadioSpecial.__init__(self, name)
        RadioBody.__init__(self, janskies, mfreq=mfreq, index=index)

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/ 

class SrcCatalog(ant.SrcCatalog):
    """A class adding simulation capability to SrcCatalog"""
    def get_vecs(self):
        """Return the 4 arrays needed by sim(): the 
        equatorial vector locations, fluxes, spectral indices, and the 
        frequencies at which the fluxes were measured."""
        srcs = self.values()
        s_eqs = n.array([s.eq for s in srcs]).transpose()
        fluxes = n.array([s.janskies for s in srcs])
        indices = n.array([s.index for s in srcs])
        mfreqs = n.array([s.mfreq for s in srcs])
        return s_eqs, fluxes, indices, mfreqs


#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam2DGaussian(ant.Beam):
    """A 2D Gaussian beam pattern, with default setting for a flat beam."""
    def __init__(self, freqs, xwidth=n.Inf, ywidth=n.Inf):
        ant.Beam.__init__(self, freqs)
        self.update(xwidth=xwidth, ywidth=ywidth)
    def update(self, xwidth=None, ywidth=None):
        """Set the width in the x and y directions of the gaussian beam."""
        if not xwidth is None: self.xwidth = xwidth
        if not ywidth is None: self.ywidth = ywidth
    def response(self, xyz):
        """Return the beam response across the active band for the specified
        topocentric coordinates (with z = up, x = east). 2nd axis should be 
        multiple coordinates.  Returns 'x' pol (rotate pi/2 for 'y')."""
        x,y,z = xyz
        x,y = n.arcsin(x)/self.xwidth, n.arcsin(y)/self.ywidth
        resp = n.sqrt(n.exp(-(x**2 + y**2)))
        resp = n.resize(resp, (self.afreqs.size, resp.size))
        return resp

class BeamPolynomial(ant.Beam):
    """A Beam model that uses a 2D polynomial in cos(2*n*az) for first axis,
    and in freq**n for second axis."""
    def __init__(self, freqs, poly_azfreq=n.array([[.5]])):
        self.poly = poly_azfreq
        ant.Beam.__init__(self, freqs)
        self.update(poly_azfreq)
    def select_chans(self, active_chans):
        ant.Beam.select_chans(self, active_chans)
        self.update()
    def update(self, poly_azfreq=None):
        """Set the width in the x and y directions of the gaussian beam."""
        if poly_azfreq is None: poly_azfreq = self.poly
        elif len(poly_azfreq.shape) == 1: poly_azfreq.shape = self.poly.shape
        self.poly = poly_azfreq
        f = n.resize(self.afreqs, (self.poly.shape[1], self.afreqs.size))
        f = f**n.array([range(self.poly.shape[1])]).transpose()
        self.sigma = n.dot(self.poly, f)
    def response(self, top):
        """Return the beam response across the active band for the specified
        topocentric coordinates (with z = up, x = east). 2nd axis should be 
        multiple coordinates.  Returns 'x' pol (rotate pi/2 for 'y')."""
        az,alt = coord.top2azalt(top)
        zang = n.pi/2 - alt
        if zang.size == 1:
            zang = n.array([zang]); zang.shape = (1,)
            az = n.array([az]); az.shape = (1,)
        a = 2 * n.arange(self.poly.shape[0], dtype=n.float)
        a.shape = (1,) + a.shape; az.shape += (1,); zang.shape += (1,)
        a = n.cos(n.dot(az, a))
        a[:,0] = 0.5
        s = n.dot(a, self.sigma)
        return n.sqrt(n.exp(-(zang/s)**2)).transpose()

class BeamCosSeries(ant.Beam):
    def __init__(self, freqs, poly_cos=n.array([[0., 1.]]), 
            poly_wid=n.array([0., 1.])):
        self.poly_cos = poly_cos
        self.poly_wid = poly_wid
        ant.Beam.__init__(self, freqs)
        self.update(poly_cos, poly_wid)
    def select_chans(self, active_chans):
        ant.Beam.select_chans(self, active_chans)
        self.update()
    def update(self, poly_cos=None, poly_wid=None):
        """Set the width in the x and y directions of the gaussian beam."""
        if poly_cos is None: poly_cos = self.poly_cos
        elif len(poly_cos.shape) == 1: poly_cos.shape = self.poly_cos.shape
        if poly_wid is None: poly_wid = self.poly_wid
        self.poly_cos = poly_cos
        self.poly_wid = poly_wid
    def response(self, top):
        az,alt = coord.top2azalt(top)
        zang = n.pi/2 - alt
        if zang.size == 1:
            zang = n.array([zang]); zang.shape = (1,)
            az = n.array([az]); az.shape = (1,)
        wid = 2 * n.arange(self.poly_wid.shape[0], dtype=n.float)
        wid.shape = (1,) + wid.shape
        az.shape += (1,)
        wid = n.dot(n.cos(n.dot(az, wid)), self.poly_wid)
        x = n.cos(wid * zang)**2
        a = 2 * n.arange(self.poly_cos.shape[0], dtype=n.float)
        a.shape = (1,) + a.shape; zang.shape += (1,)
        p = n.dot(n.cos(n.dot(az, a)), self.poly_cos)
        rv = n.polyval(p.transpose(), x.transpose())
        rv.shape = (1,) + rv.shape
        return rv.clip(0, n.Inf)

class BeamAlm(ant.Beam):
    def __init__(self, freqs, lmax=10, mmax=10, coeffs=None, nside=32):
        self.alm = healpix.Alm(lmax,mmax)
        self.hmap = healpix.HealpixMap(nside, ordering='RING')
        ant.Beam.__init__(self, freqs)
        self.update(coeffs)
    def select_chans(self, active_chans):
        ant.Beam.select_chans(self, active_chans)
        self.update()
    def update(self, coeffs=None):
        if coeffs is None: coeffs = self.alm.get_data()
        self.alm.set_data(coeffs)
        self.hmap.from_alm(self.alm)
    def response(self, top):
        #x,y,z = top
        #top = (-n.abs(x), -n.abs(y), z)
        rv = self.hmap[n.array(top).transpose()]
        rv.shape = (1,) + rv.shape
        return rv.clip(0, n.Inf)
    

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
        beam:       Object with function 'response(xyz=None, azalt=None)'
        delay:      Cable/systematic delay in ns
        offset:     Frequency-independent phase offset
        bp:         Decimated sampling of passband. Default is flat.
        pointing:   Antenna pointing=(az, alt).  Default is zenith"""
        ant.Antenna.__init__(self, x,y,z, beam=beam, delay=delay, offset=offset)
        # Implement a flat passband of ones if no bp is provided
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
        rot = coord.rot_m(twist, z)
        rot = n.dot(rot, coord.rot_m(alt-n.pi/2, y))
        rot = n.dot(rot, coord.rot_m(-az, z))
        self.rot_pol_x = rot
        self.rot_pol_y = n.dot(coord.rot_m(-n.pi/2, z), rot)
    def response(self, top, pol='x'):
        """Return the total antenna response to the specified topocentric 
        coordinates (with z = up, x = east).  This includes beam response and
        per-frequency gain.  1st axis should be xyz, 2nd axis should be 
        multiple coordinates."""
        top = {'x':top, 'y':n.dot(self.rot_pol_y, top)}[pol]
        beam_resp = self.beam.response(top)
        gain = self.gain
        if len(beam_resp.shape) == 2: gain = n.reshape(gain, (gain.size, 1))
        return beam_resp * gain

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(ant.AntennaArray):
    """A class which adds simulation functionality to AntennaArray."""
    def sim(self, i, j, s_eqs, fluxes, 
            indices=0., mfreqs=n.array([.150]), pol='xx'):
        """Simulate visibilites for the (i,j) baseline based on source
        locations (in equatorial coordinates), fluxes, spectral indices,
        the frequencies at which fluxes were measured, and the polarization."""
        assert(pol in ('xx','yy','xy','yx'))
        # Get topocentric coordinates of all srcs
        m = coord.eq2top_m(-self.sidereal_time(), self.lat)
        src_top = n.dot(m, s_eqs)
        # Throw out everthing that is below the horizon
        valid = n.logical_and(src_top[2,:] > 0, fluxes > 0)
        if n.all(valid == 0): return n.zeros_like(self.ants[i].beam.afreqs)
        fluxes = fluxes.compress(valid)
        indices = indices.compress(valid)
        mfreqs = mfreqs.compress(valid)
        src_top = src_top.compress(valid, axis=1)
        s_eqs = s_eqs.compress(valid, axis=1)
        # Get antenna source-dependent gains
        GAi_sf = self.ants[i].response(src_top, pol=pol[0]).transpose()
        GAj_sf = self.ants[j].response(src_top, pol=pol[1]).transpose()
        # Get src fluxes vs. freq
        fluxes.shape = (fluxes.size, 1)
        mfreqs.shape = (mfreqs.size, 1)
        indices.shape = (indices.size, 1)
        freqs = n.resize(self.ants[i].beam.afreqs, 
            (fluxes.size, self.ants[i].beam.afreqs.size))
        I_sf = fluxes * (freqs / mfreqs)**indices
        # Get the phase of each src vs. freq
        E_sf = n.conjugate(self.gen_phs(s_eqs.transpose(), i, j))
        # Combine and sum over sources
        GBIE_sf = GAi_sf * GAj_sf * I_sf**2 * E_sf
        Vij_f = GBIE_sf.sum(axis=0)
        return Vij_f
