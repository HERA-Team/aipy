'''This is a calibration file for data collected at PAPER in Boolardy
on JD 2454303.'''

import aipy as a, numpy as n

class BeamNoFlaps(ant.Beam):
    """A specific beam model for the PAPER experiment.  This model is for
    a single dipole element with no flaps."""
    def __init__(self, freqs, active_chans=None, **kwargs):
        """The axes of the Cs matrices are polynomials in frequency going
        to the right, and polys in cos(2*az) going down."""
        CsAm = n.array([        # N. Gugliucci 08/07
            [ 2.3541326  ,-0.0039133512 , 1.6055088e-05,-2.7468911e-08],
            [-0.46909345 , 0.0084471178 ,-5.1260711e-05, 1.0299793e-07],
            [ 0.32753617 ,-0.0081176326 , 5.8822952e-05,-1.3303273e-07],
            [-0.046844105, 0.00076223627,-3.5474502e-06, 4.4132035e-09],
            [-0.073523813, 0.0018151892 ,-1.3435102e-05, 3.1225928e-08],
            [ 0.047340855,-0.00085097424, 4.9799933e-06,-9.5164123e-09],
        ])
        CsXo = n.array([        # N. Gugliucci 08/07
           [-121.29224, 1.9851554 ,-0.011876889  , 2.5222526e-05],
           [ 76.969303,-1.3947796 , 0.0085644354 ,-1.7448153e-05],
           [-36.638691, 0.93699466,-0.0068616164 , 1.5544311e-05],
           [ 10.189859,-0.18212180, 0.00098309486,-1.6152395e-06],
           [ 5.9997050,-0.15737420, 0.0012090764 ,-2.8862905e-06],
           [-5.6561847, 0.10468756,-0.00063126068, 1.2444705e-06],
        ])
        CsSd = n.array([        # N. Gugliucci 08/07
            [ 143.84525,-1.1088605  , 0.0048397670 ,-7.1054741e-06],
            [-104.00886, 1.9980993  ,-0.013304344  , 2.8955473e-05],
            [ 28.304230,-0.75088201 , 0.0056338561 ,-1.2898564e-05],
            [-8.7149717, 0.16108215 ,-0.00090283393, 1.5386691e-06],
            [-3.4672940, 0.091929321,-0.00071502397, 1.7311496e-06],
            [ 3.4123240,-0.063083812, 0.00038093617,-7.5356570e-07],
        ])
        mhz_freqs = 1e3 * freqs # GHz -> MHz
        mhz_freqs = n.array([
           n.ones_like(mhz_freqs),
           mhz_freqs,
           mhz_freqs**2,
           mhz_freqs**3,
        ])
        self.BAm = n.dot(CsAm, mhz_freqs)
        self.BXo = n.dot(CsXo, mhz_freqs)
        self.BSd = n.dot(CsSd, mhz_freqs)
        ant.Beam.__init__(self, freqs, active_chans)
    def response(self, xyz):
        """Return the beam response across the active band for the specified
        topocentric coordinates (with z = up, x = east). 2nd axis should be 
        multiple coordinates.  Returns 'x' pol (rotate pi/2 for 'y')."""
        az,alt = coord.top2azalt(xyz)
        zang = n.pi/2 - alt
        if zang.size == 1:
            zang = n.array([zang]); zang.shape = (1,)
            az = n.array([az]); az.shape = (1,)
        a = n.array([0,2,4,6,8,10],dtype=n.float)
        a.shape = (1,) + a.shape; az.shape += (1,); zang.shape += (1,)
        a = n.cos(n.dot(az, a))
        a[:,0] = 0.5
        a1 = n.dot(a, self.BAm_sel)
        a2 = n.dot(a, self.BXo_sel)
        a3 = n.dot(a, self.BSd_sel)
        z = (180*zang/n.pi - a2) / a3
        rv = n.sqrt(a1 * n.exp(-z**2/2))
        return rv.transpose()
    def get_params(self, prm_list):
        return {}
    def set_params(self, prm_list):
        return

prms = {
    'loc': ('-26:44:12.74', '116:39:59.33'), # Boolardy, Australia
    'antpos':
        [[    0.,    0.,    0.],
         [-100.67, 138.14,-198.28],
         [ -68.62, 381.81,-133.88],
         [  59.4, 465.3, 120.0],],
    'delays': [0., -0.355, -0.400, -1.72],
    'offsets': [0., 0, 0., 0],
    'amps': [1.27e-2,1.27e-2,1.27e-2,1.27e-2],
    'passbands': n.array([
        [1,] * 4
    ]),
    'beam': BeamNoFlaps,
}

def get_aa(freqs):
    '''Return the AntennaArray to be used fro simulation.'''
    beam = prms['beam'](freqs)
    try: beam.set_params(prms)
    except(AttributeError): pass
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    assert(len(prms['delays']) == nants and len(prms['offsets']) == nants \
        and len(prms['amps']) == nants and len(prms['passbands']) == nants)
    for pos, dly, off, amp, bp in zip(prms['antpos'], prms['delays'],
            prms['offsets'], prms['amps'], prms['passbands']):
        antennas.append(
            a.fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, bp=bp)
        )
    aa = a.fit.AntennaArray(prms['loc'], antennas)
    aa.set_params(prms['aa'])
    return aa

