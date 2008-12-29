'''This is a calibration file for data collected at PAPER in Boolardy
on JD 2454303.'''

import aipy as a, numpy as n

class AntennaArray(a.fit.AntennaArray):
    def __init__(self, *args, **kwargs):
        self.blpolys = {
            a.miriad.ij2bl(0,1): [-104796371139.4473, 160441814654.44785, -110182538631.84428, 44694468216.905701, -11858579128.368996, 2150321743.3557014, -269862279.81853288, 23143997.116391312, -1298089.2221271475, 42994.932040217027, -638.59608561083041] ,
            a.miriad.ij2bl(1,2): [-63740636591.484711, 96885738955.039749, -66034569398.731644, 26574386662.96114, -6992258572.4814806, 1256832254.7677145, -156281492.34156352, 13273518.018142868, -736908.42666358687, 24146.529490444853, -354.60739558008208] ,
            a.miriad.ij2bl(0,3): [-63682359362.109413, 97105366117.894302, -66403097441.965538, 26814720468.010704, -7080870701.9061899, 1277552414.0722442, -159485576.05700821, 13601911.851185348, -758442.50680062058, 24966.75381468124, -368.43799666599625] ,
            a.miriad.ij2bl(0,2): [-45588496962.026039, 69493557585.217636, -47517797605.718414, 19191536744.091511, -5069824217.0035839, 915280728.20160222, -114357136.33411066, 9763406.0424792115, -545095.19123534171, 17969.814784303493, -265.61684338818304] ,
            a.miriad.ij2bl(2,3): [-75075689969.707001, 115586048317.9888, -79809122012.314484, 32542957939.526997, -8677729578.2758389, 1581065076.2872937, -199324540.45551682, 17168266.634772703, -966847.79879398574, 32146.384893290269, -479.18050079945351] ,
            a.miriad.ij2bl(1,3): [-76848353942.03508, 117642473000.98462, -80752823843.033478, 32728959328.150063, -8673117413.04319, 1570135131.1837056, -196648478.39832252, 16823764.697559636, -940905.72888705588, 31062.451719870805, -459.66636751364967] ,
        }
        a.fit.AntennaArray.__init__(self, *args, **kwargs)
    def update(self):
        a.fit.AntennaArray.update(self)
        afreqs = self.ants[0].beam.afreqs
        self.blgain = {}
        for bl in self.blpolys:
            self.blgain[bl] = n.polyval(self.blpolys[bl], afreqs)
            self.blgain[bl] = n.clip(self.blgain[bl], 0, n.Inf)
    def passband(self, i, j):
        bl = self.ij2bl(i,j)
        return a.fit.AntennaArray.passband(self, i, j) * self.blgain[bl]

class BeamNoFlaps(a.ant.Beam):
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
        a.ant.Beam.__init__(self, freqs, active_chans)
    def response(self, xyz):
        """Return the beam response across the active band for the specified
        topocentric coordinates (with z = up, x = east). 2nd axis should be 
        multiple coordinates.  Returns 'x' pol (rotate pi/2 for 'y')."""
        az,alt = a.coord.top2azalt(xyz)
        zang = n.pi/2 - alt
        if zang.size == 1:
            zang = n.array([zang]); zang.shape = (1,)
            az = n.array([az]); az.shape = (1,)
        A = n.array([0,2,4,6,8,10],dtype=n.float)
        A.shape = (1,) + A.shape; az.shape += (1,); zang.shape += (1,)
        A = n.cos(n.dot(az, A))
        A[:,0] = 0.5
        a1 = n.dot(A, self.BAm_sel)
        a2 = n.dot(A, self.BXo_sel)
        a3 = n.dot(A, self.BSd_sel)
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
        [[   0.,    0.,    0.],
         [-105.39, 139.44,-201.57],
         [ -68.28, 381.97,-134.25],
         [  62.76, 462.75, 119.78]],
    'delays': [0., 5.077, -0.223, -3.564],
    #'delays': [ 0.0, 0.0, 0.0, 0.0],
    #'offsets': [0.0, 0.013,-0.150, -0.323],
    #'offsets': [0.0, 0.0, 0.0, 0.941],
    'offsets': [0.0, 0.039,-0.035, -0.049],
    #'amps': [2.68e-3,2.92e-3,2.68e-3,2.75e-3],
    'amps': [1., 1., 1., 1.],
    'bp_r': n.array([
        [1.],
        [1.],
        [1.],
        [1.],
        #[-9.77733e+05,  6.05955e+05, -1.39787e+05,  1.42159e+04, -5.35700e+02],
        #[-9.77733e+05,  6.05955e+05, -1.39787e+05,  1.42159e+04, -5.35700e+02],
        #[-9.77733e+05,  6.05955e+05, -1.39787e+05,  1.42159e+04, -5.35700e+02],
        #[-9.77733e+05,  6.05955e+05, -1.39787e+05,  1.42159e+04, -5.35700e+02],
    ]),
    'bp_i': n.array([
        [0.000],
        [0.000],
        [0.000],
        [0.000],
    ]),
    #'beam': BeamNoFlaps,
    'beam': a.fit.BeamPolynomial,
    'bm_poly': n.array([        # N. Gugliucci 08/07
        [ 143.84525,-1.1088605  , 0.0048397670 ,-7.1054741e-06],
        [-104.00886, 1.9980993  ,-0.013304344  , 2.8955473e-05],
        [ 28.304230,-0.75088201 , 0.0056338561 ,-1.2898564e-05],
        [-8.7149717, 0.16108215 ,-0.00090283393, 1.5386691e-06],
        [-3.4672940, 0.091929321,-0.00071502397, 1.7311496e-06],
        [ 3.4123240,-0.063083812, 0.00038093617,-7.5356570e-07],
    ]),
}

def get_aa(freqs):
    '''Return the AntennaArray to be used fro simulation.'''
    beam = prms['beam'](freqs)
    try: beam.set_params(prms)
    except(AttributeError): pass
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    assert(len(prms['delays']) == nants and len(prms['amps']) == nants \
        and len(prms['bp_r']) == nants and len(prms['bp_i']) == nants)
    for pos, dly, off, amp, bp_r, bp_i in zip(prms['antpos'], prms['delays'],
            prms['offsets'], prms['amps'], prms['bp_r'], prms['bp_i']):
        antennas.append(
            a.fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, bp_r=bp_r, bp_i=bp_i)
        )
    aa = AntennaArray(prms['loc'], antennas)
    return aa

src_prms = {
        #'cyg': {
        #    'str': 10500,
        #    'index':-0.69
        #},
        #'cas': {
        #    #'str': 9150,
        #    'str': 12200,
        #    'index':-0.73,
        #    #'index': -0.93
        #},
        'Sun': {
            'str': 72900,
            'index':2.51,
            'angsize':0.00573,
        },
        #'vir': {
        #    #'str': 1446,
        #    'str': 1700,
        #    #'index': -0.86,
        #    'index': -0.75
        #},
        #'crab': {
        #    #'str':  1838,
        #    'str': 1400,
        #    #'index': -0.30,
        #    'index': -0.29
        #},
}
