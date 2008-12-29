'''This is a calibration file for data collected at PAPER in Boolardy
on JD 2454303.'''

import aipy as a, numpy as n

class AntennaArray(a.fit.AntennaArray):
    def __init__(self, *args, **kwargs):
        self.blpolys = {
            a.miriad.ij2bl(0,1): [-605821024356.96509, 929327909307.20325, -639459462462.92554, 259893778099.61926, -69088708767.947281, 12551566807.948889, -1578126290.8351059, 135589226.64294812, -7618297.2864271775, 252763.07267842517, -3760.4444961469708] ,
            a.miriad.ij2bl(1,2): [-363950777552.90387, 551299949563.00793, -374450804986.1629, 150166447205.28934, -39373210391.942276, 7052081956.5478392, -873738998.22397363, 73937559.130901143, -4089407.1546051479, 133482.11906468633, -1952.4558572234528] ,
            a.miriad.ij2bl(0,3): [-314400863696.44781, 475591742451.27417, -322626614858.6228, 129239594502.42848, -33853433151.757805, 6058476863.01019, -750134823.4077388, 63445823.632708505, -3507908.504287676, 114480.33574243639, -1674.4753298044598] ,
            a.miriad.ij2bl(0,2): [-320523940797.76831, 490784507932.59937, -337114830020.28723, 136786033418.61569, -36305448645.863884, 6585943658.0970211, -826896530.1481384, 70950328.350660875, -3981373.5992053226, 131933.53181877444, -1960.4797399176532] ,
            a.miriad.ij2bl(2,3): [-549824022273.61963, 844460153112.9895, -581769954209.30811, 236732885152.59372, -63007167261.236298, 11460311654.46904, -1442617910.4162033, 124091388.95223087, -6980358.4309154879, 231864.4113735327, -3453.4937543959531] ,
            a.miriad.ij2bl(1,3): [-437464873820.13055, 662973515098.25513, -450534294799.39746, 180780502542.15082, -47429654063.030296, 8500930209.4925852, -1054060540.6233016, 89273844.940607563, -4942442.8821271891, 161502.8750442483, -2365.2396877128194] ,
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
        #self.CsAm = n.array([        # N. Gugliucci 08/07
        #    [ 2.3541326  ,-0.0039133512 , 1.6055088e-05,-2.7468911e-08],
        #    [-0.46909345 , 0.0084471178 ,-5.1260711e-05, 1.0299793e-07],
        #    [ 0.32753617 ,-0.0081176326 , 5.8822952e-05,-1.3303273e-07],
        #    [-0.046844105, 0.00076223627,-3.5474502e-06, 4.4132035e-09],
        #    [-0.073523813, 0.0018151892 ,-1.3435102e-05, 3.1225928e-08],
        #    [ 0.047340855,-0.00085097424, 4.9799933e-06,-9.5164123e-09],
        #])
        #self.CsXo = n.array([        # N. Gugliucci 08/07
        #   [-121.29224, 1.9851554 ,-0.011876889  , 2.5222526e-05],
        #   [ 76.969303,-1.3947796 , 0.0085644354 ,-1.7448153e-05],
        #   [-36.638691, 0.93699466,-0.0068616164 , 1.5544311e-05],
        #   [ 10.189859,-0.18212180, 0.00098309486,-1.6152395e-06],
        #   [ 5.9997050,-0.15737420, 0.0012090764 ,-2.8862905e-06],
        #   [-5.6561847, 0.10468756,-0.00063126068, 1.2444705e-06],
        #])
        self.CsSd = n.array([        # N. Gugliucci 08/07
            [ 143.84525,-1.1088605  , 0.0048397670 ,-7.1054741e-06],
            [-104.00886, 1.9980993  ,-0.013304344  , 2.8955473e-05],
            [ 28.304230,-0.75088201 , 0.0056338561 ,-1.2898564e-05],
            [-8.7149717, 0.16108215 ,-0.00090283393, 1.5386691e-06],
            [-3.4672940, 0.091929321,-0.00071502397, 1.7311496e-06],
            [ 3.4123240,-0.063083812, 0.00038093617,-7.5356570e-07],
        ])
        a.ant.Beam.__init__(self, freqs, active_chans)
    def select_chans(self, active_chans):
        a.ant.Beam.select_chans(self, active_chans)
        mhz_freqs = 1e3 * self.afreqs # GHz -> MHz
        mhz_freqs = n.array([
           n.ones_like(mhz_freqs),
           mhz_freqs,
           mhz_freqs**2,
           mhz_freqs**3,
        ])
        #self.BAm = n.dot(self.CsAm, mhz_freqs)
        #self.BXo = n.dot(self.CsXo, mhz_freqs)
        self.BSd = n.dot(self.CsSd, mhz_freqs)
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
        #a1 = n.dot(A, self.BAm)
        #a2 = n.dot(A, self.BXo)
        a3 = n.dot(A, self.BSd)
        #z = (180*zang/n.pi - a2) / a3
        z = (180*zang/n.pi) / a3
        #rv = n.sqrt(a1 * n.exp(-z**2/2))
        rv = n.sqrt(n.exp(-z**2/2))
        return rv.transpose()
    def get_params(self, prm_list):
        return {}
    def set_params(self, prm_list):
        return

class Antenna(a.fit.Antenna):
    def bm_response(self, top, pol='x'):
        """Return response of beam for specified polarization."""
        top = n.array(top)
        top = {'y':top, 'x':n.dot(self.rot_pol_y, top)}[pol]
        x,y,z = top
        return self.beam.response((x,y,z))

prms = {
    'loc': ('-26:44:12.74', '116:39:59.33'), # Boolardy, Australia
    'antpos':
        #[[   0.,    0.,    0.],
        # [-105.39, 139.44,-201.57],
        # [ -68.28, 381.97,-134.25],
        # [  62.76, 462.75, 119.78]],
        #---------------------------
        [[  -2.59,  -0.07,  -0.54],
         [-103.07, 138.19,-198.50],
         [ -71.14, 381.77,-134.54],
         [  56.90, 463.28, 117.57]],
        #---------------------------
    'delays': [2.440, 1.311,  2.103, 2.245],
    'offsets': [0.003, 0.073, 0.005, 0.003],
    'amps': [1., 1., 1., 1.],
    'bp_r': n.array([
        [1.],
        [1.],
        [1.],
        [1.],
    ]),
    'bp_i': n.array([
        [0.000],
        [0.000],
        [0.000],
        [0.000],
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
    assert(len(prms['delays']) == nants and len(prms['amps']) == nants \
        and len(prms['bp_r']) == nants and len(prms['bp_i']) == nants)
    for pos, dly, off, amp, bp_r, bp_i in zip(prms['antpos'], prms['delays'],
            prms['offsets'], prms['amps'], prms['bp_r'], prms['bp_i']):
        antennas.append(
            a.fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
            #Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, bp_r=bp_r, bp_i=bp_i)
        )
    aa = AntennaArray(prms['loc'], antennas)
    return aa

src_prms = {
        'cyg': {
            'str': 8490,
            #'index':-0.09
        },
        #'cas': {
        #    #'str': 9150,
        #    'str': 12200,
        #    'index':-0.73,
        #    #'index': -0.93
        #},
        'Sun': {
            'str': 32550,
            'index':3.06,
            'angsize':0.00564,
        },
        'vir': {
            'str': 730,
            #'index': -0.86,
        },
        'crab': {
            'str': 1050,
        },
        'pic': {
            'str': 190,
        },
        'for': {
            'str': 180,
        },
        'hyd': {
            'str': 150,
        },
        'cen': {
            'str': 770,
            'angsize': 0.00450,
        },
}
