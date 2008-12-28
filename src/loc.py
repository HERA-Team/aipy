"""Latest location-specific information."""

import numpy as n, ant, sim, fit, os, coord

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
    def select_chans(self, active_chans):
        if active_chans is None: active_chans = n.arange(self.freqs.size)
        ant.Beam.select_chans(self, active_chans)
        self.BAm_sel = self.BAm.take(active_chans, axis=1)
        self.BXo_sel = self.BXo.take(active_chans, axis=1)
        self.BSd_sel = self.BSd.take(active_chans, axis=1)
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

#  _                    _   _                 
# | |    ___   ___ __ _| |_(_) ___  _ __  ___ 
# | |   / _ \ / __/ _` | __| |/ _ \| '_ \/ __|
# | |__| (_) | (_| (_| | |_| | (_) | | | \__ \
# |_____\___/ \___\__,_|\__|_|\___/|_| |_|___/

locations = {
    'pwa303': {
        'loc': ('-26:44:12.74', '116:39:59.33'), # Boolardy, Australia
        'antpos': 
            #[[    0.,    0.,    0.],
            # [-100.7, 139.1,-197.7],
            # [ -68.7, 383.4,-132.6],
            # [  59.4, 465.3, 120.7],],
            [[    0.,    0.,    0.],
             [-100.67, 138.14,-198.28],
             [ -68.62, 381.81,-133.88],
             [  59.4, 465.3, 120.0],],
        #'delays': [0., -1.62, -1.41, -1.67],
        'delays': [0., -0.355, -0.400, -1.72],
        #'offsets': [0., 0.16,-0.20, 0.28],
        #'offsets':[0., .65, 0., -.31],
        'offsets':[0., 0, 0., 0],
        #'amps': [3.7e-7, 3.2e-7, 3.4e-7, 3.3e-7],
        'amps': [1.27e-2,1.27e-2,1.27e-2,1.27e-2],
        'passbands': [None, None, None, None,],
        'beam': BeamNoFlaps,
    },
    'pgb371': {
        'loc':
            ( '38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
        'antpos':
            #-Full Bandwidth----------------
            [[   0.00,    0.00,    0.00],
             [ 214.27, -136.52, -260.96],
             [ 196.53, -810.27, -253.55],
             [-255.16, -674.34,  308.17],
             [-285.31, -463.01,  350.48],
             [-277.44, -361.92,  342.78],
             [-174.60, -102.56,  218.51],
             [ -76.55,  -20.95,   98.04]],
            #-Channels 16_20----------------
            #[[   0.00,    0.00,    0.00],
            # [ 214.33, -136.53, -260.86],
            # [ 196.59, -810.22, -253.46],
            # [-255.06, -674.31,  308.20],
            # [-285.19, -462.93,  350.53],
            # [-277.39, -361.83,  342.84],
            # [-174.49, -102.49,  218.60],
            # [ -76.55,  -20.95,   98.04]],
            #-From Miriad-------------------
            #[[  -8.5, 455.3,   9.8],
            # [ 205.5, 319.5,-251.7],
            # [ 187.1,-353.0,-232.6],
            # [-262.7,-219.1, 318.7],
            # [-293.4,  -7.7, 360.2],
            # [-286.0,  93.2, 352.2],
            # [-182.7, 353.2, 227.6],
            # [ -84.2, 434.5, 107.2],],
        'delays':
            # Full Bandwidth
            [0.000,-7.318, 6.996, 8.309, 9.467, 2.992, 3.311, 0.000],
            # Channels 16_20:
            #[0.000,-7.287, 6.994, 8.307, 9.472, 2.913, 3.310, 0.000],
        'offsets':
            [0., 0., 0., 0., 0., 0., 0., 0.],
        'amps':
            # BeamPolynomial amplitudes:
            [ 0.00403,0.00395,0.00446,0.00484,0.00415,0.00363,0.00397,0],
            # BeamAlm amplitudes:
            #[1.84e-3,1.80e-3,2.04e-3,2.21e-3,1.89e-3,1.66e-3,1.82e-3,0],
        'passbands': n.array([
[1.802171753471304, 3.0186259708950001, 1.0960425238684581, -3.9870318256746744, 6.8264344645435893, 1.7216886128894933, -7.9148956235162498, 12.465724304735865, -5.9591673812243133, 0.99880841937299358, 2.1640733575225459, 5.495135638552556, -10.936437033313222, 13.141365226309759, -10.48412089384102, 7.2682339557625157]
            ] * 8),
        'beam':
            fit.BeamAlmSymm,
            #fit.BeamPolynomial,
        'bm_poly':
            n.reshape(n.array(
[1.2046020370915436, -2.6095786668890595, 2.5925890523798021, -0.080803171031906151, 0.57181294299205443, 1.4389047479023493, -0.13582588411664062, 0.4091040646302877, 1.6046483772054452, -0.15521824743919593, 0.71629841628670787, 0.66484890913174732, -0.13185943216178289, 0.73602913691468719, 0.38414058012426966, -0.053400811201470243, 0.2543327969236715, -0.07065113699453246]
#[0.89418135994583947, 0.060316723460104658, -0.025467121762497075, -0.0079476086198618338, 0.0045975567943226502, -0.00400680525014051720]
), (6,3)),
        'bm_alm_r':
            n.array(
#symm[-157.42163055730106, 179.68485987960415, -188.12965192074896, -91.181323654366906, -3538.6721276536782, -2.4060147092585322, 21.755281677401666, 2467.3794791284031, 0.57926173051600083, 1898.2846250007442, 0.52686763619333321, 0.098384478053203883, 0.21499805987779197, 0.2711207710159248, 0.12232923060914505, -0.41079125965142715, 0.20515817284434473, -0.15735378929100907, -0.55462568835925175, -0.54687363391154853, 0.11352024768740611]
[-156.65755374951397, 179.68831584914614, 5.5653376277512123, -91.346400400027392, 42.56805290955279, 2.7092897401498597, 21.684652478272604, 68.293123959114368, -1.5920911942298033, 25.619750025903279]
#symm[-182.10558546757943, 4750873.5085385637, 207.53712465200925, -15359.801959513701, -216427.35055944428, 4146.8211150462866, -104.69663733602167, 938334.66346961982, 40544.709072240745, 94455.682240035909, 1.8269692380704008, 0.22549283308792439, 24.367608890891319, -332408.08259872429, -49107.753493655298, -1040171.9054042659, -1.1274397197045212, -0.1095578784275146, 38018.21079111805, 37067.877460584714]
            ),
    }
}

def get_freqs(sdf, sfreq, nchan):
    return n.arange(nchan, dtype=n.float) * sdf + sfreq

def lookup_key(loc_key):
    global locations
    if os.path.exists(loc_key): locations.update(eval(open(loc_key).read()))
    return locations[loc_key]

def get_aa(loc_key, sdf, sfreq, nchan, use_bp=True, use_ants=None):
    """Return an antenna array with the latest-greatest fit parameters for
    a given location.
        loc_key: The name of the antenna array.  Currently: 'pwa303','pgb371'
        sdf, sfreq, nchan: Their usual Miriad-defined meanings
        use_bp: Use bandpass information, or set to 1 (which is what you want
            if you've applied a bp calibration already).
        use_ants: A list of ant indices you want in the array if 
            you don't want all of them."""
    dat = lookup_key(loc_key)
    freqs = get_freqs(sdf, sfreq, nchan)
    beam = dat['beam'](freqs)
    try: beam.set_params(dat)
    except(AttributeError): pass
    location = dat['loc']
    antennas = []
    nants = len(dat['antpos'])
    assert(len(dat['delays']) == nants and len(dat['offsets']) == nants \
        and len(dat['amps']) == nants and len(dat['passbands']) == nants)
    for pos, dly, off, amp, dec_bp in zip(dat['antpos'], dat['delays'], 
            dat['offsets'], dat['amps'], dat['passbands']):
        if not use_bp:
            amp = 1
            dec_bp = None
        antennas.append(
            fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, bp=dec_bp)
        )
    if not use_ants is None: antennas = [antennas[i] for i in use_ants]
    return fit.AntennaArray(antennas, dat['loc'])

def get_loc(loc_key):
    """Return an array location (with information about lat and long)."""
    return ant.ArrayLocation(lookup_key(loc_key)['loc'])
