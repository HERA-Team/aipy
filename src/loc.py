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
    'pgb220': {
        'loc': ( '38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
        'antpos': 
            [[  -8.5, 455.3,   9.8],
             [ 205.3, 319.6,-250.8],
             [-262.1,-218.8, 319.6],
             [-292.3,  -7.5, 360.6],],
        'delays': [0.,  3.86,-14.51,  2.87],
        'offsets': [0., 0.07, 0.8, 0.0085],
        'amps': [3e-5, 3e-5, 3e-5, 3e-5],
        'passbands': [None, None, None, None,],
        'beam': BeamNoFlaps,
    },
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
            [[   0.00,    0.00,    0.00],
             [ 214.35, -136.52, -260.94],
             [ 196.65, -810.26, -253.52],
             [-254.93, -674.28,  308.22],
             [-285.12, -462.90,  350.54],
             [-277.39, -361.81,  342.78],
             [-174.46, -102.48,  218.57],
             [ -76.55,  -20.95,   98.04]],
            #[[  -7.65, 455.45,   9.16],
            # [ 207.13, 318.92,-251.30],
            # [ 188.60,-354.63,-244.97],
            # [-261.59,-218.96, 318.93],
            # [-289.35,  -8.19, 364.48],
            # [-282.25,  93.35, 355.49],
            # [-182.52, 352.99, 227.22],
            # [ -84.2,  434.5,  107.2],],
            #[[  -8.5, 455.3,   9.8],
            # [ 205.5, 319.5,-251.7],
            # [ 187.1,-353.0,-232.6],
            # [-262.7,-219.1, 318.7],
            # [-293.4,  -7.7, 360.2],
            # [-286.0,  93.2, 352.2],
            # [-182.7, 353.2, 227.6],
            # [ -84.2, 434.5, 107.2],],
        'delays':
            [0.000,-7.249, 7.001, 8.222, 9.431, 2.959, 3.316, 0.000],
            #[0.000,-7.959, 7.913, 6.512, 3.755,-1.198, 4.103, 0.000],
            #[-0.077,-8.065, 7.721, 6.315, 3.348,-1.317, 3.656, 0.],
        'offsets':
            [0., 0., 0., 0., 0., 0., 0., 0.],
        'amps':
            [1.87e-3,1.73e-3,1.94e-3,2.24e-3,1.93e-3,1.74e-3,1.79e-3],
        'passbands': n.array([
[0.92119178272251157, 0.70497609872015965, 1.213526285569519, -1.4818472395827489, 5.7261552261104267, -2.949995608568341, 1.065070150093482, 2.9842814038554168, 1.020055815261081, -1.1386388014485318, 3.0543822393482918, 1.4548931649749526, -1.3614683369296627, 3.2394754402047674, -0.2104689321640455, 1.558789743003282, 0.33640002574916428, 1.7019925021202851, 0.65038791604530655, 0.68362881815210019, 2.2990169585786182, -1.2988081697133809, 3.7291631117304949, -1.9701750509709797, 3.3457765335522653, -0.51929454836969557, 0.17408396844198432, 0.15209545054936813, -0.43296126160353721, 0.79276982046468736, 1.0812217729759541, 1.9379225585315418]
            ] * 8),
        'beam':
            fit.BeamAlm,
        'bm_alm':
            n.array(
[-181.53833842631752, 1706811.9585676491, 207.56333528582263, 7224.1791133664701, -161528.35278112761, 23489.296631910096, -105.2948013170477, -816080.89344837982, 46077.100205250637, 83787.666568546265, 1.9075559723296482, 0.060466706544603094, 24.663307473598323, -93658.474087108683, -115926.6440020618, -399954.1085421748, -1.0104477695504377, -0.060896525507602911, 14700.834366869836, 23969.279351218662]
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
