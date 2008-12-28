'''This is a calibration file for data collected at PAPER in Green Bank
on JD 2454371.'''

import aipy as a, numpy as n

prms = {
    'loc': ( '38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
    'antpos':
        [[   0.00,    0.00,    0.00],
         [ 214.27, -136.52, -260.96],
         [ 196.53, -810.27, -253.55],
         [-255.16, -674.34,  308.17],
         [-285.31, -463.01,  350.48],
         [-277.44, -361.92,  342.78],
         [-174.60, -102.56,  218.51],
         [ -76.55,  -20.95,   98.04]],
    'delays': [0.000,-7.318, 6.996, 8.309, 9.467, 2.992, 3.311, 0.000],
    'offsets': [0., 0., 0., 0., 0., 0., 0., 0.],
    'amps': [ 0.00403,0.00395,0.00446,0.00484,0.00415,0.00363,0.00397,0],
    'passbands': n.array([
        [1,] * 8
    ]),
    'beam': a.fit.BeamPolynomial,
    'bm_poly': n.reshape(n.array(
        [1.2046020370915436, -2.6095786668890595, 2.5925890523798021, -0.080803171031906151, 0.57181294299205443, 1.4389047479023493, -0.13582588411664062, 0.4091040646302877, 1.6046483772054452, -0.15521824743919593, 0.71629841628670787, 0.66484890913174732, -0.13185943216178289, 0.73602913691468719, 0.38414058012426966, -0.053400811201470243, 0.2543327969236715, -0.07065113699453246]
        ), (6,3)),
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

