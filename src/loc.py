"""Latest location-specific information."""
import numpy, ant, eor, fit, os

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
    },
    'pwa303': {
        'loc': ('-26:44:12.74', '116:39:59.33'), # Boolardy, Australia
        'antpos': 
            #[[    0.,    0.,    0.],
            # [-100.7, 139.1,-197.7],
            # [ -68.7, 383.4,-132.6],
            # [  59.4, 465.3, 120.7],],
            [[    0.,    0.,    0.],
             [-100.7, 139.1,-197.8],
             [ -68.7, 383.4,-133.0],
             [  59.4, 465.3, 120.0],],
        #'delays': [0., -1.62, -1.41, -1.67],
        'delays': [0., -1.56, -1.46, -1.72],
        #'offsets': [0., 0.16,-0.20, 0.28],
        'offsets':[0., .65, 0., -.31],
        #'amps': [3.7e-7, 3.2e-7, 3.4e-7, 3.3e-7],
        'amps': [10.0e-5,10.0e-5, 9.3e-5, 9.9e-5],
        'passbands': [None, None, None, None,],
    },
    'pgb371': {
        'loc': ( '38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
        'antpos':
            [[  -8.5, 455.3,   9.8],
             [ 205.8, 318.6,-252.0],
             [ 189.0,-355.2,-241.6],
             [-262.7,-219.2, 319.2],
             [-292.2,  -7.9, 366.0],
             [-286.0,  93.2, 356.8],
             [-182.2, 352.7, 232.8],
             [ -84.2, 434.5, 107.2],],
            #[[  -8.5, 455.3,   9.8],
            # [ 205.5, 319.5,-251.7],
            # [ 187.1,-353.0,-232.6],
            # [-262.7,-219.1, 318.7],
            # [-293.4,  -7.7, 360.2],
            # [-286.0,  93.2, 352.2],
            # [-182.7, 353.2, 227.6],
            # [ -84.2, 434.5, 107.2],],
        'delays': [0., -0.02, 4.78, 0.05, -2.42, 0.01, -0.51, 0.],
        #'delays': [0., 0., -1.16, 1.38, -1.42, 0., -0.14, 0.],
        'offsets': [0., 0., -0.04, 0.02, 0.01, 0., 0., 0.],
        'amps': [1., 1., 1., 1., 1., 1., 1., 1.],
        'passbands': [None, None, None, None, None, None, None, None,]
    }
}

def get_freqs(sdf, sfreq, nchan):
    return numpy.arange(nchan, dtype=numpy.float) * sdf + sfreq

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
    beam = eor.Beam(freqs)  # This may need to be changed...
    location = dat['loc']
    antennas = []
    for pos, dly, off, amp, dec_bp in zip(dat['antpos'], dat['delays'], 
            dat['offsets'], dat['amps'], dat['passbands']):
        if not use_bp:
            amp = 1
            dec_bp = None
        antennas.append(
            fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, dec_bp=dec_bp)
        )
    if not use_ants is None: antennas = [antennas[i] for i in use_ants]
    return fit.AntennaArray(antennas, dat['loc'])

def get_loc(loc_key):
    """Return an array location (with information about lat and long)."""
    return ant.ArrayLocation(lookup_key(loc_key)['loc'])
