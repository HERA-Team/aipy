"""Latest location-specific information."""
import numpy, aipy.fit, aipy.eor

locations = {
    'pwa303': ('-26:44:12.74', '116:39:59.33'), # Boolardy, Australia
    'pgb371': ( '38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
}

antpos = {
    'pwa303':
        [[    0.,    0.,    0.],
         [-100.7, 139.1,-197.7],
         [ -68.7, 383.4,-132.6],
         [  59.4, 465.3, 120.7],],
    'pgb371':
        [[  -8.5, 455.3,   9.8],
         [ 205.5, 319.5,-251.7],
         [ 187.1,-353.0,-232.6],
         [-262.7,-219.1, 318.7],
         [-293.4,  -7.7, 360.2],
         [-286.0,  93.2, 352.2],
         [-182.7, 353.2, 227.6],
         [ -84.2, 434.5, 107.2],],
}

delays = {
    'pwa303': [0., -1.62, -1.41, -1.67],
    'pgb371': [0., 0., 0., 0., 0., 0., 0., 0.],
}
        
offsets = {
    'pwa303': [0., 0.16, -0.20, 0.28],
    'pgb371': [0., 0., 0., 0., 0., 0., 0., 0.],
}

amplitudes = {
    'pwa303': [3.7e-7, 3.2e-7, 3.4e-7, 3.3e-7],
    'pgb371': [1., 1., 1., 1., 1., 1., 1., 1.],
}

pwa303_poly = [  2.14249664e+14,  -3.27824728e+14,   2.25052468e+14,
                -9.12789456e+13,   2.42212141e+13,  -4.39353531e+12,
                 5.51696241e+11,  -4.73521292e+10,   2.65851097e+09,
                -8.81582278e+07,   1.31113443e+06]
passbands = {
    'pwa303':
        [pwa303_poly,
         pwa303_poly,
         pwa303_poly,
         pwa303_poly,],
    'pgb371':
        [[1.],
         [1.],
         [1.],
         [1.],
         [1.],
         [1.],
         [1.],
         [1.],]
}

def get_freqs(sdf, sfreq, nchan):
    return numpy.arange(nchan, dtype=numpy.float) * sdf + sfreq

def get_aa(loc_key, sdf, sfreq, nchan, use_bp=True, use_ants=None):
    """Return an antenna array with the latest-greatest fit parameters for
    a given location.
        loc_key: The name of the antenna array.  Currently: 'pwa303','pgb371'
        sdf, sfreq, nchan: Their usual Miriad-defined meanings
        use_bp: Use bandpass information, or set to 1 (which is what you want
            if you've applied a bp calibration already).
        use_ants: A list of ant indices you want in the array if 
            you don't want all of them."""
    location = locations[loc_key]
    freqs = get_freqs(sdf, sfreq, nchan)
    beam = aipy.eor.Beam(freqs)
    antennas = []
    for pos, dly, off, amp, gpoly in zip(antpos[loc_key], delays[loc_key], 
            offsets[loc_key], amplitudes[loc_key], passbands[loc_key]):
        if not use_bp:
            amp = 1
            gpoly = [1]
        antennas.append(
            aipy.fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, gain_poly=gpoly)
        )
    if not use_ants is None: antennas = [antennas[i] for i in use_ants]
    return aipy.fit.AntennaArray(antennas, location)
