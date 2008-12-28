"""Latest location-specific information."""
import numpy, ant, eor, fit

locations = {
    'pgb220': ( '38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
    'pwa303': ('-26:44:12.74', '116:39:59.33'), # Boolardy, Australia
    'pgb371': ( '38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
}

antpos = {
    'pgb220':
        [[  -8.5, 455.3,   9.8],
         [ 205.3, 319.6,-250.8],
         [-262.1,-218.8, 319.6],
         [-292.3,  -7.5, 360.6],],
    'pwa303':
        #[[    0.,    0.,    0.],
        # [-100.7, 139.1,-197.7],
        # [ -68.7, 383.4,-132.6],
        # [  59.4, 465.3, 120.7],],
        [[    0.,    0.,    0.],
         [-100.7, 139.1,-197.8],
         [ -68.7, 383.4,-133.0],
         [  59.4, 465.3, 120.0],],
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
    'pgb220': [0.,  3.86,-14.51,  2.87],
    #'pwa303': [0., -1.62, -1.41, -1.67],
    'pwa303': [0., -1.56, -1.46, -1.72],
    'pgb371': [0., 0., 0., 0., 0., 0., 0., 0.],
}
        
offsets = {
    'pgb220': [0., 0.07, 0.8, 0.0085],
    #'pwa303': [0., 0.16,-0.20, 0.28],
    'pwa303':[0., .65, 0., -.31],
    'pgb371': [0., 0., 0., 0., 0., 0., 0., 0.],
}

amplitudes = {
    'pgb220': [3e-5, 3e-5, 3e-5, 3e-5],
    #'pwa303': [3.7e-7, 3.2e-7, 3.4e-7, 3.3e-7],
    'pwa303': [10.0e-5,10.0e-5, 9.3e-5, 9.9e-5],
    'pgb371': [1., 1., 1., 1., 1., 1., 1., 1.],
}

#pwa303_poly = [  2.14249664e+14,  -3.27824728e+14,   2.25052468e+14,
#                -9.12789456e+13,   2.42212141e+13,  -4.39353531e+12,
#                 5.51696241e+11,  -4.73521292e+10,   2.65851097e+09,
#                -8.81582278e+07,   1.31113443e+06]
pwa303_spline = (
     numpy.array([ 0.1125,  0.1125    ,  0.1125    ,  0.1125    ,  0.1171875 ,
        0.121875  ,  0.13125   ,  0.1359375 ,  0.140625  ,  0.15      ,
        0.1546875 ,  0.159375  ,  0.1640625 ,  0.16875   ,  0.1734375 ,
        0.178125  ,  0.18742676,  0.18742676,  0.18742676,  0.18742676]),
    numpy.array([ 0.20047258, 0.13553328,  0.47578084,  1.54538574,  1.36744419,
        1.28911202,  1.18573184,  1.13628592,  1.01282591,  1.04695018,
        0.95181596,  1.06482914,  1.08785934,  0.86663611,  0.99864105,
        0.69673599]),
    3
)

passbands = {
    'pgb220': [None, None, None, None,],
    'pwa303': [pwa303_spline, pwa303_spline, pwa303_spline, pwa303_spline,],
    'pgb371': [None, None, None, None, None, None, None, None,]
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
    beam = eor.Beam(freqs)
    antennas = []
    for pos, dly, off, amp, spline in zip(antpos[loc_key], delays[loc_key], 
            offsets[loc_key], amplitudes[loc_key], passbands[loc_key]):
        if not use_bp:
            amp = 1
            gspline = [1., 1.]
        antennas.append(
            fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, spline=spline)
        )
    if not use_ants is None: antennas = [antennas[i] for i in use_ants]
    return fit.AntennaArray(antennas, location)

def get_loc(loc_key):
    """Return an array location (with information about lat and long).
    loc_key can be 'pwa303' or 'pgb371'."""
    location = locations[loc_key]
    return ant.ArrayLocation(location)
