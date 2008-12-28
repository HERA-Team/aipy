"""Latest location-specific information."""
import numpy, aipy.fit, aipy.eor

def get_freqs(nchan=64):
    sfreq = 0.1126
    sdf = 0.00234 / 2
    freqs = numpy.arange(nchan, dtype=numpy.float) * sdf + sfreq
    return freqs

# Boolardy JD2454303
bp_data = numpy.array([
#    CH     ANT0         ANT1       ANT2        ANT3        FIT
    [ 0,   2.69e-9,    2.48e-9,    1.82e-9,    2.55e-9,    1.82e-14],
    [ 5,   2.10e-8,    2.06e-8,    2.32e-8,    2.26e-8,    1.41e-13],
    [10,   4.21e-8,    3.57e-8,    3.86e-8,    3.77e-8,    9.74e-13],
    [15,   3.98e-8,    3.68e-8,    3.67e-8,    3.60e-8,    6.38e-13],
    [20,   3.80e-8,    3.27e-8,    3.21e-8,    3.17e-8,    4.00e-13],
    [30,   3.23e-8,    2.69e-8,    2.73e-8,    2.63e-8,    1.34e-13],
    [35,   3.00e-8,    2.40e-8,    2.56e-8,    2.55e-8,    9.90e-14],
    [40,   2.89e-8,    2.52e-8,    2.48e-8,    2.54e-8,    9.33e-14],
    [45,   2.88e-8,    2.47e-8,    2.61e-8,    2.55e-8,    8.23e-14],
    [50,   2.76e-8,    2.33e-8,    2.79e-8,    2.50e-8,    6.94e-14],
    [55,   2.72e-8,    2.40e-8,    2.77e-8,    2.41e-8,    6.10e-14],
    [60,   2.62e-8,    2.40e-8,    2.43e-8,    2.31e-8,    4.66e-14],
    [63,   2.27e-8,    2.09e-8,    2.14e-8,    2.29e-8,    4.28e-14],
])

freqs = get_freqs()
amps = bp_data[:,1:5].sum(axis=0)  # Latest estimate of antenna amplitudes
mfreqs = freqs.take(bp_data[:,0].astype(numpy.int))
avg_bp = bp_data[:,1:5] / amps 
avg_bp = numpy.average(avg_bp, axis=1)
bp_poly = numpy.polyfit(mfreqs, avg_bp, 10) # Latest estimate of passbands

location = (-0.466646451963, 2.03621421241) # Boolardy
b = aipy.eor.Beam(freqs)

def get_aa(use_bp=True):
    global amps
    global bp_poly
    if not use_bp:
        amps = numpy.ones_like(amps)
        bp_poly = [1]
    ants = [
        aipy.fit.Antenna(    0.,    0.,    0.,b,delay=+0.0, offset=+0.0,
            amp=amps[0],
            gain_poly=bp_poly),
        aipy.fit.Antenna(-100.7, 139.1,-197.7,b,delay=-1.62,offset=+0.16,
            amp=amps[1],
            gain_poly=bp_poly),
        aipy.fit.Antenna( -68.7, 383.4,-132.6,b,delay=-1.41,offset=-0.20,
            amp=amps[2],
            gain_poly=bp_poly),
        aipy.fit.Antenna(  59.4, 465.3, 120.7,b,delay=-1.67,offset=+0.28,
            amp=amps[3],
            gain_poly=bp_poly),
    ]
    return aipy.fit.AntennaArray(ants, location)
