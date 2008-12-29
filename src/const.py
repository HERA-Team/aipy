"""
Module for keeping track of physical constants.  All constants in cgs units.
See description() for a dictionary of constants and their descriptions.
"""

pi = 3.141592653589793      # Pi
c = 2.99792458e10           # Speed of Light
G = 6.6726e-8               # Newton's Gravitational Constant
h = 6.6260755e-27           # Planck constant
e = 4.8032068e-10           # Elementary Charge
m_e = 9.1093898e-28         # Electron Mass
m_p = 1.6726231e-24         # Proton Mass
k = 1.380658e-16            # Boltzmann constant
s = 5.67051e-5              # Stefan-Boltzmann constant
au = 1.496e13               # cm in 1 AU
r_sun = 6.96e10             # radius of sun (cm)
m_sun = 1.98892e33          # mass of sun (g)
ev = 1.60217733e-12         # Electron Volts
pc = 3.09e18                # cm in 1 pc
km = 1e5                    # cm in 1 km
s_per_day = 24.*60*60       # seconds in a solar day
s_per_yr = 3.1556926e7      # seconds in a year
sidereal_day = 86164.0908   # seconds
len_ns = c * 1e-9           # length of a nanosecond in cm
deg = pi / 180.             # degrees in radians
sq_deg = 3.04617e-4         # square degree in steradians
arcmin = deg / 60.          # arcminute in radians
arcsec = arcmin / 60.       # arcsecond in radians
ft = 30.48006096012         # length of a foot in cm

def description():
    return {
        'pi'          : "Pi",
        'c'           : "Speed of Light",
        'G'           : "Newton's Gravitational Constant",
        'h'           : "Planck constant",
        'e'           : "Elementary Charge",
        'm_e'         : "Electron Mass",
        'm_p'         : "Proton Mass",
        'k'           : "Boltzmann constant",
        's'           : "Stefan-Boltzmann constant",
        'au'          : "cm in 1 AU",
        'r_sun'       : "radius of sun (cm)",
        'm_sun'       : "mass of sun (g)",
        'ev'          : "Electron Volts",
        'pc'          : "cm in 1 pc",
        'km'          : "cm in 1 km",
        's_per_day'   : "seconds in a solar day",
        's_per_yr'    : "seconds in a year",
        'sidereal_day': "seconds in a sidereal day",
        'len_ns'      : "length of a nanosecond in cm",
        'deg'         : "degrees in radians",
        'sq_deg'      : "square degree in steradians",
        'arcmin'      : "arcminute in radians",
        'arcsec'      : "arcsecond in radians",
        'ft'          : "length of a foot in cm",
    }
