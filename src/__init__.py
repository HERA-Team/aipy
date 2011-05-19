"""
This package collects together tools for radio astronomical interferometry.
In addition to pure-python phasing, calibration, imaging, and
deconvolution code, this package includes interfaces to MIRIAD (a Fortran
interferometry package) and HEALPix (a package for representing spherical
data sets), and some math/fitting routines from SciPy.
All code provided is released under the GNU General Public License
(see LICENSE.txt).

Author: Aaron Parsons
"""

import phs, const, coord, deconv
import ephem, fit, healpix, img 
import interp, cal, map, miriad
import optimize, rfi, amp, scripting, src, _src, utils
import dsp
