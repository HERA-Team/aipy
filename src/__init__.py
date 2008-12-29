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

import ant, const, coord, deconv, fit, healpix, img, interp, loc, map, miriad
import optimize, rfi, sim, scripting, src, utils
