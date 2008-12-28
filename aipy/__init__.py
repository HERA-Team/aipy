"""
This package collects together tools for radio astronomical interferometry.  In
addition to pure-python phasing, calibration, imaging, and
deconvolution code, this package includes interfaces to MIRIAD (a Fortran
interferometry package), HEALPix (a package for representing spherical data
sets), routines from SciPy for fitting, and the PyFITS and PyEphem packages
verbatim.  All code provided is released under the GNU General Public License
(see LICENSE.txt).

Author: Aaron Parsons
Date: 2/12/07
Revisions: lots
"""

__version__ = '0.4.1'
import ant, const, coord, deconv, eor, fit, img, loc, sim, src
import miriad, interp, optimize, pyephem, pyfits, utils, healpix
