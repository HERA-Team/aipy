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

from __future__ import print_function, division, absolute_import

from . import phs, const, coord, deconv
from . import fit, healpix, img
from . import interp, cal, map, miriad
import scipy.optimize as optimize
from . import rfi, amp, scripting, src, _src, utils
from . import dsp
from . import pol, twodgauss #added by dfm
import ephem

try:
    from .__gitlog__ import __gitlog__
    from .__branch__ import __branch__
    from .__version__ import __version__
except ImportError:
    __gitlog__ = None
    __branch__ = None
    fh = open('VERSION', 'r')
    __version__ = fh.read()
    fh.close()
