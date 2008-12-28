"""
A package for viewing, calibrating, and imaging interferometric data.
Includes interfacing to MIRIAD data files and the LBFGSB constrained
function minimization algorithm (though this minimizer doesn't work so well
and may be scrapped from the package someday).  This package also contains
some general-purpose deconvolution algorithms.

All code provided (with the exception of LBFGSB, which has its own usage
policy--see documentation) is released under the GNU General Public License
(see LICENSE.txt).

Author: Aaron Parsons
Date: 2/12/07
Revisions: lots
"""

__version__ = '0.3.0'
import ant, constants, deconv, eor, fit, img, loc, sim, src
import miriad, interpolate, optimize, pyephem
