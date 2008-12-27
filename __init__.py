"""
A package for viewing and calibrating interferometric data.

Author: Aaron Parsons
Email: aparsons at astron.berkeley.edu
Date: 2/12/07
"""

__all__ = ['antennas', 'constants', 'params', 'rfi', 'fit']
__version__ = '0.1.0'
from numpy._import_tools import PackageLoader
pl = PackageLoader()
pl(verbose=False, postpone=True)
__doc__ += pl.get_pkgdocs()
