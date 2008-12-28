from distutils.core import setup, Extension
import numpy, os, glob

__version__ = '0.4.1'

def indir(path, files):
    return [os.path.join(path, f) for f in files]

setup(name = 'aipy',
    version = __version__,
    description = 'Astronomical Interferometry in PYthon',
    long_description = \
"""
This package collects together tools for radio astronomical interferometry.  In
addition to pure-python phasing, calibration, imaging, and
deconvolution code, this package includes interfaces to MIRIAD (a Fortran
interferometry package), HEALPix (a package for representing spherical data
sets), and fitting routines from SciPy.
""",
    license = 'GPL',
    author = 'Aaron Parsons',
    author_email = 'aparsons@astron.berkeley.edu',
    url = 'http://setiathome.berkeley.edu/~aparsons/aipy',
    classifiers = [
        'Development Status :: 2 - Immature',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    packages = ['aipy', 'aipy.optimize'],
    ext_modules = [
        Extension('aipy._healpix',
            ['aipy/_healpix/healpix_wrap.cpp',
            'aipy/_healpix/cxx/Healpix_cxx/healpix_base.cc'],
            include_dirs = [numpy.get_include(), 'aipy/_healpix/cxx/cxxsupport',
                'aipy/_healpix/cxx/Healpix_cxx']),
        Extension('aipy._miriad', ['aipy/_miriad/miriad_wrap.cpp'] + \
            indir('aipy/_miriad/mir', ['uvio.c','hio.c','pack.c','bug.c',
                'dio.c','headio.c','maskio.c']),
            include_dirs = [numpy.get_include(),
                'aipy/_miriad', 'aipy/_miriad/mir']),
        Extension('aipy.utils', ['aipy/utils/utils.cpp'],
            include_dirs = [numpy.get_include()])
    ],
    scripts=glob.glob('scripts/*'),
    package_data = {'aipy': ['doc/*.tex', 'doc/*.png']},
)
