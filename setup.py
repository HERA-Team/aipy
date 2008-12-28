from distutils.core import setup, Extension
import os, glob, numpy

__version__ = open('VERSION').read().strip()

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
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    package_dir = {'aipy':'src', 'aipy.optimize':'src/optimize'},
    packages = ['aipy', 'aipy.optimize'],
    ext_modules = [
        Extension('aipy._healpix',
            ['src/_healpix/healpix_wrap.cpp', 
            'src/_healpix/cxx/Healpix_cxx/healpix_base.cc'],
            include_dirs = [numpy.get_include(), 'src/_healpix/cxx/cxxsupport',
                'src/_healpix/cxx/Healpix_cxx']),
        Extension('aipy._alm',
            ['src/_healpix/alm_wrap.cpp', 
            'src/_healpix/cxx/Healpix_cxx/alm_map_tools.cc',
            'src/_healpix/cxx/libfftpack/ls_fft.c',
            'src/_healpix/cxx/libfftpack/bluestein.c',
            'src/_healpix/cxx/libfftpack/fftpack.c',
            'src/_healpix/cxx/Healpix_cxx/healpix_map.cc',
            'src/_healpix/cxx/Healpix_cxx/healpix_base.cc'],
            include_dirs = [numpy.get_include(), 'src/_healpix/cxx/cxxsupport',
                'src/_healpix/cxx/Healpix_cxx']),
        Extension('aipy._miriad', ['src/_miriad/miriad_wrap.cpp'] + \
            indir('src/_miriad/mir', ['uvio.c','hio.c','pack.c','bug.c',
                'dio.c','headio.c','maskio.c']),
            include_dirs = [numpy.get_include(), 'src/_miriad', 
                'src/_miriad/mir']),
        Extension('aipy.utils', ['src/utils/utils.cpp'],
            include_dirs = [numpy.get_include()])
    ],
    scripts=glob.glob('scripts/*'),
    package_data = {'aipy': ['doc/*.tex', 'doc/*.png']},
)
