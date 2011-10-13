from distutils.core import setup, Extension

import os, glob, numpy, sys
if 'upload' in sys.argv or 'register' in sys.argv:
    from ez_setup import use_setuptools; use_setuptools()
    from setuptools import setup, Extension

__version__ = open('VERSION').read().strip()

def get_description():
    lines = [L.strip() for L in open('README').readlines()]
    d_start = None
    for cnt, L in enumerate(lines):
        if L.startswith('DESCRIPTION'): d_start = cnt + 1
        elif not d_start is None:
            if len(L) == 0: return ' '.join(lines[d_start:cnt])
    raise RuntimeError('Bad README')

def indir(path, files):
    return [os.path.join(path, f) for f in files]

setup(name = 'aipy',
    version = __version__,
    description = 'Astronomical Interferometry in PYthon',
    long_description = get_description(),
    license = 'GPL',
    author = 'Aaron Parsons',
    author_email = 'aparsons@astron.berkeley.edu',
    url = 'http://setiathome.berkeley.edu/~aparsons/aipy/aipy.cgi',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    setup_requires = ['numpy>=1.2'],
    install_requires = ['pyephem>=3.7.3.2', 'pyfits>=2.1', 'numpy>=1.2'],
    dependency_links = [
        'http://www.stsci.edu/resources/software_hardware/pyfits'
    ],
    package_dir = {'aipy':'src', 'aipy.optimize':'src/optimize', 'aipy._src':'src/_src'},
    packages = ['aipy', 'aipy.optimize','aipy._src'],
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
        Extension('aipy._deconv', ['src/_deconv/deconv.cpp'],
            include_dirs = [numpy.get_include()]),
        #Extension('aipy._img', ['src/_img/img.cpp'],
        #    include_dirs = [numpy.get_include()]),
        Extension('aipy._dsp', ['src/_dsp/dsp.c', 'src/_dsp/grid/grid.c'],
            include_dirs = [numpy.get_include(), 'src/_dsp', 'src/_dsp/grid']),
        Extension('aipy.utils', ['src/utils/utils.cpp'],
            include_dirs = [numpy.get_include()]),
        Extension('aipy._cephes',
            ['src/_cephes/_cephesmodule.c', 'src/_cephes/ufunc_extras.c'] + \
            glob.glob('src/_cephes/cephes/*.c') + \
            glob.glob('src/_cephes/c_misc/*.c'),
            include_dirs = [numpy.get_include()]),
    ],
    scripts=glob.glob('scripts/*'),
    package_data = {'aipy': ['_src/*.txt']},
)
