#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

from setuptools import setup, Extension

import os, glob, numpy, subprocess

print("Generating aipy_src/__version__.py: ", end='')
__version__ = open('VERSION').read().strip()
print(__version__)
open('aipy_src/__version__.py','w').write('__version__="%s"'%__version__)

#read the latest git status out to an installed file
try:
#    gitbranch = subprocess.check_output('git symbolic-ref -q HEAD',shell=True, cwd='.').strip().split('/')[-1]
    gitbranch = os.popen('git symbolic-ref -q HEAD').read().strip()
    print("Generating aipy_src/__branch__.py")
#    gitlog = subprocess.check_output('git log -n1 --pretty="%h%n%s%n--%n%an%n%ae%n%ai"',shell=True, cwd='.').strip()
    gitlog = os.popen('git log -n1 --pretty="%h%n%s%n--%n%an%n%ae%n%ai"').read().strip()
    print("Generating aipy_src/__gitlog__.py.")
    print(gitlog)
except:
    gitbranch = "unknown branch"
    gitlog = "git log not found"
open('aipy_src/__branch__.py','w').write('__branch__ = \"%s\"'%gitbranch)
open('aipy_src/__gitlog__.py','w').write('__gitlog__ = \"\"\"%s\"\"\"'%gitlog)


def get_description():
    lines = [L.strip() for L in open('README.md').readlines()]
    d_start = None
    for cnt, L in enumerate(lines):
        if L.startswith('## Description'): d_start = cnt + 1
        elif not d_start is None:
            if len(L) == 0: return ' '.join(lines[d_start:cnt])
    raise RuntimeError('Bad README')

def indir(path, files):
    return [os.path.join(path, f) for f in files]

global_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

setup(name = 'aipy',
    version = __version__,
    description = 'Astronomical Interferometry in PYthon',
    long_description = get_description(),
    license = 'GPL',
    author = 'Aaron Parsons',
    author_email = 'aparsons@berkeley.edu',
    url = 'http://github.com/HERA-Team/aipy',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    setup_requires = ['numpy>=1.2'],
    install_requires = ['pyephem>=3.7.3.2', 'astropy>=1.0 <2.0', 'numpy>=1.2', 'scipy>=0.19'],
    package_dir = {'aipy':'aipy_src', 'aipy._src':'aipy_src/_src'},
    packages = ['aipy', 'aipy._src'],
    ext_modules = [
        Extension('aipy._healpix',
            ['aipy_src/_healpix/healpix_wrap.cpp',
            'aipy_src/_healpix/cxx/Healpix_cxx/healpix_base.cc'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy_src/_healpix/cxx/cxxsupport',
                            'aipy_src/_healpix/cxx/Healpix_cxx', 'aipy_src/_common']),
        Extension('aipy._alm',
            ['aipy_src/_healpix/alm_wrap.cpp',
            'aipy_src/_healpix/cxx/Healpix_cxx/alm_map_tools.cc',
            'aipy_src/_healpix/cxx/libfftpack/ls_fft.c',
            'aipy_src/_healpix/cxx/libfftpack/bluestein.c',
            'aipy_src/_healpix/cxx/libfftpack/fftpack.c',
            'aipy_src/_healpix/cxx/Healpix_cxx/healpix_map.cc',
            'aipy_src/_healpix/cxx/Healpix_cxx/healpix_base.cc'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy_src/_healpix/cxx/cxxsupport',
                'aipy_src/_healpix/cxx/Healpix_cxx', 'aipy_src/_common']),
        Extension('aipy._miriad', ['aipy_src/_miriad/miriad_wrap.cpp'] + \
            indir('aipy_src/_miriad/mir', ['uvio.c','hio.c','pack.c','bug.c',
                'dio.c','headio.c','maskio.c']),
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy_src/_miriad',
                'aipy_src/_miriad/mir', 'aipy_src/_common']),
        Extension('aipy._deconv', ['aipy_src/_deconv/deconv.cpp'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy_src/_common']),
        #Extension('aipy._img', ['aipy_src/_img/img.cpp'],
        #    include_dirs = [numpy.get_include()]),
        Extension('aipy._dsp', ['aipy_src/_dsp/dsp.c', 'aipy_src/_dsp/grid/grid.c'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy_src/_dsp', 'aipy_src/_dsp/grid', 'aipy_src/_common']),
        Extension('aipy.utils', ['aipy_src/utils/utils.cpp'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy_src/_common']),
    ],
    scripts=glob.glob('scripts/*'),
    package_data = {'aipy': ['aipy_src/_src/*.txt']},
    include_package_data = True,
    zip_safe = False,
    test_suite = "tests.aipy_test.TestSuite",
)
