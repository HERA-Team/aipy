#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

from setuptools import setup, Extension

import os, glob, numpy, subprocess, sys

PY2 = sys.version_info.major < 3
if PY2:
    MATPLOTLIB_DEP = 'matplotlib<3'
    ASTROPY_DEP = 'astropy>=1.0'
else:
    MATPLOTLIB_DEP = 'matplotlib'
    ASTROPY_DEP = 'astropy>=3.0'

def get_description():
    def get_description_lines():
        seen_desc = False

        with open('README.md') as f:
            for line in f:
                if seen_desc:
                    if line.startswith('##'):
                        break
                    line = line.strip()
                    if len(line):
                        yield line
                elif line.startswith('## Description'):
                    seen_desc = True

    return ' '.join(get_description_lines())

def indir(path, files):
    return [os.path.join(path, f) for f in files]

global_macros = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]

setup(
    name = 'aipy',
    use_setuptools_scm={
        "write_to": "aipy/_version.py",
        "parentdir_prefix_version": "aipy-",
        "fallback_version": "0.0.0",
    },
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
    setup_requires = [
        'numpy>=1.2',
        'setuptools_scm',
    ],
    install_requires = [
        ASTROPY_DEP,
        'astropy-healpix',
        MATPLOTLIB_DEP,
        'numpy>=1.2',
        'ephem>=3.7.3.2',
        'scipy>=0.19',
    ],
    extras_require = {
        'dev': [
            'pytest',
            'pytest-cov'
        ]
    },
    package_dir = {'aipy':'aipy', 'aipy._src':'aipy/_src'},
    packages = ['aipy', 'aipy._src'],
    ext_modules = [
        Extension('aipy._miriad', ['aipy/_miriad/miriad_wrap.cpp'] + \
            indir('aipy/_miriad/mir', ['uvio.c','hio.c','pack.c','bug.c',
                'dio.c','headio.c','maskio.c']),
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy/_miriad',
                'aipy/_miriad/mir', 'aipy/_common']),
        Extension('aipy._deconv', ['aipy/_deconv/deconv.cpp'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy/_common']),
        #Extension('aipy._img', ['aipy/_img/img.cpp'],
        #    include_dirs = [numpy.get_include()]),
        Extension('aipy._dsp', ['aipy/_dsp/dsp.c', 'aipy/_dsp/grid/grid.c'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy/_dsp', 'aipy/_dsp/grid', 'aipy/_common']),
        Extension('aipy.utils', ['aipy/utils/utils.cpp'],
            define_macros = global_macros,
            include_dirs = [numpy.get_include(), 'aipy/_common']),
    ],
    scripts = glob.glob('scripts/*'),

    include_package_data = True,
    zip_safe = False,
    test_suite = "tests.aipy_test.TestSuite",
)
