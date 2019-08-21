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

print("Generating aipy/__version__.py: ", end='')
__version__ = open('VERSION').read().strip()
print(__version__)
open('aipy/__version__.py','w').write('__version__="%s"'%__version__)

#read the latest git status out to an installed file
try:
#    gitbranch = subprocess.check_output('git symbolic-ref -q HEAD',shell=True, cwd='.').strip().split('/')[-1]
    gitbranch = os.popen('git symbolic-ref -q HEAD').read().strip()
    print("Generating aipy/__branch__.py")
#    gitlog = subprocess.check_output('git log -n1 --pretty="%h%n%s%n--%n%an%n%ae%n%ai"',shell=True, cwd='.').strip()
    gitlog = os.popen('git log -n1 --pretty="%h%n%s%n--%n%an%n%ae%n%ai"').read().strip()
    print("Generating aipy/__gitlog__.py.")
    print(gitlog)
except:
    gitbranch = "unknown branch"
    gitlog = "git log not found"
open('aipy/__branch__.py','w').write('__branch__ = \"%s\"'%gitbranch)
open('aipy/__gitlog__.py','w').write('__gitlog__ = \"\"\"%s\"\"\"'%gitlog)


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

    setup_requires = [
        'numpy>=1.2'
    ],

    install_requires = [
        ASTROPY_DEP,
        'healpy>=1.11',
        MATPLOTLIB_DEP,
        'numpy>=1.2',
        'pyephem>=3.7.3.2',
        'scipy>=0.19',
    ],

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
