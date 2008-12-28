import ez_setup, os, glob, numpy
ez_setup.use_setuptools()

from setuptools import setup, Extension

__version__ = '0.4.2'

def indir(path, files):
    return [os.path.join(path, f) for f in files]

#class NumpyExtension(Extension):
#    """Clever workaround for waiting until numpy dependency is resolved
#    before importing the numpy include structure needed to build extension"""
#    def __init__(self, *args, **kwargs):
#        Extension.__init__(self, *args, **kwargs)
#        self._include_dirs = self.include_dirs
#        del self.include_dirs
#    @property
#    def include_dirs(self):
#        import numpy
#        return self._include_dirs + [numpy.get_include()]

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
    install_requires = ['pyephem>=3.7.2.3', 'pyfits>=1.1',
        'matplotlib>=0.91'], #'basemap==0.9.1'],
    dependency_links = [
        'http://www.stsci.edu/resources/software_hardware/pyfits/pyfits-1.1.tar.gz'],
    packages = ['aipy', 'aipy.optimize'],
    ext_modules = [
        #NumpyExtension('aipy._healpix',
        Extension('aipy._healpix',
            ['aipy/_healpix/healpix_wrap.cpp',
            'aipy/_healpix/cxx/Healpix_cxx/healpix_base.cc'],
            include_dirs = [numpy.get_include(), 'aipy/_healpix/cxx/cxxsupport',
                'aipy/_healpix/cxx/Healpix_cxx']),
        #NumpyExtension('aipy._miriad', ['aipy/_miriad/miriad_wrap.cpp'] + \
        Extension('aipy._miriad', ['aipy/_miriad/miriad_wrap.cpp'] + \
            indir('aipy/_miriad/mir', ['uvio.c','hio.c','pack.c','bug.c',
                'dio.c','headio.c','maskio.c']),
            include_dirs = [numpy.get_include(), 'aipy/_miriad', 'aipy/_miriad/mir']),
        #NumpyExtension('aipy.utils', ['aipy/utils/utils.cpp'],
        Extension('aipy.utils', ['aipy/utils/utils.cpp'],
            include_dirs = [numpy.get_include()])
    ],
    scripts=glob.glob('scripts/*'),
    package_data = {'aipy': ['doc/*.tex', 'doc/*.png']},
)
