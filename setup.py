def configuration(parent_package='', top_path=None):
    import glob
    from numpy.distutils.misc_util import Configuration
    config = Configuration('aipy', parent_package, top_path,
        version='0.4.0',
        author='Aaron Parsons',
        author_email='aparsons at astron.berkeley.edu',
        url='http://setiathome.berkeley.edu/~aparsons',
        license='GPL',
        description='Astronomical Interferometry in PYthon',
        long_description=
"""Tools for radio astronomical interferometry, including 
pure-python phasing, calibration, imaging, and deconvolution code, 
interfaces to MIRIAD, HEALPix, fitting routines from SciPy, and 
the PyFITS and PyEphem packages verbatim.""",
    )
    config.add_subpackage('optimize')
    config.add_subpackage('interpolate')
    config.add_subpackage('pyephem')
    config.add_subpackage('pyfits')
    config.add_subpackage('miriad')
    config.add_subpackage('healpix')
    config.add_subpackage('utils')
    config.add_data_dir('data')
    #config.add_data_dir('doc')
    config.add_scripts(glob.glob('scripts/*'))
    config.add_data_files(('.', 'LICENSE.txt'))
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
