import sys

def configuration(parent_package='',top_path=None):
    if not hasattr(sys, 'version_info') or sys.version_info < (2,3,0,'alpha',0):
        raise SystemExit, "Python 2.3 or later required to build pyfits."
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyfits', parent_package, top_path)
    config.add_extension('pyfits',
        sources=['*.py'],
    )
    config.add_data_files(('.', ['LICENSE.txt']))
    return config

if __name__ == '__main__':
    from distutils.core import setup
    setup(configuration=configuration)
        
