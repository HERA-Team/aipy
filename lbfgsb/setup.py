def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('lbfgsb', parent_package, top_path)
    config.add_extension('_lbfgsb',
        sources = ['routines.pyf', 'lbfgsb.py', 'routines.f'],
    )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
