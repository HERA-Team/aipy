def configuration(parent_package='', top_path=None):
    import glob
    from numpy.distutils.misc_util import Configuration
    config = Configuration('aipy', parent_package, top_path,
        version='0.3.0',
        author='Aaron Parsons',
        author_email='aparsons at astron.berkeley.edu',
        url='http://setiathome.berkeley.edu/~aparsons',
        license='GPL',
        description='Astronomical Interferometry in PYthon',
    )
    config.add_subpackage('optimize')
    config.add_subpackage('interpolate')
    config.add_subpackage('special')
    config.add_subpackage('pyephem')
    config.add_subpackage('miriad')
    config.add_data_dir('data')
    #config.add_data_dir('doc')
    config.add_scripts(glob.glob('scripts/*'))
    config.add_data_files(('.', 'LICENSE.txt'))
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
