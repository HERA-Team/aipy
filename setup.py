def configuration(parent_package='', top_path=None):
    import glob
    from numpy.distutils.misc_util import Configuration
    config = Configuration('aipy', parent_package, top_path,
        version='0.2.2')
    config.add_subpackage('lbfgsb')
    config.add_subpackage('miriad')
    config.add_data_dir('data')
    config.add_data_dir('doc')
    config.add_scripts(glob.glob('scripts/*'))
    config.add_data_files(('.', 'LICENSE.txt'))
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
