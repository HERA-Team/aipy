libastro_version = '3.7.2'

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyephem', parent_package, top_path)
    config.add_extension('ephem',
        sources=['libastro-%s/*.c' % libastro_version, 'ephem.c'],
        include_dirs=['libastro-' + libastro_version],
    )
    config.add_data_files(('.', ['COPYING', 'README']))
    return config

if __name__ == '__main__':
    from distutils.core import setup
    setup(configuration=configuration)
