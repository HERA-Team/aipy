def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    import numpy

    config = Configuration('miriad', parent_package, top_path)
    # To update the python package to the latest MIRIAD release, you need to
    # copy files from $MIRINC and $MIRSRC/subs to ./mirsrc
    config.add_include_dirs(['mirsrc'])
    config.add_extension('_miruv',
        sources=['miruv.i', 'miriad.py', 'wrap_miruv_swig.c', 
            'mirsrc/uvio.c', 'mirsrc/hio.c', 'mirsrc/pack.c', 'mirsrc/bug.c', 
            'mirsrc/dio.c', 'mirsrc/headio.c', 'mirsrc/maskio.c'],
        include_dirs=['./', 'mirsrc', numpy.get_include(), 
            numpy.__path__[0]+'/doc/swig'],
    )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
