def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os, numpy

    env = os.environ
    if not env.has_key('MIRSRC') or not env.has_key('MIRLIB'):
        raise Exception('Please source your miriad environment.')

    config = Configuration('miriad', parent_package, top_path)
    config.add_include_dirs([env['MIRSRC']+'/subs'])
    config.add_extension('_miruv',
        sources=['miruv.i', 'miriad.py', 'wrap_miruv_swig.c'],
        include_dirs=['./', env['MIRSRC']+'/subs', numpy.get_include(), 
            numpy.__path__[0]+'/doc/swig'],
        library_dirs=[env['MIRLIB']],
        extra_objects=[env['MIRLIB']+'/libmir.a'],
        libraries=['m'],
    )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
