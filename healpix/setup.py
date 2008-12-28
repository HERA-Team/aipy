def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    import numpy
    def indir(dir, files): return [dir+f for f in files]
    src_files = ['healpix_base_wrap.cpp']
    src_files += indir('cxx/Healpix_cxx/', ['healpix_base.cc'])
    config = Configuration('healpix', parent_package, top_path)
    config.add_extension('_healpix_base',
        sources=src_files + ['healpix.py'],
        include_dirs=[numpy.get_include(), 'cxx/cxxsupport', 'cxx/Healpix_cxx'],
        #libraries=['m'],
    )
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
