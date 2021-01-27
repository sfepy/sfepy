#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    import os.path as op
    from numpy.distutils.misc_util import Configuration

    from sfepy import Config

    site_config = Config()
    system = site_config.system()
    os_flag = {'posix' : 0, 'windows' : 1}[system]

    auto_dir = op.dirname(__file__)
    auto_name = op.split(auto_dir)[-1]
    config = Configuration(auto_name, parent_package, top_path)

    inline = 'inline' if system == 'posix' else '__inline'
    defines = [('SFEPY_PLATFORM', os_flag),
               ('inline', inline)]
    if '-DDEBUG_FMF' in site_config.debug_flags():
        defines.append(('DEBUG_FMF', None))
    if '-DDEBUG_MESH' in site_config.debug_flags():
        defines.append(('DEBUG_MESH', None))

    common_src = ['fmfield.c', 'refmaps.c', 'geommech.c', 'common_python.c']

    config.add_library('sfepy_common',
                       sources=common_src,
                       extra_compiler_args=site_config.compile_flags(),
                       extra_link_args=site_config.link_flags(),
                       include_dirs=[auto_dir, site_config.python_include()],
                       macros=[('SFEPY_PLATFORM', os_flag),
                               ('inline', inline)])

    src = ['_fmfield.pyx']
    config.add_extension('_fmfield',
                         sources=src,
                         libraries=['sfepy_common'],
                         depends=common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['mappings.pyx']
    config.add_extension('mappings',
                         sources=src,
                         libraries=['sfepy_common'],
                         depends=common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['assemble.pyx']
    config.add_extension('assemble',
                         sources=src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['cmesh.pyx', 'geomtrans.c', 'mesh.c', 'meshutils.c', 'sort.c',
           'common_python.c']
    config.add_extension('cmesh',
                         sources=src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['crefcoors.pyx', 'refcoors.c', 'geomtrans.c', 'mesh.c']
    config.add_extension('crefcoors',
                         sources=src,
                         libraries=['sfepy_common'],
                         depends=common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['_geommech.pyx']
    config.add_extension('_geommech',
                         sources=src,
                         libraries=['sfepy_common'],
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    # Include *.pxd files in distribution tarball and install them along
    # with the extension modules.
    pxd_files = ['cmesh.pxd', 'mappings.pxd', 'types.pxd',
                 '_fmfield.pxd', '_geommech.pxd']
    config.add_data_files(('', pxd_files))

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
