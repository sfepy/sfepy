#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    import os.path as op
    from numpy.distutils.misc_util import Configuration

    from sfepy import Config

    site_config = Config()
    os_flag = {'posix' : 0, 'windows' : 1}

    auto_dir = op.dirname(__file__)
    auto_name = op.split(auto_dir)[-1]
    config = Configuration(auto_name, parent_package, top_path)

    defines = [('__SDIR__', "'\"%s\"'" % auto_dir),
               ('SFEPY_PLATFORM', os_flag[site_config.system()])]
    if '-DDEBUG_FMF' in site_config.debug_flags():
        defines.append(('DEBUG_FMF', None))

    common_src = ['fmfield.c', 'refmaps.c', 'geommech.c', 'common_python.c']

    config.add_library('sfepy_common',
                       sources=common_src,
                       extra_compiler_args=site_config.compile_flags(),
                       extra_link_args=site_config.link_flags(),
                       include_dirs=[auto_dir, site_config.python_include()],
                       macros=defines)

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

    src = ['bases.pyx', 'lagrange.c']
    config.add_extension('bases',
                         sources=src,
                         libraries=['sfepy_common'],
                         depends=common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['cmesh.pyx', 'geomtrans.c', 'meshutils.c', 'sort.c',
           'common_python.c']
    config.add_extension('cmesh',
                         sources=src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['lobatto_bases.pyx', 'lobatto.c', 'lobatto1d.c']
    config.add_extension('lobatto_bases',
                         sources=src,
                         libraries=['sfepy_common'],
                         depends=common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    # Include *.pxd files in distribution tarball and install them along
    # with the extension modules.
    pxd_files = ['mappings.pxd', 'types.pxd', '_fmfield.pxd']
    config.add_data_files(('', pxd_files))

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
