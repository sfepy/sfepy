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
    config.add_data_files(('sfepy/fem/extmods', ('version.h.in',)))

    defines = [('__SDIR__', "'\"%s\"'" % auto_dir),
               ('SFEPY_PLATFORM', os_flag[site_config.system()])]
    if '-DDEBUG_FMF' in site_config.debug_flags():
        defines.append(('DEBUG_FMF', None))

    src = ['fmfield.c', 'geometry.c', 'geometry.i', 'geommech.c',
           'common_python.c']
    config.add_extension('_geometry',
                         sources=src,
                         depends=['array.i', 'common.i', 'fmfield.i'],
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    common_src = ['fmfield.c', 'geometry.c', 'geommech.c', 'common_python.c']

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

    src = ['bases.pyx']
    config.add_extension('bases',
                         sources=src,
                         libraries=['sfepy_common'],
                         depends=common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['mesh.pyx', 'geomtrans.c', 'meshutils.c', 'sort.c',
           'common_python.c']
    config.add_extension('mesh',
                         sources=src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
