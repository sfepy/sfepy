#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    import os.path as op
    import sys
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

    src = ['fmfield.c', 'fem.c', 'fem.i', 'geommech.c', 'sort.c',
           'common_python.c']
    config.add_extension('_fem',
                         sources=src,
                         depends=['array.i', 'common.i', 'fmfield.i'],
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['fmfield.c', 'geometry.c', 'geometry.i', 'geommech.c',
           'common_python.c']
    config.add_extension('_geometry',
                         sources=src,
                         depends=['array.i', 'common.i', 'fmfield.i'],
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
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['mesh.pyx', 'geomtrans.c', 'meshutils.c', 'common_python.c']
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
