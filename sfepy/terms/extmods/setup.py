#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    import os.path as op
    import sys
    import glob
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

    fem_src = ['common_python.c', 'fmfield.c', 'geommech.c']
    fem_src = [op.join('../../fem/extmods', ii) for ii in fem_src]

    csrc = [op.split(ii)[1] for ii in glob.glob('sfepy/terms/extmods/*.c')]
    try:
        csrc.remove('terms_wrap.c')
    except ValueError:
        pass

    common_src = ['fmfield.c', 'geommech.c', 'common_python.c']
    common_src = [op.join('../../fem/extmods', ii) for ii in common_src]

    src = ['terms.i']

    depends=['array.i', 'common.i', 'fmfield.i']
    depends = [op.join('../../fem/extmods', ii) for ii in depends]

    config.add_library('sfepy_terms',
                       sources=csrc,
                       depends=common_src,
                       extra_compiler_args=site_config.compile_flags(),
                       extra_link_args=site_config.link_flags(),
                       include_dirs=[auto_dir, '../../fem/extmods',
                                     site_config.python_include()],
                       macros=defines)

    config.add_extension('_terms',
                         sources=src,
                         libraries=['sfepy_terms', 'sfepy_common'],
                         depends=depends + csrc + common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir, '../../fem/extmods'],
                         define_macros=defines)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
