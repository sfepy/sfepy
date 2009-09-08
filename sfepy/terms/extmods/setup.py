#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    import distutils.sysconfig as sysc
    from numpy.distutils.misc_util import Configuration
    import os.path as op
    import glob

    import sys;
    if 'script' not in sys.path:
        sys.path.append('script')
    from config import Config
    system = Config().system()
    os_flag = {'posix' : 0, 'windows' : 1}

    auto_dir = op.dirname(__file__)
    auto_name = op.split(auto_dir)[-1]
    config = Configuration(auto_name, parent_package, top_path)

    defines = [('__SDIR__', "'\"%s\"'" % auto_dir),
               ('DEBUGFMF', None),
               ('SFEPY_PLATFORM', os_flag[system])]

    fem_src = ['common_python.c', 'fmfield.c', 'geommech.c']
    fem_src = [op.join('../../fem/extmods', ii) for ii in fem_src]

    src = [op.split(ii)[1] for ii in glob.glob('sfepy/terms/extmods/*.c')]
    src += ['terms.i']
    try:
        src.remove('terms_wrap.c')
    except ValueError:
        pass

    depends=['array.i', 'common.i', 'fmfield.i']
    depends = [op.join('../../fem/extmods', ii) for ii in depends]

    config.add_extension('_terms',
                         sources=src + fem_src,
                         depends=depends,
                         extra_compile_args=['-O2'],
                         include_dirs=[auto_dir, '../../fem/extmods'],
                         define_macros=defines)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
