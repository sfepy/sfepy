#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    import distutils.sysconfig as sysc
    from numpy.distutils.misc_util import Configuration
    import os.path as op
##     from numpy.distutils.command import build_src
##     import Cython
##     import Cython.Compiler.Main
##     build_src.Pyrex = Cython
##     build_src.have_pyrex = True

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

    src = ['rcm.c', 'rcm.i', 'common_python.c']
    config.add_extension('_rcm',
                         sources=src,
                         depends=['array.i', 'common.i'],
                         extra_compile_args=['-O2'],
                         include_dirs=[auto_dir],
                         define_macros=defines)

##     src = ['ct.pyx']
##     config.add_extension('_ct',
##                          sources=src,
##                          extra_compile_args=['-O2'],
##                          include_dirs=[auto_dir],
##                          define_macros=defines)
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
