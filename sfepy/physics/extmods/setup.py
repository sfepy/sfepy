#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    import distutils.sysconfig as sysc
    from numpy.distutils.misc_util import Configuration
    import os.path as op

    auto_dir = op.dirname(__file__)
    auto_name = op.split(auto_dir)[-1]
    config = Configuration(auto_name, parent_package, top_path)

    defines = [('__SDIR__', "'\"%s\"'" % auto_dir),
               ('DEBUGFMF', None)]

    src = ['dft.c', 'dft.i']
    config.add_extension('_dft',
                         sources=src,
                         depends=[],
                         include_dirs=[auto_dir, '../../fem/extmods'],
                         define_macros=defines)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
