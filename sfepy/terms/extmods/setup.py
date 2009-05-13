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

    fem_src = ['common_python.c', 'fmfield.c', 'geommech.c']
    fem_src = [op.join('../../fem/extmods', ii) for ii in fem_src]
    src = ['formSDCC.c', 'terms.c', 'termsBasic.c', 'termsElectric.c',
           'termsMass.c', 'termsNavierStokes.c', 'termsBiot.c', 'termsLaplace.c',
           'termsLinElasticity.c', 'termsHyperElasticity.c', 'termsPiezo.c',
           'termsSurface.c', 'termsVolume.c', 'terms.i']
    config.add_extension('_terms',
                         sources=src + fem_src,
                         depends=[],
                         include_dirs=[auto_dir, '../../fem/extmods'],
                         define_macros=defines)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
