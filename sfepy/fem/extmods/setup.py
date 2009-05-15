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

    src = ['fmfield.c', 'fem.c', 'fem.i', 'geommech.c', 'sort.c',
           'common_python.c']
    config.add_extension('_fem',
                         sources=src,
                         depends=[],
                         extra_compile_args=['-O2'],
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['geomtrans.c', 'meshutils.c', 'meshutils.i', 'common_python.c',
           'sort.c']
    config.add_extension('_meshutils',
                         sources=src,
                         depends=[],
                         extra_compile_args=['-O2'],
                         include_dirs=[auto_dir],
                         define_macros=defines)

    src = ['fmfield.c', 'geometry.c', 'geometry.i', 'geommech.c',
           'common_python.c']
    config.add_extension('_geometry',
                         sources=src,
                         depends=[],
                         extra_compile_args=['-O2'],
                         include_dirs=[auto_dir],
                         define_macros=defines)

    
## gcc -g -O2 -fPIC -DPIC -D__SDIR__='"tests"' -DDEBUG_FMF -Iexamples -Iinput
## -Isfepy -Isfepy/applications -Isfepy/base -Isfepy/eldesc -Isfepy/fem
## -Isfepy/fem/extmods -Isfepy/geom -Isfepy/homogenization -Isfepy/mechanics
## -Isfepy/solvers -Isfepy/terms -Isfepy/terms/extmods -Isfepy/physics
## -Isfepy/physics/extmods -Itests -I/usr/include/python2.5
## -I/home/share/software/usr/lib/python/site-packages/numpy/core/include -Wall
## -c -DISRELEASE sfepy/fem/extmods/fmfield.c -o sfepy/fem/extmods/fmfield.o

## swig -python -Iexamples -Iinput -Isfepy -Isfepy/applications -Isfepy/base
## -Isfepy/eldesc -Isfepy/fem -Isfepy/fem/extmods -Isfepy/geom
## -Isfepy/homogenization -Isfepy/mechanics -Isfepy/solvers -Isfepy/terms
## -Isfepy/terms/extmods -Isfepy/physics -Isfepy/physics/extmods -Itests
## -I/usr/include/python2.5
## -I/home/share/software/usr/lib/python/site-packages/numpy/core/include
## -DISRELEASE -o sfepy/fem/extmods/meshutils_wrap.c
## sfepy/fem/extmods/meshutils.i

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
