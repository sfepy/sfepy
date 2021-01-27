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

    common_path = '../../common/extmods'

    common_src = ['fmfield.c', 'geommech.c', 'common_python.c']
    common_src = [op.join(common_path, ii) for ii in common_src]

    src = ['igac.pyx', 'nurbs.c']
    config.add_extension('igac',
                         sources=src + common_src,
                         depends=common_src,
                         extra_compile_args=site_config.compile_flags(),
                         extra_link_args=site_config.link_flags(),
                         include_dirs=[auto_dir, common_path],
                         define_macros=defines)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
