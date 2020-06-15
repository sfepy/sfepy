import sys
import os
import shutil
import sysconfig

msg_unknown_os = """could not determine operating system!
try setting it in site_cfg.py manually, see site_cfg_template.py"""

msg_numpydoc = """could not find numpydoc!
If it is installed in a non-standard location, try setting it in
site_cfg.py manually."""

if not os.path.exists('site_cfg.py'):
    try:
        shutil.copyfile('site_cfg_template.py', 'site_cfg.py')

    except:
        pass

try:
    import site_cfg

except ImportError:
    site_cfg = None

has_attr = lambda obj, attr: obj and hasattr(obj, attr)

class Config(object):
    def python_version(self):
        if has_attr(site_cfg, 'python_version'):
            if '*' in site_cfg.python_version:
                return "%d.%d" % tuple(sys.version_info[:2])
            else:
                return site_cfg.python_version
        else:
            return "%d.%d" % tuple(sys.version_info[:2])

    def python_include(self):
        if (has_attr(site_cfg, 'python_include')
            and (site_cfg.python_include != 'auto')):
            return site_cfg.python_include

        else:
            return sysconfig.get_config_var('INCLUDEPY')


    def system(self):
        if has_attr(site_cfg, 'system') and site_cfg.system is not None:
            return site_cfg.system
        else:
            if os.name in ['posix']:
                return 'posix'
            elif os.name in ['nt']:
                return 'windows'
            else:
                raise ValueError(msg_unknown_os)

    def compile_flags(self):
        if has_attr(site_cfg, 'compile_flags'):
            flags = site_cfg.compile_flags

        else:
            flags = '-g -O2'

        return flags.split()

    def link_flags(self):
        if has_attr(site_cfg, 'link_flags'):
            flags =  site_cfg.link_flags

        else:
            flags = ''

        return flags.split()

    def debug_flags(self):
        if has_attr(site_cfg, 'debug_flags'):
            return site_cfg.debug_flags
        else:
            return ''

    def numpydoc_path(self):
        if (has_attr(site_cfg, 'numpydoc_path') and
            (site_cfg.numpydoc_path is not None)):
            return site_cfg.numpydoc_path

        else:
            try:
                import numpydoc
            except ImportError:
                raise ValueError(msg_numpydoc)

    def is_release(self):
        if has_attr(site_cfg, 'is_release'):
            return site_cfg.is_release
        else:
            return ''

    def tetgen_path(self):
        if has_attr(site_cfg, 'tetgen_path'):
            return site_cfg.tetgen_path
        else:
            return '/usr/bin/tetgen'

    def refmap_memory_factor(self):
        if has_attr(site_cfg, 'refmap_memory_factor'):
            return site_cfg.refmap_memory_factor
        else:
            return None
