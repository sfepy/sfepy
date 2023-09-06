import sys
import os
import shutil
import sysconfig
from warnings import warn

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


def compose_system_compile_flags(is_posix: bool) -> list:
    """
    Provides a list of compile flags that tries to be as similar as possible
    to what Python itself was built with. This has been done historically by
    numpy.distutils (now deprecated) and a squeezed version of it is brought
    over to here.
    """
    if not is_posix:
        return []

    cflags, configure_cppflags, configure_cflags = sysconfig.get_config_vars(
        'CFLAGS', 'CONFIGURE_CPPFLAGS', 'CONFIGURE_CFLAGS')

    return (cflags + ' ' + configure_cppflags + ' ' + configure_cflags).split()


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
            if isinstance(flags, str):
                warn('Compile flags should be given as a list of strings.'
                     ' Space-separated strings may be removed in the near future.',
                     DeprecationWarning, stacklevel=2)

                flags = flags.split()

        else:
            flags = ['-g', '-O2']

        return flags + compose_system_compile_flags(self.system() == 'posix')

    def debug_flags(self) -> list:
        if has_attr(site_cfg, 'debug_flags'):
            return site_cfg.debug_flags
        else:
            return []

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
