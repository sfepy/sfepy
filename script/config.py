#! /usr/bin/env python
"""
SfePy configuration script. It is used in the main Makefile to detect the
correct python paths and other options that can be specified in
site_cfg.py. Prints results to stdout.

options:

python_version

    The default Python version that is installed on the system.

system

    The operating system (posix or windows).

compile_flags

    Extra compile flags added to the flags supplied by distutils.

link_flags

    Extra linker flags added to the flags supplied by distutils.

debug_flags

    Debugging flags.

numpydoc_path

    The path to numpydoc (required for the sphinx documentation).

is_release

    If set, the version is a release.

tetgen_path

    Tetgen executable path.

New options should be added both to site_cfg_template.py and Config class below.

Examples:

$ ./config.py python_version
2.5
$

$ ./config.py archlib
lib
$
"""
import sys
sys.path.append( '.' )
import os
import shutil

msg_unknown_os = """could not determine operating system!
try setting it in site_cfg.py manually, see site_cfg_template.py"""

msg_numpydoc = """could not find numpydoc!
If it is installed in a non-standard location, try setting it in
site_cfg.py manually."""

try:
    import site_cfg
except:
    try:
        shutil.copyfile( 'site_cfg_template.py', 'site_cfg.py' )
        import site_cfg
    except:
        site_cfg = None

has_attr = lambda obj, attr: obj and hasattr( obj, attr )

def find_in_path( filename ):
    pathlist=sys.path
    outpath=''
    for curpath in pathlist:
        temppath=os.path.join(curpath,filename)
        if os.path.exists(temppath):
            outpath=temppath
            break
    return outpath 

class Config( object ):
    def python_version( self ):
        if has_attr( site_cfg, 'python_version' ):
            if site_cfg.python_version == 'auto':
                return "%d.%d" % tuple(sys.version_info[:2])
            else:
                return site_cfg.python_version
        else:
            return "%d.%d" % tuple(sys.version_info[:2])

    def system( self ):
        if has_attr( site_cfg, 'system' ) and site_cfg.system is not None:
            return site_cfg.system
        else:
            if os.name in ['posix']:
                return 'posix'
            elif os.name in ['nt']:
                return 'windows'
            else:
                raise ValueError(msg_unknown_os)

    def compile_flags( self ):
        if has_attr(site_cfg, 'compile_flags'):
            flags = site_cfg.compile_flags

        else:
            flags = '-g -O2'

        return flags.split()

    def link_flags( self ):
        if has_attr( site_cfg, 'link_flags' ):
            flags =  site_cfg.link_flags

        else:
            flags = ''

        return flags.split()

    def debug_flags( self ):
        if has_attr( site_cfg, 'debug_flags' ):
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

    def tetgen_path( self ):
        if has_attr( site_cfg, 'tetgen_path' ):
            return site_cfg.tetgen_path
        else:
            return '/usr/bin'

usage = """Usage: %s option"""

def main():
    try:
        mode = sys.argv[1]
    except:
        print usage % sys.argv[0]
        return

    config = Config()
    try:
        print getattr( config, mode )()
    except:
        raise

if __name__ == '__main__':
    main()

