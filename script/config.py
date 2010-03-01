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

archlib

    'lib' or 'lib64' depending on your architecture (32bit or 64bit)

tetgen_path

    Tetgen executable path.

numpy_include

    Full path to the numpy headers.  This path should end in numpy/core/include.

opt_flags

    Compiler optimization flags.

link_flags

    Linker flags.

debug_flags

    Debugging flags.

numpydoc_path

    The path to numpydoc (required for the sphinx documentation).

is_release

    If set, the version is a release.

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

    def archlib( self ):
        if has_attr( site_cfg, 'archlib' ):
            return site_cfg.archlib
        else:
            return 'lib'

    def tetgen_path( self ):
        if has_attr( site_cfg, 'tetgen_path' ):
            return site_cfg.tetgen_path
        else:
            return '/usr/bin'

    def numpy_include( self ):
        if has_attr( site_cfg, 'numpy_include' ) and (site_cfg.numpy_include is not None):
            return site_cfg.numpy_include
        else:
            numpypath = find_in_path( 'numpy' )
            numpyfullpath = os.path.join( numpypath, 'core', 'include' )
            if not os.path.exists(numpyfullpath):
                print('could not find core/include in '+ numpypath)
            return numpyfullpath

    def opt_flags( self ):
        if has_attr( site_cfg, 'opt_flags' ):
            return site_cfg.opt_flags
        else:
            return '-g -O2 -fPIC -DPIC'

    def link_flags( self ):
        if has_attr( site_cfg, 'link_flags' ):
            return site_cfg.link_flags
        else:
            return '-shared -fPIC -DPIC'

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

