#! /usr/bin/env python
"""
SfePy configuration script. It is used in the main Makefile to detect the
correct python paths and other options that can be specified in
site_cfg.py. Prints results to stdout.

options:

python_version

    The default Python version that is installed on the system.

archlib

    'lib' or 'lib64' depending on your architecture (32bit or 64bit)

tetgen_path

    Tetgen executable path.

numpy_prefix

    Installation prefix. Ignore if you do not install numpy/scipy on local
    account.  Numpy headers should be in
    <numpy_prefix>/usr/<archlib>/python<python_version>/site-packages/numpy/core/include

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
import shutil

try:
    import site_cfg
except:
    try:
        shutil.copyfile( 'site_cfg_template.py', 'site_cfg.py' )
        import site_cfg
    except:
        site_cfg = None

has_attr = lambda obj, attr: obj and hasattr( obj, attr )

class Config( object ):
    def python_version( self ):
        if has_attr( site_cfg, 'python_version' ):
            if site_cfg.python_version == 'auto':
                return "%d.%d" % tuple(sys.version_info[:2])
            else:
                return site_cfg.python_version
        else:
            return "%d.%d" % tuple(sys.version_info[:2])
    
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

    def numpy_prefix( self ):
        if has_attr( site_cfg, 'numpy_prefix' ):
            return site_cfg.numpy_prefix
        else:
            return '/'
    
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
        pass

if __name__ == '__main__':
    main()
