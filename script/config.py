#! /usr/bin/env python
import sys
sys.path.append( '.' )
import shutil

try:
    import site_cfg
except:
    shutil.copyfile( 'site_cfg_template.py', 'site_cfg.py' )
    try:
        import site_cfg
    except:
        site_cfg = None

has_attr = lambda obj, attr: obj and hasattr( obj, attr )

class Config( object ):
    def python_version( self ):
        """
        Detects the default Python version that is installed on the system and prints
        it to stdout. This is used in the main Makefile to detect the correct python
        paths.

        Examples:

            $ ./python_version.py
            2.5
            $

        on some other system:

            $ ./python_version.py
            2.4
            $

        """
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
