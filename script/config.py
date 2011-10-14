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

    True for a release, False otherwise. If False, current git commit hash
    is appended to version string, if the sources are in a repository.

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
sys.path.append('.')

from sfepy import Config

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

