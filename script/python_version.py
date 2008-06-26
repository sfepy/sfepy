#! /usr/bin/env python

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

import sys
ver = "%d.%d" % tuple(sys.version_info[:2])
print ver
