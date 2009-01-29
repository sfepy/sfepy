#!/usr/bin/env python
import sys
sys.path.append('.')

from optparse import OptionParser
from sfepy.fem import Mesh

usage = """%prog [options] filename_in filename_out

Convert a mesh file from one SfePy-supported format to another.
"""

def main():
    parser = OptionParser( usage = usage )
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        sys.exit(1)

    filename_in, filename_out = args
    
    mesh = Mesh.from_file(filename_in)
    mesh.write(filename_out, io='auto')

if __name__ == '__main__':
    main()
