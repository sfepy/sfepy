#!/usr/bin/env python
from optparse import OptionParser

import sfepy
from sfepy.base.base import pause
from sfepy.postprocess import Viewer

usage = """%prog [options] filename"""

help = {
    'is_3d' :
    '3d plot mode',
}

def main():
    parser = OptionParser(usage=usage, version="%prog " + sfepy.__version__)
    parser.add_option( "--3d",
                       action = "store_true", dest = "is_3d",
                       default = False, help = help['is_3d'] )
    options, args = parser.parse_args()

    if (len(args) == 1):
        filename = args[0]
    else:
        parser.print_help(),
        return

    view = Viewer(filename)
    view(is_3d=options.is_3d)

if __name__ == '__main__':
    main()
