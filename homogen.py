#!/usr/bin/env python
from __future__ import absolute_import
from argparse import ArgumentParser

import sfepy
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.homogenization.homogen_app import HomogenizationApp

helps = {
    'debug':
    'automatically start debugger when an exception is raised',
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
}

def main():
    parser = ArgumentParser()
    parser.add_argument("--version", action="version",
                        version="%(prog)s " + sfepy.__version__)
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        default=False, help=helps['debug'])
    parser.add_argument("-o", metavar='filename', action="store",
                        dest="output_filename_trunk",
                        default=None, help=helps['filename'])
    parser.add_argument('filename_in')
    options = parser.parse_args()

    if options.debug:
        from sfepy.base.base import debug_on_error; debug_on_error()

    filename_in = options.filename_in

    required, other = get_standard_keywords()
    required.remove('equations')

    conf = ProblemConf.from_file(filename_in, required, other)

    app = HomogenizationApp(conf, options, 'homogen:')
    opts = conf.options
    if hasattr(opts, 'parametric_hook'):  # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
