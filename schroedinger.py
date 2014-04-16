#!/usr/bin/env python
"""
Solver for basic electronic structure problems.

Examples
--------

Compute the 2D predefined problems::

  $ ./schroedinger.py --well
  $ ./schroedinger.py --oscillator
  $ ./schroedinger.py --hydrogen
  $ ./schroedinger.py --boron

See ``examples/quantum/`` directory for the problem description files
corresponding to the above examples.
"""
import os
from optparse import OptionParser

import sfepy
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.physics.schroedinger_app import SchroedingerApp

def fix_path(filename):
    return os.path.join(sfepy.data_dir, filename)

usage = """%prog [options] [filename_in]\n""" + __doc__.rstrip()

help = {
    'conf' :
    'override problem description file items, written as python'
    ' dictionary without surrounding braces',
    'options' : 'override options item of problem description,'
    ' written as python dictionary without surrounding braces',
    'filename' :
    'basename of output file(s) [default: <basename of input file mesh>]',
    'well' :
    'solve infinite potential well (particle in a box) problem',
    'oscillator' :
    'solve spherically symmetric linear harmonic oscillator '
    '(1 electron) problem',
    'hydrogen' :
    'solve the hydrogen atom',
    'boron' :
    'solve the boron atom with 1 electron',
    'n_eigs' :
    'number of eigenpairs to compute [default: as set in the examples]',
    'tau' :
    'target value of the Pysparse eigenvalue solver. Eigenvalues in the'
    ' vicinity of tau will be computed [default: as set in the examples]',
}

def main():
    parser = OptionParser(usage=usage, version='%prog ' + sfepy.__version__)
    parser.add_option('-c', '--conf', metavar='"key : value, ..."',
                      action='store', dest='conf', type='string',
                      default=None, help= help['conf'])
    parser.add_option('-O', '--options', metavar='"key : value, ..."',
                      action='store', dest='app_options', type='string',
                      default=None, help=help['options'])
    parser.add_option('-o', '', metavar='filename',
                      action='store', dest='output_filename_trunk',
                      default=None, help=help['filename'])
    parser.add_option('--oscillator',
                      action='store_true', dest='oscillator',
                      default=False, help=help['oscillator'])
    parser.add_option('--well',
                      action='store_true', dest='well',
                      default=False, help=help['well'])
    parser.add_option('--hydrogen',
                      action='store_true', dest='hydrogen',
                      default=False, help=help['hydrogen'])
    parser.add_option('--boron',
                      action='store_true', dest='boron',
                      default=False, help=help['boron'])
    parser.add_option('-n', '--n-eigs', type='int', metavar='int',
                      action='store', dest='n_eigs',
                      default=None, help=help['n_eigs'])
    parser.add_option('-t', '--tau', type='float', metavar='float',
                      action='store', dest='tau',
                      default=None, help=help['tau'])

    options, args = parser.parse_args()

    if len(args) == 1:
        filename_in = args[0];

    elif len(args) == 0:
        if options.oscillator:
            filename_in = fix_path("examples/quantum/oscillator.py")

        elif options.well:
            filename_in = fix_path("examples/quantum/well.py")

        elif options.hydrogen:
            filename_in = fix_path("examples/quantum/hydrogen.py")

        elif options.boron:
            filename_in = fix_path("examples/quantum/boron.py")

        else:
            parser.print_help()
            return

    else:
        parser.print_help()
        return

    define_args = {}

    if options.n_eigs is not None:
        define_args['n_eigs'] = options.n_eigs

    if options.tau is not None:
        define_args['tau'] = options.tau

    required, other = get_standard_keywords()
    conf = ProblemConf.from_file_and_options(filename_in, options,
                                             required, other,
                                             define_args=define_args)

    app = SchroedingerApp(conf, options, 'schroedinger:')
    opts = conf.options
    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
