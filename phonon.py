#!/usr/bin/env python
# 12.01.2007, c
from __future__ import absolute_import
from argparse import ArgumentParser

import sfepy
from sfepy.base.base import output
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.homogenization.band_gaps_app import AcousticBandGapsApp
from sfepy.base.plotutils import plt

helps = {
    'debug':
    'automatically start debugger when an exception is raised',
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'detect_band_gaps' :
    'detect frequency band gaps',
    'analyze_dispersion' :
    'analyze dispersion properties (low frequency domain)',
    'plot' :
    'plot frequency band gaps, assumes -b',
    'phase_velocity' :
    'compute phase velocity (frequency-independet mass only)'
}

def main():
    parser = ArgumentParser()
    parser.add_argument("--version", action="version",
                        version="%(prog)s " + sfepy.__version__)
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        default=False, help=helps['debug'])
    parser.add_argument("-o", metavar='filename',
                        action="store", dest="output_filename_trunk",
                        default=None, help=helps['filename'])
    parser.add_argument("-b", "--band-gaps",
                        action="store_true", dest="detect_band_gaps",
                        default=False, help=helps['detect_band_gaps'])
    parser.add_argument("-d", "--dispersion",
                        action="store_true", dest="analyze_dispersion",
                        default=False, help=helps['analyze_dispersion'])
    parser.add_argument("-p", "--plot",
                        action="store_true", dest="plot",
                        default=False, help=helps['plot'])
    parser.add_argument("--phase-velocity",
                        action="store_true", dest="phase_velocity",
                        default=False, help=helps['phase_velocity'])
    parser.add_argument("filename_in")
    options = parser.parse_args()

    if options.debug:
        from sfepy.base.base import debug_on_error; debug_on_error()

    if options.plot:
        if plt is None:
            output('matplotlib.pyplot cannot be imported, ignoring option -p!')
            options.plot = False
        elif options.analyze_dispersion == False:
            options.detect_band_gaps = True

    required, other = get_standard_keywords()
    required.remove('equations')
    if not options.analyze_dispersion:
        required.remove('solver_[0-9]+|solvers')
    if options.phase_velocity:
        required = [ii for ii in required if 'ebc' not in ii]
    conf = ProblemConf.from_file(options.filename_in, required, other)

    app = AcousticBandGapsApp(conf, options, 'phonon:')
    opts = conf.options
    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
