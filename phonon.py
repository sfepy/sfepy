#!/usr/bin/env python
# 12.01.2007, c
from optparse import OptionParser

import sfepy
from sfepy.base.base import output
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.homogenization.band_gaps_app import AcousticBandGapsApp
from sfepy.base.plotutils import plt

usage = """%prog [options] filename_in"""

help = {
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
    parser = OptionParser(usage=usage, version="%prog " + sfepy.__version__)
    parser.add_option("-o", "", metavar='filename',
                      action="store", dest="output_filename_trunk",
                      default=None, help=help['filename'])
    parser.add_option("-b", "--band-gaps",
                      action="store_true", dest="detect_band_gaps",
                      default=False, help=help['detect_band_gaps'])
    parser.add_option("-d", "--dispersion",
                      action="store_true", dest="analyze_dispersion",
                      default=False, help=help['analyze_dispersion'])
    parser.add_option("-p", "--plot",
                      action="store_true", dest="plot",
                      default=False, help=help['plot'])
    parser.add_option("--phase-velocity",
                      action="store_true", dest="phase_velocity",
                      default=False, help=help['phase_velocity'])

    options, args = parser.parse_args()
    if options.plot:
        if plt is None:
            output('matplotlib.pyplot cannot be imported, ignoring option -p!')
            options.plot = False
        elif options.analyze_dispersion == False:
            options.detect_band_gaps = True

    if (len(args) == 1):
        filename_in = args[0];
    else:
        parser.print_help(),
        return

    required, other = get_standard_keywords()
    required.remove('equations')
    if not options.analyze_dispersion:
        required.remove('solver_[0-9]+|solvers')
    if options.phase_velocity:
        required.remove('ebc_[0-9]+|ebcs')
    conf = ProblemConf.from_file(filename_in, required, other)

    app = AcousticBandGapsApp(conf, options, 'phonon:')
    opts = conf.options
    if hasattr(opts, 'parametric_hook'): # Parametric study.
        parametric_hook = conf.get_function(opts.parametric_hook)
        app.parametrize(parametric_hook)
    app()

if __name__ == '__main__':
    main()
