#!/usr/bin/env python
# 12.01.2007, c 
import os.path as op
from optparse import OptionParser

import sfepy
from sfepy.base.base import *
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.fem import MeshIO, ProblemDefinition

usage = """%prog [options] <input file> <results file>

Probe the data in the results file corresponding to the problem defined in the
input file. The input file options must contain 'gen_probes' and 'probe_hook'
keys, pointing to proper functions accessible from the input file scope.

For each probe returned by gen_probes() a data plot figure is saved, see the
options below."""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'output_format' :
    'output file format (supported by the matplotlib backend used) '\
    '[default: %default]',
    'auto_dir' :
    'the directory of the results file is determined automatically using the '\
    '"output_dir" option in input file options',
    'same_dir' :
    'store the probe figures in the directory of the results file',
    'only_names' :
    'probe only named data',
}

def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option( "-o", "", metavar = 'filename',
                       action = "store", dest = "output_filename_trunk",
                       default = None, help = help['filename'] )
    parser.add_option( "", "--auto-dir",
                       action = "store_true", dest = "auto_dir",
                       default = False, help = help['auto_dir'] )
    parser.add_option( "", "--same-dir",
                       action = "store_true", dest = "same_dir",
                       default = False, help = help['same_dir'] )
    parser.add_option("-f", "--format", metavar='format',
                      action="store", dest="output_format",
                      default="png", help=help['output_format'])
    parser.add_option("--only-names", metavar='list of names',
                      action="store", dest="only_names",
                      default=None, help=help['only_names'])

    options, args = parser.parse_args()
#    print options; pause()

    if (len( args ) == 2):
        filename_input, filename_results = args
    else:
        parser.print_help(),
        return

    if options.only_names is not None:
        options.only_names = options.only_names.split(',')
    
    output.prefix = 'probe:'

    required, other = get_standard_keywords()
    conf = ProblemConf.from_file( filename_input, required, other )
    opts = conf.options

    if options.auto_dir:
        output_dir = get_default_attr( opts, 'output_dir', '.' )
        filename_results = os.path.join(output_dir, filename_results)

    output('results in: %s' % filename_results)

    io = MeshIO.any_from_filename(filename_results)
    all_data = io.read_data(0)
    output('loaded:', all_data.keys())

    if options.only_names is None:
        data = all_data
    else:
        data = {}
        for key, val in all_data.iteritems():
            if key in options.only_names:
                data[key] = val

    problem = ProblemDefinition.from_conf(conf, options,
                                          init_equations=False,
                                          init_solvers=False)

    gen_probes = getattr(conf.funmod, conf.options.gen_probes)
    probe_hook = getattr(conf.funmod, conf.options.probe_hook)
    
    probes, labels = gen_probes(problem)

    if options.output_filename_trunk is None:
	    options.output_filename_trunk = problem.ofn_trunk

    filename_template = options.output_filename_trunk \
                        + ('_%%d.%s' % options.output_format)
    if options.same_dir:
        filename_template = os.path.join(os.path.dirname(filename_results),
                                         filename_template)

    for ip, probe in enumerate(probes):
        output(ip, probe.name)
        fig = probe_hook(data, probe, labels[ip], problem)

        if fig is not None:
            filename = filename_template % ip
            fig.savefig(filename)
            output('->', os.path.normpath(filename))

if __name__ == '__main__':
    main()
