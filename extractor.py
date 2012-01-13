#!/usr/bin/env python
"""
Examples
--------

$ ./extractor.py -e "p e 0 1999" bone.h5
$ ./extractor.py -e "p e 0 1999" bone.h5 -a
$ ./extractor.py -e "p e 0 1999" bone.h5 -o extracted.h5
$ ./extractor.py -e "p e 0 1999" bone.h5 -o extracted.h5 -a
"""
import os
from optparse import OptionParser

import sfepy
from sfepy.base.base import nm, dict_to_struct, get_default, Struct
from sfepy.base.ioutils import get_trunk
import sfepy.postprocess.time_history as th

def create_problem(filename):
    from sfepy.fem import ProblemDefinition

    problem = ProblemDefinition.from_conf_file(filename,
                                               init_equations=False,
                                               init_solvers=False)
    return problem

def parse_linearization(linearization):
    out = {}
    for item in linearization.split(','):
        key, val = item.split(':')
        if key == 'eps':
            val = float(val)
        elif key in ('min_level', 'max_level'):
            val = int(val)
        elif key == 'kind':
            pass
        else:
            raise ValueError('wrong linearization option key! (%s)'
                             % key)
        out[key] = val

    return dict_to_struct(out)

usage = """%prog [options] [<input file>] <results file>

Extract information from a SfePy multi-time-step results file (HDF5
format) and/or linearize results with stored higher order DOFs.

For the linearization, the original input (problem description) file must
be specified as the first argument. Use the option --linearization below
to override linearization parameters defined in the input file. The
linearization forces --dump option, i.e., output to VTK files.
"""

help = {
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'dump' :
    'dump to sequence of VTK files',
    'same_dir' :
    'store the dumped VTK files in the directory of filename_in',
    'linearization' :
    'linearization options.'
    " Example: 'kind:adaptive,min_level:0,max_level:2,eps:1e-2'"
    ' [default: %default]',
    'times' :
    'extract and print times of individual time steps',
    'from' :
    'start dumping from time step ii [default: %default]',
    'to' :
    'stop dumping at time step ii [default: <last step>]',
    'step' :
    'use every ii-th step for dumping [default: %default]',
    'extract' :
    'extract variables according to extraction list.'
    " Example: 'u n 10 15, p e 0' means variable 'u' in nodes 10, 15"
    " and variable 'p' in element 0",
    'average' :
    'average vertex variable into cells ("e" extraction mode)'
}

def main():
    parser = OptionParser(usage=usage, version="%prog " + sfepy.__version__)
    parser.add_option('-o', '', metavar='filename',
                      action='store', dest='output_filename_trunk',
                      default=None, help=help['filename'])
    parser.add_option('-d', '--dump', action='store_true', dest='dump',
                       default=False, help=help['dump'])
    parser.add_option('', '--same-dir', action='store_true', dest='same_dir',
                      default=False, help=help['same_dir'])
    parser.add_option('-l', '--linearization', metavar='options',
                      action='store', dest='linearization',
                      default=None, help=help['linearization'])
    parser.add_option('', '--times', action='store_true', dest='times',
                      default=False, help=help['times'])
    parser.add_option('-f', '--from', type=int, metavar='ii',
                      action='store', dest='step_from',
                      default=0, help=help['from'])
    parser.add_option('-t', '--to', type=int, metavar='ii',
                      action='store', dest='step_to',
                      default=None, help=help['to'])
    parser.add_option('-s', '--step', type=int, metavar='ii',
                      action='store', dest='step_by',
                      default=1, help=help['step'])
    parser.add_option('-e', '--extract', metavar='list',
                      action='store', dest='extract',
                      default=None, help=help['extract'])
    parser.add_option('-a', '--average', action='store_true', dest='average',
                      default=False, help=help['average'])

    (options, args) = parser.parse_args()

    nargs = len(args)
    if nargs == 1:
        filename_results = args[0]
        linearize = False

    elif nargs == 2:
        filename_in, filename_results = args
        linearize = True
        options.dump = True

    else:
        parser.print_help()
        return

    if options.times:
        times, nts, dts = th.extract_times(filename_results)
        for step, time in enumerate(times):
            print '%d %e %e %e' % (step, time, nts[step], dts[step])

    if options.dump:
        trunk = get_default(options.output_filename_trunk,
                            get_trunk(filename_results))
        if options.same_dir:
            trunk = os.path.join(os.path.dirname(filename_results),
                                 os.path.basename(trunk))

        args = {}
        if linearize:
            problem = create_problem(filename_in)

            if options.linearization is not None:
                linearization = parse_linearization(options.linearization)

            else:
                linearization = problem.conf.options.get('linearization', None)
                if linearization is not None:
                    linearization = dict_to_struct(linearization)

                else:
                    linearization = Struct(kind='adaptive', min_level=0,
                                           max_level=2, eps=1e-2)

            args.update({'fields' : problem.fields,
                         'linearization' : linearization})

        if options.step_to is None:
            args.update({'step0' : options.step_from})

        else:
            args.update({'steps' : nm.arange(options.step_from,
                                             options.step_to + 1,
                                             options.step_by, dtype=nm.int)})

        th.dump_to_vtk(filename_results, output_filename_trunk=trunk, **args)

    if options.extract:
        ths, ts = th.extract_time_history(filename_results, options.extract)

        if options.average:
            ths = th.average_vertex_var_in_cells(ths)

        if options.output_filename_trunk:
            th.save_time_history(ths, ts, options.output_filename_trunk + '.h5')

        else:
            print dict_to_struct(ths, flag=(1, 1, 1)).str_all()

if __name__ == '__main__':
    main()
