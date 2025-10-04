#!/usr/bin/env python
"""
Extract information from a SfePy multi-time-step results file (HDF5
format) and/or linearize results with stored higher order DOFs.

For the linearization, the original input (problem description) file must
be specified as the first argument. Use the option --linearization below
to override linearization parameters defined in the input file. The
linearization forces --dump option, i.e., output to VTK files.

Examples
--------

Extract variables according to an extraction list::

  python3 sfepy/scripts/extractor.py -e "p e 0 1999" bone.h5
  python3 sfepy/scripts/extractor.py -e "p e 0 1999" bone.h5 -a
  python3 sfepy/scripts/extractor.py -e "p e 0 1999" bone.h5 -o extracted.h5
  python3 sfepy/scripts/extractor.py -e "p e 0 1999" bone.h5 -o extracted.h5 -a
"""
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import sfepy
from sfepy.base.base import nm, dict_to_struct, get_default, Struct
from sfepy.base.ioutils import get_trunk
import sfepy.postprocess.time_history as th

def create_problem(filename):
    from sfepy.discrete import Problem

    problem = Problem.from_conf_file(filename,
                                     init_equations=False, init_solvers=False)
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

helps = {
    'debug':
    'automatically start debugger when an exception is raised',
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'dump' :
    'dump to sequence of VTK files',
    'same_dir' :
    'store the dumped VTK files in the directory of filename_in',
    'linearization' :
    'linearization options. Default values apply if neither command'
    ' line nor input file options are set.'
    " [default: 'kind:adaptive,min_level:0,max_level:2,eps:1e-2']",
    'times' :
    'extract and print times of individual time steps',
    'from' :
    'start dumping from time step ii [default: %(default)s]',
    'to' :
    'stop dumping at time step ii [default: <last step>]',
    'step' :
    'use every ii-th step for dumping [default: %(default)s]',
    'extract' :
    'extract variables according to extraction list.'
    " Example: 'u n 10 15, p e 0' means variable 'u' in nodes 10, 15"
    " and variable 'p' in element 0",
    'average' :
    'average vertex variable into cells ("e" extraction mode)'
}

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + sfepy.__version__)
    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        default=False, help=helps['debug'])
    parser.add_argument('-o', metavar='filename',
                        action='store', dest='output_filename_trunk',
                        default=None, help=helps['filename'])
    parser.add_argument('-d', '--dump', action='store_true', dest='dump',
                        default=False, help=helps['dump'])
    parser.add_argument('--same-dir', action='store_true', dest='same_dir',
                        default=False, help=helps['same_dir'])
    parser.add_argument('-l', '--linearization', metavar='options',
                        action='store', dest='linearization',
                        default=None, help=helps['linearization'])
    parser.add_argument('--times', action='store_true', dest='times',
                        default=False, help=helps['times'])
    parser.add_argument('-f', '--from', type=int, metavar='ii',
                        action='store', dest='step_from',
                        default=0, help=helps['from'])
    parser.add_argument('-t', '--to', type=int, metavar='ii',
                        action='store', dest='step_to',
                        default=None, help=helps['to'])
    parser.add_argument('-s', '--step', type=int, metavar='ii',
                        action='store', dest='step_by',
                        default=1, help=helps['step'])
    parser.add_argument('-e', '--extract', metavar='list',
                        action='store', dest='extract',
                        default=None, help=helps['extract'])
    parser.add_argument('-a', '--average', action='store_true',
                        dest='average', default=False, help=helps['average'])
    parser.add_argument('input_file', nargs='?', default=None)
    parser.add_argument('results_file')
    options = parser.parse_args()

    if options.debug:
        from sfepy.base.base import debug_on_error; debug_on_error()

    filename_in = options.input_file
    filename_results = options.results_file

    if filename_in is None:
        linearize = False
    else:
        linearize = True
        options.dump = True

    if options.times:
        steps, times, nts, dts = th.extract_times(filename_results)
        for ii, time in enumerate(times):
            step = steps[ii]
            print('%d %e %e %e' % (step, time, nts[ii], dts[ii]))

    if options.dump:
        trunk = get_default(options.output_filename_trunk,
                            get_trunk(filename_results))
        if options.same_dir:
            trunk = os.path.join(os.path.dirname(filename_results),
                                 os.path.basename(trunk))

        args = {}
        if linearize:
            problem = create_problem(filename_in)

            linearization = Struct(kind='adaptive', min_level=0,
                                   max_level=2, eps=1e-2)
            aux = problem.conf.options.get('linearization', None)
            linearization.update(aux)

            if options.linearization is not None:
                aux = parse_linearization(options.linearization)
                linearization.update(aux)

            args.update({'fields' : problem.fields,
                         'linearization' : linearization})

        if options.step_to is None:
            args.update({'step0' : options.step_from})

        else:
            args.update({'steps' : nm.arange(options.step_from,
                                             options.step_to + 1,
                                             options.step_by, dtype=nm.int32)})

        th.dump_to_vtk(filename_results, output_filename_trunk=trunk, **args)

    if options.extract:
        ths, ts = th.extract_time_history(filename_results, options.extract)

        if options.average:
            ths = th.average_vertex_var_in_cells(ths)

        if options.output_filename_trunk:
            th.save_time_history(ths, ts, options.output_filename_trunk + '.h5')

        else:
            print(dict_to_struct(ths, flag=(1, 1, 1)).str_all())

if __name__ == '__main__':
    main()
