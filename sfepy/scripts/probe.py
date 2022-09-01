#!/usr/bin/env python
# 12.01.2007, c
"""
Probe finite element solutions in points defined by various geometrical probes.

In the examples below it is supposed that sfepy is installed. When using the
in-place build, replace ``sfepy-probe`` by ``python3 sfepy/scripts/probe.py``.

Generation mode
---------------
sfepy-probe [generation options] <input file> <results file>

Probe the data in the results file corresponding to the problem defined in the
input file. The input file options must contain 'gen_probes' and 'probe_hook'
keys, pointing to proper functions accessible from the input file scope.

For each probe returned by `gen_probes()` a data plot figure and a text
file with the data plotted are saved, see the options below.

Generation options
------------------
-o, --auto-dir, --same-dir, -f, --only-names, -s

Postprocessing mode
-------------------
sfepy-probe [postprocessing options] <probe file> <figure file>

Read a previously probed data from the probe text file, re-plot them,
and integrate them along the probe.

Postprocessing options
----------------------
--postprocess, --radial, --only-names

Notes
-----
For extremely thin hexahedral elements the Newton's iteration for finding the
reference element coordinates might converge to a spurious solution outside
of the element. To obtain some values even in this case, try increasing the
--close-limit option value.
"""
from __future__ import absolute_import
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import numpy as nm

import sfepy
from sfepy.base.base import output, assert_
from sfepy.base.ioutils import edit_filename
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.discrete import Problem
from sfepy.discrete.fem import MeshIO
from sfepy.discrete.probes import write_results, read_results
import six

helps = {
    'debug':
    'automatically start debugger when an exception is raised',
    'filename' :
    'basename of output file(s) [default: <basename of input file>]',
    'output_format' :
    'output figure file format (supported by the matplotlib backend used) '\
    '[default: %(default)s]',
    'auto_dir' :
    'the directory of the results file is determined automatically using the '\
    '"output_dir" option in input file options',
    'same_dir' :
    'store the probe figures/data in the directory of the results file',
    'only_names' :
    'probe only named data',
    'step' :
    'probe the given time step',
    'close_limit' :
    'maximum limit distance of a point from the closest element allowed'
    ' for extrapolation. [default: %(default)s]',
    'postprocess' :
    'postprocessing mode',
    'radial' :
    'assume radial integration',
}

def generate_probes(filename_input, filename_results, options,
                    conf=None, problem=None, probes=None, labels=None,
                    probe_hooks=None):
    """
    Generate probe figures and data files.
    """
    if conf is None:
        required, other = get_standard_keywords()
        conf = ProblemConf.from_file(filename_input, required, other)

    opts = conf.options

    if options.auto_dir:
        output_dir = opts.get_('output_dir', '.')
        filename_results = os.path.join(output_dir, filename_results)

    output('results in: %s' % filename_results)

    io = MeshIO.any_from_filename(filename_results)
    step = options.step if options.step >= 0 else io.read_last_step()
    all_data = io.read_data(step)
    output('loaded:', list(all_data.keys()))
    output('from step:', step)

    if options.only_names is None:
        data = all_data
    else:
        data = {}
        for key, val in six.iteritems(all_data):
            if key in options.only_names:
                data[key] = val

    if problem is None:
        problem = Problem.from_conf(conf,
                                    init_equations=False, init_solvers=False)

    if probes is None:
        gen_probes = conf.get_function(conf.options.gen_probes)
        probes, labels = gen_probes(problem)

    if probe_hooks is None:
        probe_hooks = {None : conf.get_function(conf.options.probe_hook)}

    if options.output_filename_trunk is None:
            options.output_filename_trunk = problem.ofn_trunk

    filename_template = options.output_filename_trunk \
                        + ('_%%d.%s' % options.output_format)
    if options.same_dir:
        filename_template = os.path.join(os.path.dirname(filename_results),
                                         filename_template)

    output_dir = os.path.dirname(filename_results)

    for ip, probe in enumerate(probes):
        output(ip, probe.name)

        probe.set_options(close_limit=options.close_limit)

        for key, probe_hook in six.iteritems(probe_hooks):

            out = probe_hook(data, probe, labels[ip], problem)
            if out is None: continue
            if isinstance(out, tuple):
                fig, results = out
            else:
                fig = out

            if key is not None:
                filename = filename_template % (key, ip)

            else:
                filename = filename_template % ip

            if fig is not None:
                if isinstance(fig, dict):
                    for fig_name, fig_fig in six.iteritems(fig):
                        fig_filename = edit_filename(filename,
                                                     suffix='_' + fig_name)
                        fig_fig.savefig(fig_filename)
                        output('figure ->', os.path.normpath(fig_filename))

                else:
                    fig.savefig(filename)
                    output('figure ->', os.path.normpath(filename))

            if results is not None:
                txt_filename = edit_filename(filename, new_ext='.txt')

                write_results(txt_filename, probe, results)

                output('data ->', os.path.normpath(txt_filename))

def integrate_along_line(x, y, is_radial=False):
    """
    Integrate numerically (trapezoidal rule) a function :math:`y=y(x)`.

    If is_radial is True, multiply each :math:`y` by :math:`4 \pi x^2`.
    """
    dx = nm.diff(x)
    ay = 0.5 * (y[:-1] + y[1:])

    if is_radial:
        ax = 0.5 * (x[:-1] + x[1:])
        val = 4.0 * nm.pi * nm.sum(ay * dx * (ax**2))

    else:
        val = nm.sum(ay * dx)

    return val

def postprocess(filename_input, filename_results, options):
    """
    Postprocess probe data files - replot, integrate data.
    """
    from matplotlib import pyplot as plt

    header, results = read_results(filename_input,
                                   only_names=options.only_names)
    output(header)

    fig = plt.figure()
    for name, result in six.iteritems(results):
        pars, vals = result[:, 0], result[:, 1]

        ii = nm.where(nm.isfinite(vals))[0]
        # Nans only at the edges.
        assert_(nm.diff(ii).sum() == (len(ii)-1))

        val = integrate_along_line(pars[ii], vals[ii], options.radial)

        label = r'%s: $\int\ %s' % (name, name)
        if options.radial:
            label += ' (r)'
        label += '$ = %.5e'% val

        plt.plot(pars, vals, label=label, lw=0.2, marker='+', ms=1)
        plt.ylabel('probed data')
        plt.xlabel('probe coordinate')

        output(label)

    plt.legend()

    fig.savefig(filename_results)

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
    parser.add_argument('--auto-dir',
                        action='store_true', dest='auto_dir',
                        default=False, help=helps['auto_dir'])
    parser.add_argument('--same-dir',
                        action='store_true', dest='same_dir',
                        default=False, help=helps['same_dir'])
    parser.add_argument('-f', '--format', metavar='format',
                        action='store', dest='output_format',
                        default='png', help=helps['output_format'])
    parser.add_argument('--only-names', metavar='list of names',
                        action='store', dest='only_names',
                        default=None, help=helps['only_names'])
    parser.add_argument('-s', '--step', type=int, metavar='step',
                        action='store', dest='step',
                        default=0, help=helps['step'])
    parser.add_argument('-c', '--close-limit', type=float, metavar='distance',
                        action='store', dest='close_limit',
                        default=0.1, help=helps['close_limit'])
    parser.add_argument('-p', '--postprocess',
                        action='store_true', dest='postprocess',
                        default=False, help=helps['postprocess'])
    parser.add_argument('--radial',
                        action='store_true', dest='radial',
                        default=False, help=helps['radial'])
    parser.add_argument('filename_in')
    parser.add_argument('filename_out')
    options = parser.parse_args()

    if options.debug:
        from sfepy.base.base import debug_on_error; debug_on_error()

    filename_input = options.filename_in
    filename_results = options.filename_out

    if options.only_names is not None:
        options.only_names = options.only_names.split(',')

    output.prefix = 'probe:'

    if options.postprocess:
        postprocess(filename_input, filename_results, options)
    else:
        generate_probes(filename_input, filename_results, options)

if __name__ == '__main__':
    main()
