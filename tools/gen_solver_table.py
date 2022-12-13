#!/usr/bin/env python
"""
Generate available solvers table for ReST documentation.
"""
from __future__ import absolute_import
import os.path as op
import sys

from argparse import ArgumentParser

sys.path.append('.')

import sfepy
from sfepy.base.base import load_classes
from sfepy.solvers import (NonlinearSolver, TimeSteppingSolver,
                           TimeStepController, LinearSolver, EigenvalueSolver,
                           QuadraticEVPSolver, OptimizationSolver)
from sfepy.solvers.auto_fallback import AutoFallbackSolver

solver_by_type_table = [
    [[AutoFallbackSolver], "Virtual Solvers with Automatic Fallback"],
    [[TimeSteppingSolver], "Time-Stepping Solvers"],
    [[TimeStepController], "Time Step Controllers"],
    [[NonlinearSolver], "Nonlinear Solvers"],
    [[LinearSolver], "Linear Solvers"],
    [[EigenvalueSolver], "Eigenvalue Problem Solvers"],
    [[QuadraticEVPSolver], "Quadratic Eigenvalue Problem Solvers"],
    [[OptimizationSolver], "Optimization Solvers"]
]

for i in enumerate(solver_by_type_table):
    solver_by_type_table[i[0]][0] = \
        load_classes(sfepy.solvers.solver_files,
                     solver_by_type_table[i[0]][0],
                     package_name='sfepy.solvers')


def trim(docstring):
    """Trim and split (doc)string."""
    if not docstring:
        return ''
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxsize
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxsize:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)
    # Return a splitted string:
    return trimmed


def typeset_solvers_table(fd, solver_table):
    """
    Generate solvers table ReST output.
    """
    rest_tag_start = '.. <%s>\n'
    rest_tag_end = '.. </%s>\n'

    for solver_type in solver_table:
        fd.write(rest_tag_start % solver_type[1])
        for name, cls in sorted(solver_type[0].items()):
            fd.write('- :class:`%s <%s.%s>`: ' %
                     (name, cls.__module__, cls.__name__))
            fd.write('%s\n' % trim(cls.__doc__)[0])
        fd.write(rest_tag_end % solver_type[1])
        fd.write('\n')


def typeset(fd):
    """
    Utility function called by Sphinx.
    """
    fd = open(fd, 'w')
    typeset_solvers_table(fd, solver_by_type_table)
    fd.close()


def gen_solver_table(app):
    typeset(op.join(app.builder.srcdir, 'solver_table.rst'))


def setup(app):
    app.connect('builder-inited', gen_solver_table)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s " + sfepy.__version__
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="output_filename",
        action="store",
        dest="output_filename",
        default="solver_table.rst",
    )
    options = parser.parse_args()

    typeset(options.output_filename)


if __name__ == '__main__':
    main()
