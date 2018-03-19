#!/usr/bin/env python

"""
Generate table of available solvers for the Sphinx documentation.
"""
from __future__ import absolute_import
import os.path as op
import sys

from argparse import ArgumentParser

sys.path.append('.')

import sfepy
from sfepy.base.base import load_classes
from sfepy.solvers import NonlinearSolver, TimeSteppingSolver, LinearSolver, \
    EigenvalueSolver, OptimizationSolver

solver_files = sfepy.get_paths('sfepy/solvers/*.py')
remove = ['setup.py', 'solvers.py']
solver_files = [name for name in solver_files
                if op.basename(name) not in remove]

solvers_by_type_table = [
    [[TimeSteppingSolver], "Time-Stepping Solvers"],
    [[LinearSolver], "Linear Solvers"],
    [[NonlinearSolver], "Non-linear Solvers"],
    [[EigenvalueSolver], "Eigen Value Solvers"],
    [[OptimizationSolver], "Optimization Solvers"]
]

for i in enumerate(solvers_by_type_table):
    solvers_by_type_table[i[0]][0] = \
        load_classes(solver_files,
                     solvers_by_type_table[i[0]][0],
                     package_name='sfepy.solvers')


def typeset_solver_tables(fd, solver_tables):
    """
    Generate solver tables ReST output.
    """

    doc_sec_level = '"'

    for solver_type in solver_tables:
        doc_sec_label = "%s" % solver_type[1]
        fd.write(''.join([doc_sec_label,
                          '\n',
                          doc_sec_level * len(doc_sec_label),
                          '\n'])
                 )
        for name, cls in solver_type[0].items():
            fd.write('*%s*' % name)
            fd.write(cls.__doc__)
            fd.write('\n')
        fd.write('\n')


def typeset(fd):
    """
    Utility function called by Sphinx.
    """

    fd = open(fd, 'w')
    typeset_solver_tables(fd, solvers_by_type_table)
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
