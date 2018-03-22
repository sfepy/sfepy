#!/usr/bin/env python

"""
Generate table of available solvers for the Sphinx documentation.
"""
from __future__ import absolute_import, generators
import os.path as op
import sys

from argparse import ArgumentParser

sys.path.append('.')

import sfepy
from sfepy.base.base import load_classes
from sfepy.solvers import NonlinearSolver, TimeSteppingSolver, LinearSolver, \
    EigenvalueSolver, OptimizationSolver

solver_by_type_table = [
    [[LinearSolver], "Linear Solvers"],
    [[NonlinearSolver], "Non-linear Solvers"],
    [[TimeSteppingSolver], "Time-Stepping Solvers"],
    [[EigenvalueSolver], "Eigen Value Solvers"],
    [[OptimizationSolver], "Optimization Solvers"]
]

for i in enumerate(solver_by_type_table):
    solver_by_type_table[i[0]][0] = \
        load_classes(sfepy.solvers.solver_files,
                     solver_by_type_table[i[0]][0],
                     package_name='sfepy.solvers')


def paragraphs(fileobj, separator='\n'):
    """
    Read iterable string by paragraphs.
    :param fileobj: iterable text object
    :param separator: paragraph separator
    :return: generator object
    """
    if separator[-1:] != '\n':
        separator += '\n'
    paragraph = []
    for line in fileobj:
        if line == separator:
            if paragraph:
                yield ''.join(paragraph)
                paragraph = []
        else:
            paragraph.append(line)
    if paragraph:
        yield ''.join(paragraph)


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
            fd.write('*%s*\n' % name)
            fd.write(next(paragraphs(cls.__doc__, '\n')))
            fd.write('\n')
        fd.write('\n')


def typeset(fd):
    """
    Utility function called by Sphinx.
    """

    fd = open(fd, 'w')
    typeset_solver_tables(fd, solver_by_type_table)
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
