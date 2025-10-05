#!/usr/bin/env python
"""
Generate available fields table for ReST documentation.
"""
import os.path as op
import sys

from argparse import ArgumentParser

sys.path.append('.')

import sfepy
from sfepy import get_paths
from sfepy.base.base import load_classes
from sfepy.discrete import Field

field_files = [ii for ii
               in get_paths('sfepy/discrete/fem/fields*.py')
               if 'fields_base.py' not in ii]
field_files += get_paths('sfepy/discrete/structural/fields*.py')
field_files += get_paths('sfepy/discrete/iga/fields*.py')
field_files += get_paths('sfepy/discrete/dg/fields.py')
field_table = load_classes(field_files, [Field], ignore_errors=True,
                           name_attr='family_name')


header = r"""
.. tabularcolumns:: |l|l|l|p{0.55\linewidth}|
.. list-table:: Fields
   :widths: 5 15 15 65
   :header-rows: 1

   * - space
     - basis
     - region kind
     - description
"""

table_row = """   * - %s
     - %s
     - %s
     - %s
"""

class_link = ':class:`{short} <{full}>`'

translate_rtype = {'volume' : 'cell', 'surface' : 'facet'}

def typeset_field_table(fd, field_table):

    fd.write('.. _field_table:\n')
    fd.write(header)

    rows = {}
    for key, cls in field_table.items():
        rtype, space, *basis = key.split('_')
        rtype = translate_rtype[rtype]
        basis = '_'.join(basis)
        doc = cls.__doc__
        if doc is not None:
            desc = doc.strip().splitlines()[0]

        else:
            desc = ''

        sdata = rows.setdefault(space, {})
        row = sdata.setdefault(basis, ([], desc))

        name = cls.__module__ + '.' + cls.__name__
        row[0].append(class_link.format(short=rtype, full=name))
        if row[1] is None:
            row[1] = desc

    for space, sdata in rows.items():
        for basis, (rtypes, desc) in sorted(sdata.items(), key=lambda x: x[0]):
            rtype = ', '.join(rtypes)
            fd.write(table_row % (space, basis, rtype, desc))

    fd.write('\n')


def typeset(fd):
    """
    Utility function called by Sphinx.
    """
    fd = open(fd, 'w')
    typeset_field_table(fd, field_table)
    fd.close()


def gen_field_table(app):
    typeset(op.join(app.builder.srcdir, 'field_table.rst'))


def setup(app):
    app.connect('builder-inited', gen_field_table)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s ' + sfepy.__version__)
    parser.add_argument('-o', '--output', metavar='output_filename',
                        action='store', dest='output_filename',
                        default='field_table.rst')
    options = parser.parse_args()

    typeset(options.output_filename)


if __name__ == '__main__':
    main()
