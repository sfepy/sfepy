#!/usr/bin/env python
"""
Generate the table of all terms for the sphinx documentation.
"""
import os
import sys
from optparse import OptionParser
import pyparsing as pp

import numpy as nm

sys.path.append('.')
import sfepy.discrete.fem # Hack: fix circular dependency, as terms.pyx imports
                          # from sfepy.discrete.fem
from sfepy.terms import term_table

def set_section(sec):
    def action(str, loc, toks):
        if toks:
            sec[0] = toks[0][1:-1]
        return toks
    return action

def to_list(slist, sec):
    def action(str, loc, toks):
        if toks:
            slist.append((sec[0], toks[0]))
        return toks
    return action

def create_parser(slist, current_section):

    colon = pp.Literal(':')

    section = pp.Combine(colon
                         + pp.Word(pp.alphas, pp.alphanums + '_ ')
                         + colon)
    section.setParseAction(set_section(current_section))
    section.setName('section')

    text = pp.SkipTo(section | pp.StringEnd())
    text.setParseAction(to_list(slist, current_section))
    text.setName('text')

    doc = pp.StringStart()\
           + pp.Optional(text) + pp.ZeroOrMore(section + text)\
           + pp.StringEnd()

    return doc

header = """
.. tabularcolumns:: |p{0.3\linewidth}|p{0.2\linewidth}|p{0.5\linewidth}|
.. list-table:: Table of all terms.
   :widths: 30 30 40
   :header-rows: 1

   * - name/class/link
     - arguments
     - definition
"""

table_row = """   * - %s

         :class:`%s`
         :mod:`%s <%s>`
     - %s
     -
%s
"""

def format_next(text, new_text, pos, can_newline, width, ispaces):
    new_len = len(new_text)

    if (pos + new_len > width) and can_newline:
        text += '\n' + ispaces + new_text
        pos = new_len
        can_newline = False

    else:
        if pos > 0:
            text += ' ' + new_text
            pos += new_len + 1

        else:
            text += new_text
            pos += new_len

        can_newline = True

    return text, pos, can_newline

def typeset_to_indent(txt, indent0, indent, width):
    if not len(txt): return txt

    txt_lines = txt.strip().split('\n')

    ispaces = ' ' * indent
    text = (' ' * indent0) + txt_lines[0] + '\n' + ispaces

    can_newline = False
    pos = indent0
    for line in txt_lines[1:]:
        for word in line.split():
            text, pos, can_newline = format_next(text, word, pos, can_newline,
                                                 width, ispaces)

    return text

def typeset_term_syntax(term_class):
    if ((len(term_class.arg_types) > 1) and not
        isinstance(term_class.arg_types[0], str)):
        arg_types = [', '.join(['<%s>' % arg for arg in arg_type])
                     for arg_type in term_class.arg_types]
        arg_types = [' ``%s``' % arg_type for arg_type in arg_types]
        text = '\n\n       '.join(arg_types)

    else:
        text = ', '.join(['<%s>' % arg for arg in term_class.arg_types])
        text = '``%s``' % text
    return text

def typeset_term_table(fd, table):
    """Terms are sorted by name without the d*_ prefix."""
    sec_list = []
    current_section = ['']
    parser = create_parser(sec_list, current_section)

    fd.write(header)

    keys = table.keys()
    sort_keys = [key[key.find('_'):] for key in keys]
    iis = nm.argsort(sort_keys)
    for ii in iis:
        key = keys[ii]
        item_class = table[key]
        doc = item_class.__doc__

        if doc is not None:
            sec_list[:] = []
            current_section[0] = ''
            parser.parseString(doc)

            dd = [x[1] for x in sec_list if x[0].lower() == 'definition']
            if len(dd):
                dd = dd[0]
            else:
                dd = ''

            dds = dd.split('\n\n')
            definition = '\n\n'.join(typeset_to_indent(dd, 7, 11, 65)
                                     for dd in dds)

            fd.write(table_row % (item_class.name,
                                  item_class.__name__,
                                  item_class.__module__,
                                  item_class.__module__,
                                  typeset_term_syntax(item_class),
                                  definition))

    fd.write('\n')

def typeset(filename):
    """Utility function called by sphinx. """
    fd = open(filename, 'w')
    typeset_term_table(fd, term_table)
    fd.close()

def gen_term_table(app):
    typeset(os.path.join(app.builder.srcdir, 'term_table.rst'))

def setup(app):
    app.connect('builder-inited', gen_term_table)

usage = '%prog [options]\n' + __doc__.rstrip()

helps = {
    'output_filename' :
    'output file name',
}

def main():

    parser = OptionParser(usage=usage, version="%prog")
    parser.add_option("-o", "--output", metavar='output_filename',
                      action="store", dest="output_filename",
                      default="term_table.rst", help=helps['output_filename'])
    (options, args) = parser.parse_args()

    typeset(options.output_filename)

if __name__ == '__main__':
    main()
