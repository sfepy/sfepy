#!/usr/bin/env python
"""
Generate the table of all terms for the sphinx documentation.
"""
import os
from sfepy.base.base import dict_from_keys_init
from sfepy.discrete.equations import parse_definition
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.ioutils import locate_files
from sfepy import get_paths
import sys
from argparse import ArgumentParser
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

newpage = r"""
.. raw:: latex

   \newpage
"""

header = r"""
.. tabularcolumns:: |p{0.15\linewidth}|p{0.10\linewidth}|p{0.6\linewidth}|p{0.15\linewidth}|
.. list-table:: %s terms
   :widths: 15 10 60 15
   :header-rows: 1
   :class: longtable

   * - name/class
     - arguments
     - definition
     - examples
"""

table_row = """   * - %s

       :class:`%s <%s.%s>`
     - %s
     - %s
     - %s
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
    if ((len(term_class.arg_types) >= 1) and not
            isinstance(term_class.arg_types[0], str)):
        is_param = len([arg for arg in term_class.arg_types[-1]
                        if arg.startswith('parameter')]) > 0
        is_vs = len([arg for arg in term_class.arg_types[0]
                     if arg in ['virtual', 'state']]) > 0

        arg_types_ = [list(k) for k in term_class.arg_types]

        if is_param and is_vs:
            at0 = arg_types_[0]
            at1 = arg_types_[-1]
            for k in range(len(at0)):
                if (at0[k] in ['virtual', 'state'] and
                        at1[k].startswith('parameter')):
                    aux = at1[k].replace('parameter', 'param')
                    at0[k] = at0[k] + '/' + aux

            arg_types_ = arg_types_[:-1]

        arg_types = [', '.join(['``<%s>``' % arg for arg in arg_type])
                     for arg_type in arg_types_]
        text = '\n\n       '.join(arg_types)
    else:
        text = ', '.join(['``<%s>``' % arg for arg in term_class.arg_types])
    return text

link_example = ':ref:`%s <%s>`'

omits = [
    'vibro_acoustic3d_mid.py',
    'its2D_5.py',
    'linear_elastic_probes.py',
    '__init__.py',
]

def typeset_examples(term_class, term_use):
    # e.g. fem-time_advection_diffusion -> tim.adv.dif.
    to_shorter_name = lambda st: '.'.join(
                    [s[:3] for s in st.split('-')[-1].split('_')])
    link_list = [(link_example % (to_shorter_name(exmpl), exmpl))
                    for exmpl in term_use[term_class.name]]
    return ', '.join(link_list)

def get_examples(table):

    term_use = dict_from_keys_init(table.keys(), set)
    required, other = get_standard_keywords()
    for filename in locate_files('*py', get_paths('sfepy/examples/')[0]):
        try:
            conf = ProblemConf.from_file(filename, required, other,
                                         verbose=False)
        except:
            continue

        ebase = filename.split('examples/')[1]
        lbase = os.path.splitext(ebase)[0]
        label = lbase.replace('/', '-')

        pyfile_name = ebase.split('/')[1]
        if pyfile_name in omits:
            continue

        use = conf.options.get('use_equations', 'equations')
        eqs_conf = getattr(conf, use)
        for key, eq_conf in eqs_conf.items():
            term_descs = parse_definition(eq_conf)
            for td in term_descs:
                term_use[td.name].add(label)

    return term_use

def typeset_term_table(fd, keys, table, title):
    """Terms are sorted by name without the d*_ prefix."""
    sec_list = []
    current_section = ['']
    parser = create_parser(sec_list, current_section)

    fd.write('.. _term_table_%s:\n' % title)
    label = 'Table of %s terms' % title
    fd.write(''.join([newpage, label, '\n', '"' * len(label), '\n']))
    fd.write(header % (title[0].upper() + title[1:]))

    term_use = get_examples(table)

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
            dds = dd.strip().split('\n\n')
            definition = '\n\n'.join(typeset_to_indent(dd, 7, 11, 65)
                                     for dd in dds)[7:]
            fd.write(table_row % (item_class.name,
                                  item_class.__name__,
                                  item_class.__module__,
                                  item_class.__name__,
                                  typeset_term_syntax(item_class),
                                  definition,
                                  typeset_examples(item_class, term_use)))

    fd.write('\n')

def typeset_term_tables(fd, table):
    """Generate tables: basic, sensitivity, special."""
    scattab = [
        ('_st_', 2),
        ('_sd_', 0),
        ('_adj_', 0),
        ('_tl_', 1),
        ('_ul_', 1),
        ('_th', 2),
        ('_eth', 2),
        ('_of_', 2),
        ('de_', 3)]

    new_tabs = [[],[],[],[]]
    for term_name in table.keys():
        for term_tag, tab_id in scattab:
            if term_tag in term_name:
                new_tabs[tab_id].append(term_name)
                break

    basic_keys = list(set(table.keys())
                      - set(new_tabs[0]) - set(new_tabs[1])
                      - set(new_tabs[2]) - set(new_tabs[3]))
    typeset_term_table(fd, basic_keys, table, 'basic')
    typeset_term_table(fd, new_tabs[0], table, 'sensitivity')
    typeset_term_table(fd, new_tabs[1], table, 'large deformation')
    typeset_term_table(fd, new_tabs[2], table, 'special')
    typeset_term_table(fd, new_tabs[3], table, 'multi-linear')
    fd.write(newpage)

def typeset(filename):
    """Utility function called by sphinx. """
    fd = open(filename, 'w')
    typeset_term_tables(fd, term_table)
    fd.close()

def gen_term_table(app):
    typeset(os.path.join(app.builder.srcdir, 'term_table.rst'))

def setup(app):
    app.connect('builder-inited', gen_term_table)

helps = {
    'output_filename' :
    'output file name',
}

def main():

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="%(prog)s")
    parser.add_argument("-o", "--output", metavar='output_filename',
                        action="store", dest="output_filename",
                        default="term_table.rst",
                        help=helps['output_filename'])
    options = parser.parse_args()

    typeset(options.output_filename)

if __name__ == '__main__':
    main()
