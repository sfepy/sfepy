#!/usr/bin/env python
"""
Show terms use in problem description files in the given directory.
"""
from __future__ import absolute_import
import sys
import six
sys.path.append('.')
import os
from argparse import ArgumentParser

from sfepy.base.base import output, dict_from_keys_init, ordered_iteritems
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.ioutils import locate_files
from sfepy.discrete.equations import parse_definition
from sfepy.terms import term_table

helps = {
    'counts' : 'show terms use counts only',
    'unused' : 'show unused terms only',
}

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-c', '--counts',
                        action='store_true', dest='counts',
                        default=False, help=helps['counts'])
    parser.add_argument('-u', '--unused',
                        action='store_true', dest='unused',
                        default=False, help=helps['unused'])
    parser.add_argument('directory')
    options = parser.parse_args()

    pdf_dir = os.path.realpath(options.directory)

    required, other = get_standard_keywords()

    terms_use = dict_from_keys_init(term_table.keys(), set)

    for filename in locate_files('*.py', pdf_dir):
        base = filename.replace(pdf_dir, '').lstrip(os.path.sep)
        output('trying "%s"...' % base)

        try:
            conf = ProblemConf.from_file(filename, required, other,
                                         verbose=False)

        except:
            output('...failed')
            continue

        use = conf.options.get('use_equations', 'equations')
        eqs_conf = getattr(conf, use)
        for key, eq_conf in six.iteritems(eqs_conf):
            term_descs = parse_definition(eq_conf)
            for td in term_descs:
                terms_use[td.name].add(base)

        output('...ok')
    output('...done')

    if options.unused:
        output('unused terms:')

        unused = [name for name in terms_use.keys()
                  if len(terms_use[name]) == 0]
        for name in sorted(unused):
            output('  ' + name)

        output('total: %d' % len(unused))

    else:
        output('terms use:')
        for name, ex_names in ordered_iteritems(terms_use):
            output('%s: %d' % (name, len(ex_names)))
            if not options.counts:
                for ex_name in sorted(ex_names):
                    output('  ' + ex_name)

if __name__ == '__main__':
    main()
