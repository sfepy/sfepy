#!/usr/bin/env python
"""
Notes on writing new test files:
--------------------------------

A test file can contain anything, but usually it is similar to a regular input
file (defining a test problem), with a mandatory Test class. This class holds
all the test_* functions, as well as the from_conf(), which serves to
initialize the test (conf is in fact the test file itself, options are
command-line options).

All variables defined in a test file are collected in 'conf' variable passed to
a Test.__init__(). For example, 'input_name' in test_input_*.py files is
accessible as 'conf.input_name'. This is usefull if the test class is defined
outside the test file, as the classes in tests_basic.py are.

The test_* functions are collected automatically by run_tests.py, with one
exception: if a certain order of their evaluation is required, a class
attribute 'test' of the Test class with a list of the test function names
should be defined (example: test_meshio.py)."""

from __future__ import absolute_import
import sys
import os
import os.path as op
import gc
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import sfepy
if sfepy.top_dir not in sys.path: sys.path.append(sfepy.top_dir)

from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.base import output
from sfepy.base.timing import Timer

class OutputFilter(object):

    def __init__(self, allowed_lines):
        self.allowed_lines = allowed_lines
        self.msg_type1 = ['...', '!!!', '+++', '---']
        self.msg_type2 = ['<<<', '>>>']
        self.start()

    def start(self):
        self.stdout = sys.stdout
        sys.stdout = self

    def write(self, msg):
        if self.stdout is not None:
            msg_type = msg[:3]
            if msg_type in self.allowed_lines:
                if msg_type in self.msg_type1:
                    msg = ''.join((msg[:3], '  ', msg[3:]))
                elif msg_type in self.msg_type2:
                    msg = msg[4:]
                self.stdout.write(msg)
                self.stdout.write('\n')

    def flush(self):
        self.stdout.flush()

    def stop(self):
        sys.stdout = self.stdout
        self.stdout = None

def run_test(conf_name, options, ifile):
    from sfepy.base.ioutils import ensure_path

    ensure_path(op.join(options.out_dir, 'any'))

    if options.filter_none or options.raise_on_error:
        of = None
    elif options.filter_less:
        of = OutputFilter(['<<<', '>>>', '...', '!!!', '+++', '---'])
    elif options.filter_more:
        of = OutputFilter(['+++', '---'])
    else:
        of = OutputFilter(['<<<', '+++', '---'])

    print('<<< [%d] %s' % (ifile, conf_name))
    orig_prefix = output.get_output_prefix()
    output.set_output_prefix('[%d] %s' % (ifile, orig_prefix))

    _required, other = get_standard_keywords()
    required = ['Test']

    num = 1
    timer = Timer('tests')
    try:
        conf = ProblemConf.from_file(conf_name, required, _required + other)
        test = conf.funmod.Test.from_conf(conf, options)
        num = test.get_number()
        ok = True
        print('>>> test instance prepared (%d test(s))' % num)
    except KeyboardInterrupt:
        print('>>> interrupted')
        sys.exit(0)
    except:
        print('--- test instance creation failed')
        if options.raise_on_error:
            raise
        ok, n_fail, n_total = False, num, num

    if ok:
        try:
            timer.start()
            output.set_output_prefix(orig_prefix)
            ok, n_fail, n_total = test.run(debug=options.raise_on_error,
                                           ifile=ifile)
            output.set_output_prefix('[%d] %s' % (ifile, orig_prefix))
            timer.stop()
        except KeyboardInterrupt:
            print('>>> interrupted')
            sys.exit(0)
        except Exception as e:
            print('>>> %s' % e.__class__)
            if options.raise_on_error:
                raise
            ok, n_fail, n_total = False, num, num

    if ok:
        print('>>> all passed in %.2f s' % timer.total)
    else:
        print('!!! %s test failed' % n_fail)

    if of is not None:
        of.stop()

    output.set_output_prefix(orig_prefix)

    return n_fail, n_total, timer.total

def wrap_run_tests(options):
    def run_tests(stats, dir_name, filenames):
        filenames = [filename for filename in sorted(filenames)
                     if (len(filename) > 8) and
                     filename[:5] == 'test_' and filename[-3:] == '.py']

        print('<<< directory: %s, test files: %d' % (dir_name, len(filenames)))

        for filename in filenames:
            conf_name = op.join(dir_name, filename)

            n_fail, n_total, test_time = run_test(conf_name, options,
                                                  stats[0] + 1)

            stats[0] += 1
            stats[1] += n_fail
            stats[2] += n_total
            stats[3] += test_time

            gc.collect()

    return run_tests

def get_dir(default):
    if sfepy.in_source_tree:
        out = default
    else:
        out = op.normpath(op.join(sfepy.data_dir, default))
    return out

helps = {
    'dir' : 'directory with tests [default: %(default)s]',
    'out_dir' : 'directory for storing test results and temporary files'
    ' [default: %(default)s]',
    'raise_on_error' : 'raise silenced exceptions to see what was wrong',
    'debug':
    'automatically start debugger when an exception is raised',
    'filter-none' : 'do not filter any messages',
    'filter-less' : 'filter output (suppress all except test messages)',
    'filter-more' : 'filter output (suppress all except test result messages)',
    'print-doc' : 'print the docstring of this file (howto write new tests)',
}

def main():
    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s " + sfepy.__version__)
    parser.add_argument("--print-doc",
                        action="store_true", dest="print_doc",
                        default=False, help=helps['print-doc'])
    parser.add_argument("-d", "--dir", metavar='directory',
                        action="store", dest="test_dir",
                        default=get_dir('tests'),
                        help=helps['dir'])
    parser.add_argument("-o", "--output", metavar='directory',
                        action="store", dest="out_dir",
                        default=get_dir('output-tests'),
                        help=helps['out_dir'])
    parser.add_argument("--raise",
                        action="store_true", dest="raise_on_error",
                        default=False, help=helps['raise_on_error'])
    parser.add_argument("--debug",
                        action="store_true", dest="debug",
                        default=False, help=helps['debug'])
    parser.add_argument("--filter-none",
                        action="store_true", dest="filter_none",
                        default=False, help=helps['filter-none'])
    parser.add_argument("--filter-less",
                        action="store_true", dest="filter_less",
                        default=False, help=helps['filter-less'])
    parser.add_argument("--filter-more",
                        action="store_true", dest="filter_more",
                        default=False, help=helps['filter-more'])
    parser.add_argument("test_filename", nargs="*", default=[])

    options = parser.parse_args()

    if options.print_doc:
        print(__doc__)
        return

    if options.debug:
        options.raise_on_error = True
        from sfepy.base.base import debug_on_error; debug_on_error()

    run_tests = wrap_run_tests(options)
    stats = [0, 0, 0, 0.0]

    if len(options.test_filename) > 0:
        for test_filename in options.test_filename:
            dirname, filename = op.split(test_filename)
            run_tests(stats, dirname, [filename])
    else:
        for dirpath, dirnames, filenames in os.walk(options.test_dir):
            run_tests(stats, dirpath, filenames)

    print('%d test file(s) executed in %.2f s, %d failure(s) of %d test(s)'
          % (stats[0], stats[3], stats[1], stats[2]))

    return stats[1]

if __name__ == '__main__':
    sys.exit(main())
