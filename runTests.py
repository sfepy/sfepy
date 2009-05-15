#!/usr/bin/env python
"""
Notes on writing new test files:
--------------------------------

A test file can contain anything, but usually it is similar to a regular input
file (defining a test problem), with a mandatory Test class. This class holds
all the test_* functions, as well as the from_conf(), which serves to initialize
the test (conf is in fact the test file itself, options are command-line
options).

All variables defined in a test file are collected in 'conf' variable passed to
a Test.__init__(). For example, 'input_name' in test_input_*.py files is
accessible as 'conf.input_name'. This is usefull if the test class is defined
outside the test file, as the classes in tests_basic.py are.

The test_* functions are collected automatically by run_tests.py, with one
exception: if a certain order of their evaluation is required, a class
attribute 'test' of the Test class with a list of the test function names
should be defined (example: test_meshio.py)."""

import sys
import time
import os
import os.path as op
from optparse import OptionParser

import sfepy
from sfepy.base.conf import ProblemConf, get_standard_keywords

##
# 05.06.2007, c
class OutputFilter( object ):

    ##
    # c: 05.06.2007, r: 05.02.2008
    def __init__( self, allowed_lines ):
        self.allowed_lines = allowed_lines
        self.msg_type1 = ['...', '!!!', '+++', '---']
        self.msg_type2 = ['<<<', '>>>']
        self.start()

    ##
    # 05.06.2007, c
    def start( self ):
        self.stdout = sys.stdout
        sys.stdout = self

    ##
    # c: 05.06.2007, r: 05.02.2008
    def write( self, msg ):
        if self.stdout is not None:
            msg_type = msg[:3]
            if msg_type in self.allowed_lines:
                if msg_type in self.msg_type1:
                    msg = ''.join( (msg[:3], '  ', msg[3:]) )
                elif msg_type in self.msg_type2:
                    msg = msg[4:]
                self.stdout.write( msg )
                self.stdout.write( '\n' )

    ##
    # 05.06.2007, c
    def stop( self ):
        sys.stdout = self.stdout
        self.stdout = None

##
# c: 04.06.2007, r: 19.02.2008
def run_test( conf_name, options ):
    try:
        os.makedirs( options.out_dir )
    except OSError, e:
        if e.errno != 17: # [Errno 17] File exists
            raise

    if options.filter_none or options.debug:
        of = None
    elif options.filter_less:
        of = OutputFilter( ['<<<', '>>>', '...', '!!!', '+++', '---'] )
    elif options.filter_more:
        of = OutputFilter( ['+++', '---'] )
    else:
        of = OutputFilter( ['<<<', '+++', '---'] )
        
    print '<<< %s' % conf_name

    _required, other = get_standard_keywords()
    required = ['Test']

    num = 1
    test_time = 0.0
    try:
        conf = ProblemConf.from_file( conf_name, required, _required + other )
        test = conf.funmod.Test.from_conf( conf, options )
        num = test.get_number()
        ok = True
        print '>>> test instance prepared (%d test(s))' % num
    except KeyboardInterrupt:
        print '>>> interrupted'
        sys.exit( 0 )
    except:
        print '--- test instance creation failed'
        if options.debug:
            raise
        ok, n_fail, n_total = False, num, num
        
    if ok:
        try:
            tt = time.clock()
            ok, n_fail, n_total = test.run( options.debug )
            test_time = time.clock() - tt
        except KeyboardInterrupt:
            print '>>> interrupted'
            sys.exit( 0 )
        except Exception, e:
            print '>>> %s' % e.__class__
            if options.debug:
                raise
            ok, n_fail, n_total = False, num, num

    if ok:
        print '>>> all passed in %.2f s' % test_time
    else:
        print '!!! %s test failed' % n_fail

    if of is not None:
        of.stop()
        
    return n_fail, n_total, test_time

##
# 30.05.2007, c
# 01.06.2007
# 04.06.2007
# 19.07.2007
def wrap_run_tests( options ):
    def run_tests( stats, dir_name, filenames ):
        filenames = [filename for filename in filenames
                     if (len( filename ) > 8) and
                     filename[:5] == 'test_' and filename[-3:] == '.py']

        print '<<< directory: %s, test files: %d' % (dir_name, len( filenames ))

        for filename in filenames:
            conf_name = op.join( dir_name, filename )

            n_fail, n_total, test_time = run_test( conf_name, options )

            stats[0] += 1
            stats[1] += n_fail
            stats[2] += n_total
            stats[3] += test_time

    return run_tests

def get_dir(default):
    if sfepy.in_source_tree:
        out = default
    else:
        share_path = '../../../../share/sfepy'
        out = op.normpath(op.join(sfepy.top_dir, share_path, default))
    return out

usage = """%prog [options] [test_filename]"""

help = {
    'dir' : 'directory with tests [default: <top_dir>/tests]',
    'out_dir' : 'directory for storing test results and temporary files'
    ' [default: %default]',
    'debug' : 'raise silenced exceptions to see what was wrong',
    'filter-none' : 'do not filter any messages',
    'filter-less' : 'filter output (suppress all except test messages)',
    'filter-more' : 'filter output (suppress all except test result messages)',
    'print-doc' : 'print the docstring of this file (howto write new tests)',
}

##
# c: 30.05.2007, r: 06.02.2008
def main():
    parser = OptionParser(usage = usage, version = "%prog " + sfepy.__version__)
    parser.add_option( "", "--print-doc",
                       action = "store_true", dest = "print_doc",
                       default = False, help = help['print-doc'] )
    parser.add_option( "-d", "--dir", metavar = 'directory',
                       action = "store", dest = "test_dir",
                       default = get_dir('tests'),
                       help = help['dir'] )
    parser.add_option( "-o", "--output", metavar = 'directory',
                       action = "store", dest = "out_dir",
                       default = get_dir('output-tests'),
                       help = help['out_dir'] )
    parser.add_option( "", "--debug",
                       action = "store_true", dest = "debug",
                       default = False, help = help['debug'] )
    parser.add_option( "", "--filter-none",
                       action = "store_true", dest = "filter_none",
                       default = False, help = help['filter-none'] )
    parser.add_option( "", "--filter-less",
                       action = "store_true", dest = "filter_less",
                       default = False, help = help['filter-less'] )
    parser.add_option( "", "--filter-more",
                       action = "store_true", dest = "filter_more",
                       default = False, help = help['filter-more'] )
    options, args = parser.parse_args()

    if options.print_doc:
        print __doc__
        return

    if len( args ) > 1:
        parser.print_help(),
        return
    elif len( args ) == 1:
        test_filename = args[0]

        stats = run_test( test_filename, options )
        print '1 test file executed in %.2f s, %d failure(s) of %d test(s)'\
              % (stats[2], stats[0], stats[1])
        return


    stats = [0, 0, 0, 0.0]
    op.walk( options.test_dir, wrap_run_tests( options ), stats )
    print '%d test file(s) executed in %.2f s, %d failure(s) of %d test(s)'\
          % (stats[0], stats[3], stats[1], stats[2])

if __name__ == '__main__':
    main()
