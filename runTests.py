#!/usr/bin/env python
"""
Notes on writing new test files:
--------------------------------

A test file can contain anything, but usually it is similar to a regular input
file (defining a test problem), with a mandatory Test class. This class holds
all the test_* functions, as well as the fromConf(), which serves to initialize
the test (conf is in fact the test file itself, options are command-line
options).

All variables defined in a test file are collected in 'conf' variable passed to
a Test.__init__(). For example, 'inputName' in test_input_*.py files is
accessible as 'conf.inputName'. This is usefull if the test class is defined
outside the test file, as the classes in testsBasic.py are.

The test_* functions are collected automatically by runTests.py, with one
exception: if a certain order of their evaluation is required, a class
attribute 'test' of the Test class with a list of the test function names
should be defined (example: test_meshio.py)."""

import sys
import time
import os
import os.path as op
from optparse import OptionParser

import init_sfe
from sfe.base.conf import ProblemConf
from sfe.solvers.generic import getStandardKeywords
##
# 05.06.2007, c
class OutputFilter( object ):

    ##
    # c: 05.06.2007, r: 05.02.2008
    def __init__( self, allowedLines ):
        self.allowedLines = allowedLines
        self.msgType1 = ['...', '!!!', '+++', '---']
        self.msgType2 = ['<<<', '>>>']
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
            msgType = msg[:3]
            if msgType in self.allowedLines:
                if msgType in self.msgType1:
                    msg = ''.join( (msg[:3], '  ', msg[3:]) )
                elif msgType in self.msgType2:
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
def runTest( confName, options ):
    try:
        os.makedirs( options.outDir )
    except OSError, e:
        if e.errno != 17: # [Errno 17] File exists
            raise

    if options.filterNone or options.debug:
        of = None
    elif options.filterLess:
        of = OutputFilter( ['<<<', '>>>', '...', '!!!', '+++', '---'] )
    elif options.filterMore:
        of = OutputFilter( ['+++', '---'] )
    else:
        of = OutputFilter( ['<<<', '+++', '---'] )
        
    print '<<< %s' % confName

    _required, other = getStandardKeywords()
    required = ['Test']

    num = 1
    testTime = 0.0
    try:
        conf = ProblemConf.fromFile( confName, required, _required + other )
        test = conf.funmod.Test.fromConf( conf, options )
        num = test.getNumber()
        ok = True
        print '>>> test instance prepared (%d test(s))' % num
    except KeyboardInterrupt:
        print '>>> interrupted'
        sys.exit( 0 )
    except:
        print '--- test instance creation failed'
        if options.debug:
            raise
        ok, nFail, nTotal = False, num, num
        
    if ok:
        try:
            tt = time.clock()
            ok, nFail, nTotal = test.run( options.debug )
            testTime = time.clock() - tt
        except KeyboardInterrupt:
            print '>>> interrupted'
            sys.exit( 0 )
        except Exception, e:
            print '>>> %s' % e.__class__
            if options.debug:
                raise
            ok, nFail, nTotal = False, num, num

    if ok:
        print '>>> all passed in %.2f s' % testTime
    else:
        print '!!! %s test failed' % nFail

    if of is not None:
        of.stop()
        
    return nFail, nTotal, testTime

##
# 30.05.2007, c
# 01.06.2007
# 04.06.2007
# 19.07.2007
def wrapRunTests( options ):
    def runTests( stats, dirName, fileNames ):
        fileNames = [fileName for fileName in fileNames
                     if (len( fileName ) > 8) and
                     fileName[:5] == 'test_' and fileName[-3:] == '.py']

        print '<<< directory: %s, test files: %d' % (dirName, len( fileNames ))

        for fileName in fileNames:
            confName = op.join( dirName, fileName )

            nFail, nTotal, testTime = runTest( confName, options )

            stats[0] += 1
            stats[1] += nFail
            stats[2] += nTotal
            stats[3] += testTime

    return runTests

usage = """%prog [options] [testFileName]"""

help = {
    'dir' : 'directory with tests [default: %default]',
    'outDir' : 'directory for storing test results and temporary files'
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

    version = open( op.join( init_sfe.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]
    parser = OptionParser( usage = usage, version = "%prog, SFE-" + version )
    parser.add_option( "", "--print-doc",
                       action = "store_true", dest = "printDoc",
                       default = False, help = help['print-doc'] )
    parser.add_option( "-d", "--dir", metavar = 'directory',
                       action = "store", dest = "testDir",
                       default = "tests", help = help['dir'] )
    parser.add_option( "-o", "--output", metavar = 'directory',
                       action = "store", dest = "outDir",
                       default = "output-tests", help = help['outDir'] )
    parser.add_option( "", "--debug",
                       action = "store_true", dest = "debug",
                       default = False, help = help['debug'] )
    parser.add_option( "", "--filter-none",
                       action = "store_true", dest = "filterNone",
                       default = False, help = help['filter-none'] )
    parser.add_option( "", "--filter-less",
                       action = "store_true", dest = "filterLess",
                       default = False, help = help['filter-less'] )
    parser.add_option( "", "--filter-more",
                       action = "store_true", dest = "filterMore",
                       default = False, help = help['filter-more'] )
    options, args = parser.parse_args()

    if options.printDoc:
        print __doc__
        return

    if len( args ) > 1:
        parser.print_help(),
        return
    elif len( args ) == 1:
        testFileName = args[0]

        stats = runTest( testFileName, options )
        print '1 test file executed in %.2f s, %d failure(s) of %d test(s)'\
              % (stats[2], stats[0], stats[1])
        return


    stats = [0, 0, 0, 0.0]
    op.walk( options.testDir, wrapRunTests( options ), stats )
    print '%d test file(s) executed in %.2f s, %d failure(s) of %d test(s)'\
          % (stats[0], stats[3], stats[1], stats[2])

if __name__ == '__main__':
    main()
