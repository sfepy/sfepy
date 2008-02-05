#!/usr/bin/env python
import sys
import time
import os
import os.path as op
from optparse import OptionParser

import init_sfe
from sfe.base.conf import ProblemConf

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
# c: 04.06.2007, r: 05.02.2008
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
    required = ['Test']
    other = ['fileName_mesh', 'field_[0-9]+', 'ebc|nbc', 'fe',
             'equations', 'solution', 'epbc', 'lcbc',
             'region_[0-9]+', 'variables', 'material_[0-9]+', 'solver_[0-9]+',
             'functions', 'modules']

    num = 1
    testTime = 0.0
    try:
        conf = ProblemConf.fromFile( confName, required, other )
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
    'debug' : 'raise silenced exceptions to see what was wrong [default:'
    ' %default]',
    'filter-none' : 'do not filter any messages [default:'
    ' %default]',
    'filter-less' : 'filter output (suppress all except test messages) [default:'
    ' %default]',
    'filter-more' : 'filter output (suppress all except test result messages)'
    ' [default: %default]',
}

##
# c: 30.05.2007, r: 05.02.2008
def main():

    version = open( op.join( init_sfe.install_dir,
                             'VERSION' ) ).readlines()[0][:-1]
    parser = OptionParser( usage = usage, version = "%prog, SFE-" + version )
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
