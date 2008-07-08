from sfepy.base.testing import TestCommon

##
# 16.07.2007, c
class Test( TestCommon ):

    ##
    # 16.07.2007, c
    def fromConf( conf, options ):
        return Test( conf = conf, options = options )
    fromConf = staticmethod( fromConf )

    ##
    # c: 16.07.2007, r: 08.07.2008
    def test_parseEquations( self ):
        from sfepy.fem.parseEq import createBNF

        testStrs = [
            """- d_volume.i1.Omega( uc )""",
            """2 * dw_term.i1.Omega( uc ) = - 3.0 * dw_term2.i1.Omega2( uc )""",
            """d_term1.Y( fluid, u, w, Nu, dcf, mode )
                 + d_term2.Omega( u, w, Nu, dcf, mode )
                 - d_another_term.Elsewhere( w, p, Nu, dcf, mode )
                 = - dw_rhs.Y3.a( u, q, Nu, dcf, mode )""",
            """noArgs() = 0""",
            """+ something( a, b, c ) = + somethingElse( c, a, d )""",
            """term_.a.a( u )"""
        ]

        nFail = 0
        termDescs = []
        for testStr in testStrs:
            termDescs[:] = []
            try:
                bnf = createBNF( termDescs, {} )
                bnf.parseString( testStr )
            except:
                self.report( 'failed: %s' % testStr )
                if self.options.debug:
                    raise
                nFail += 1
            for td in termDescs:
                print td
        self.report( '%d failure(s)' % nFail )

        if nFail:
            raise AssertionError
        return True

    ##
    # 16.07.2007, c
    # 31.07.2007
    def test_parseRegions( self ):
        from sfepy.fem.parseReg import createBNF, _testStrs

        testStrs = ['nodes of surface -n r.Omega',
                    'r.Y_2 +n copy r.Y_1',
                    'nodes in (y <= 0.00001) & (x < 0.11)',
                    'nodes in ((y <= 0.00001) & (x < 0.11))',
                    'nodes in (((y <= 0.00001) & (x < 0.11)))',
                    'nodes in (((0.00001 < y) & (x < 0.11)))',
                    'all -n nodes in (y == 0.00001)',
                    'all -n nodes of surface',
                    'all -e r.DOmega_100',
                    'r.Y_1 -n nodes of surface *e r.Z_8 *n nodes in (y > 0)',
                    'nodes of surface +n nodes by pokus( x, y, z )',
                    'elements of group 6 +e nodes by fn2_3c( x )',
                    """r.Y_1 *n (r.Y_2 +e (nodes in (y > 0) *n r.Y_32))
                    -n nodes of surface -e r.Y_5""",
                    'nodes by noargs()',
                    'nodes by extraargs( x, y, z, abc,3 )',
                    'node in r.Gamma_3',
                    'node 10']

        stack = []
        bnf = createBNF( stack )

        nFail = 0
        for testStr in testStrs:
            stack[:] = []

            try:
                out = bnf.parseString( testStr )
            except:
                self.report( 'failed: %s' % testStr )
                nFail += 1

        self.report( '%d failures' % nFail )

        if nFail:
            raise AssertionError
        return True
