from sfepy.base.testing import TestCommon

##
# 16.07.2007, c
class Test( TestCommon ):

    ##
    # 16.07.2007, c
    def from_conf( conf, options ):
        return Test( conf = conf, options = options )
    from_conf = staticmethod( from_conf )

    ##
    # c: 16.07.2007, r: 08.07.2008
    def test_parse_equations( self ):
        from sfepy.fem.parseEq import create_bnf

        test_strs = [
            """- d_volume.i1.Omega( uc )""",
            """2 * dw_term.i1.Omega( uc ) = - 3.0 * dw_term2.i1.Omega2( uc )""",
            """d_term1.Y( fluid, u, w, Nu, dcf, mode )
                 + d_term2.Omega( u, w, Nu, dcf, mode )
                 - d_another_term.Elsewhere( w, p, Nu, dcf, mode )
                 = - dw_rhs.Y3.a( u, q, Nu, dcf, mode )""",
            """no_args() = 0""",
            """+ something( a, b, c ) = + something_else( c, a, d[-1] )""",
            """term_.a.a( u )""",
            """term.i1.Omega( v, du/dt ) + term2.i2.Gamma( v, dphi/dt)""",
            """dw_jump.isurf.Gamma12_1( jump1.val, q1, p1, tr(p2) )""",
        ]

        n_fail = 0
        term_descs = []
        for test_str in test_strs:
            term_descs[:] = []
            try:
                bnf = create_bnf( term_descs, {} )
                bnf.parseString( test_str )
            except:
                self.report( 'failed: %s' % test_str )
                if self.options.debug:
                    raise
                n_fail += 1
            for td in term_descs:
                print td
        self.report( '%d failure(s)' % n_fail )

        if n_fail:
            raise AssertionError
        return True

    ##
    # c: 16.07.2007, r: 14.07.2008
    def test_parse_regions( self ):
        from sfepy.fem.parseReg import create_bnf, _test_strs

        test_strs = ['nodes of surface -n r.Omega',
                     'r.Y_2 +n copy r.Y_1',
                     'nodes in (y <= 0.00001) & (x < 0.11)',
                     'nodes in ((y <= 0.00001) & (x < 0.11))',
                     'nodes in (((y <= 0.00001) & (x < 0.11)))',
                     'nodes in (((0.00001 < y) & (x < 0.11)))',
                     'nodes of group 0',
                     """nodes of group 0 +n nodes of group 1
                     +e elements by afun( domain )""",
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
                     'node 10',
                     'elements by afun( domain )']

        stack = []
        bnf = create_bnf( stack )

        n_fail = 0
        for test_str in test_strs:
            stack[:] = []

            try:
                out = bnf.parseString( test_str )
            except:
                self.report( 'failed: %s' % test_str )
                n_fail += 1

        self.report( '%d failures' % n_fail )

        if n_fail:
            raise AssertionError
        return True
