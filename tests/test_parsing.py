import sfepy.base.testing as tst

def test_parse_equations():
    from sfepy.discrete.parse_equations import create_bnf

    test_strs = [
        """- ev_volume.i1.Omega(uc)""",
        """- 2 * dw_term.i1.Omega(uc)
         = - 3.0 * dw_term2.i1.Omega2(uc)""",
        """2 * dw_term.i1.Omega(uc)
         = 3.0 * dw_term2.i1.Omega2(uc)""",
        """- (2 + 1j - 3) * dw_term.i1.Omega(uc)
           = - (1j - 3.0 + 2.0j) * dw_term2.i1.Omega2(uc)""",
        """(2 + 1j) * dw_term.i1.Omega(uc)
           = (3.0 + 2.0j) * dw_term2.i1.Omega2(uc)""",
        """d_term1.Y(fluid, u, w, Nu, dcf, mode)
             + d_term2.Omega(u, w, Nu, dcf, mode)
             - d_another_term.Elsewhere(w, p, Nu, dcf, mode)
             = - dw_rhs.a.Y3(u, q, Nu, dcf, mode)""",
        """no_args() = 0""",
        """+ something(a, b, c) = + something_else(c, a, d[-1])""",
        """term_.a.a(u)""",
        """term.i1.Omega(v, du/dt) + term2.i2.Gamma(v, dphi/dt)""",
        """dw_jump.isurf.Gamma12_1(jump1.val, q1, p1, tr(p2))""",
        """term_with_trace.2.Gamma1(v, tr(Gamma2, u))""",
    ]

    n_fail = 0
    term_descs = []
    for test_str in test_strs:
        term_descs[:] = []
        try:
            bnf = create_bnf(term_descs)
            bnf.parseString(test_str)
        except:
            tst.report('failed: %s' % test_str)
            n_fail += 1

    tst.report('%d failure(s)' % n_fail)

    assert n_fail == 0

def test_parse_regions():
    from sfepy.discrete.parse_regions import create_bnf

    test_strs = ['vertices of surface -v r.Omega',
                 'r.Y_2 +v copy r.Y_1',
                 'vertices in (y <= 0.00001) & (x < 0.11)',
                 'vertices in ((y <= 0.00001) & (x < 0.11))',
                 'vertices in (((y <= 0.00001) & (x < 0.11)))',
                 'vertices in (((0.00001 < y) & (x < 0.11)))',
                 'vertices of group 0',
                 """vertices of group 0 +v vertices of group 1
                 +c cells by afun""",
                 'all -v vertices in (y == 0.00001)',
                 'all -v vertices of surface',
                 'all -c r.DOmega_100',
                 """r.Y_1 -v vertices of surface *c r.Z_8
                    *v vertices in (y > 0)""",
                 'vertices of surface +v vertices by pokus',
                 'cells of group 6 +c vertices by fn2_3c',
                 """r.Y_1 *v (r.Y_2 +c (vertices in (y > 0) *v r.Y_32))
                    -v vertices of surface -c r.Y_5""",
                 'copy r.ab2-b-c +v r.d12-23',
                 'vertices by afun',
                 'vertex in r.Gamma_3',
                 'vertex 10',
                 'vertex 10, 20, 30',
                 'cell 10',
                 'cell 10, 20, 30',
                 'vertex 10, 20 +v cell 30, 40',
                 '(vertex 10, 20) +v (cell 30, 40)',
                 'cell 10, 20, 30 +v vertex 10',
                 'cells by afun']

    stack = []
    bnf = create_bnf(stack)

    n_fail = 0
    for test_str in test_strs:
        stack[:] = []

        try:
            bnf.parseString(test_str)

        except:
            tst.report('failed: %s' % test_str)
            n_fail += 1

    tst.report('%d failures' % n_fail)

    assert n_fail == 0
