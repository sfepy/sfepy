# c: 02.05.2008, r: 02.05.2008
file_name_mesh = 'database/phono/cube_cylinder.mesh'

material_2 = {
    'name' : 'coef',
    'mode' : 'here',
    'region' : 'Omega',
    'val' : 1.0,
}

field_1 = {
    'name' : 'temperature',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '3_4_P1'}
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

regions = {
    'Omega' : ('all', {}),
    'Gamma_Left' : ('nodes in (x < 0.0001)', {}),
    'Gamma_Right' : ('nodes in (x > 0.999)', {}),
}

ebcs = {
    't1' : ('Gamma_Left', {'t.0' : 2.0}),
    't2' : ('Gamma_Right', {'t.0' : -2.0}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o1_d3',
}

equations = {
    'Temperature' : """dw_laplace.i1.Omega( coef.val, s, t ) = 0"""
}

solver_100 = {
    'name' : 'dls100',
    'kind' : 'ls.umfpack',
}

solver_101 = {
    'name' : 'ls101',
    'kind' : 'ls.pyamg',

    'method' : 'ruge_stuben_solver',
    'eps_a'   : 1e-12,
}

solver_102 = {
    'name' : 'ls102',
    'kind' : 'ls.pyamg',

    'method' : 'smoothed_aggregation_solver',
    'eps_a'   : 1e-12,
}

solver_110 = {
    'name' : 'ls110',
    'kind' : 'ls.scipy_iterative',

    'method' : 'cg',
    'i_max'   : 1000,
    'eps_a'   : 1e-12,
}

solver_111 = {
    'name' : 'ls111',
    'kind' : 'ls.scipy_iterative',

    'method' : 'bicgstab',
    'i_max'   : 1000,
    'eps_a'   : 1e-12,
}

solver_112 = {
    'name' : 'ls112',
    'kind' : 'ls.scipy_iterative',

    'method' : 'qmr',
    'i_max'   : 1000,
    'eps_a'   : 1e-12,
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 1,
    'eps_a'      : 1e-10,
    'eps_r'      : 1.0,
    'macheps'   : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp' : 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'is_plot'    : False,
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore i_max)
}

options = {
    'nls' : 'newton',
}

fe = {
    'chunk_size' : 100000
}

from sfepy.base.testing import TestCommon
output_name = 'test_linear_solvers_%s.vtk'

##
# c: 02.05.2008
class Test( TestCommon ):
    can_fail = ['ls.pyamg']

    ##
    # c: 02.05.2008, r: 02.05.2008
    def from_conf( conf, options ):
        from sfepy.fem.problemDef import ProblemDefinition

        problem = ProblemDefinition.from_conf( conf )
        problem.time_update()

        test = Test( problem = problem, 
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 02.05.2008, r: 02.05.2008
    def _list_linear_solvers( self, confs ):
        d = []
        for key, val in confs.iteritems():
            if val.kind.find( 'ls.' ) == 0:
                d.append( val )
        d.sort( cmp = lambda a, b: cmp( a.name, b.name ) )

        return d

    ##
    # c: 02.05.2008, r: 07.05.2008
    def test_solvers( self ):
        from sfepy.solvers.generic import solve_stationary
        from sfepy.base.base import IndexedStruct
        import os.path as op

        solver_confs = self._list_linear_solvers( self.problem.solver_confs )

        ok = True
        tt = []
        for solver_conf in solver_confs:
            if hasattr( solver_conf, 'method' ):
                method = solver_conf.method
            else:
                method = ''
            name = ' '.join( (solver_conf.name, solver_conf.kind, method) )
            self.report( name )
            self.report( 'matrix size:', self.problem.mtx_a.shape )
            self.report( '        nnz:', self.problem.mtx_a.nnz )
            status = IndexedStruct()
            try:
                self.problem.init_solvers( nls_status = status,
                                          ls_conf = solver_conf )
                state = self.problem.solve()
                failed = status.condition != 0
##                 self.problem.mtx_a.save( 'mtx_laplace_cube',
##                                         format='%d %d %.12e\n' )
            except Exception, exc:
                failed = True
                status = None

            ok = ok and ((not failed) or (solver_conf.kind in self.can_fail))

            if status is not None:
                for kv in status.time_stats.iteritems():
                    self.report( '%10s: %7.2f [s]' % kv )
                self.report( 'condition: %d, err0: %.3e, err: %.3e'\
                             % (status.condition, status.err0, status.err) )
                tt.append( [name, status.time_stats['solve']] )

                fname = op.join( self.options.out_dir,
                                op.split( self.conf.output_name )[1] ) % name
                self.problem.save_state( fname, state )
            else:
                self.report( 'solver failed:' )
                self.report( exc )
                tt.append( [name, 1e10] )


        tt.sort( cmp = lambda a, b: cmp( a[1], b[1] ) )
        self.report( 'solution times:' )
        for row in tt:
            self.report( '%.2f' % row[1], ':', row[0] )

        return ok
