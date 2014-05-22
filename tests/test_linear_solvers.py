from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/special/cube_cylinder.mesh'

if 0:
    from sfepy.discrete.fem.utils import refine_mesh
    refinement_level = 1
    filename_mesh = refine_mesh(filename_mesh, refinement_level)

material_2 = {
    'name' : 'coef',
    'values' : {'val' : 1.0},
}

field_1 = {
    'name' : 'temperature',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

variables = {
    't' : ('unknown field', 'temperature', 0),
    's' : ('test field',    'temperature', 't'),
}

regions = {
    'Omega' : 'all',
    'Gamma_Left' : ('vertices in (x < 0.0001)', 'facet'),
    'Gamma_Right' : ('vertices in (x > 0.999)', 'facet'),
}

ebcs = {
    't1' : ('Gamma_Left', {'t.0' : 2.0}),
    't2' : ('Gamma_Right', {'t.0' : -2.0}),
}

integral_1 = {
    'name' : 'i',
    'order' : 1,
}

equations = {
    'Temperature' : """dw_laplace.i.Omega( coef.val, s, t ) = 0"""
}

solvers = {
    'd00' : ('ls.umfpack',
             {'warn' : True,}
    ),
##     'd01' : ('ls.scipy_direct',
##              {'method' : 'superlu',
##               'warn' : True,}
##     ),
    'i00' : ('ls.pyamg',
             {'method' : 'ruge_stuben_solver',
              'accel' : 'cg',
              'eps_r'   : 1e-12,}
    ),
    'i01' : ('ls.pyamg',
             {'method' : 'smoothed_aggregation_solver',
              'accel' : 'cg',
              'eps_r'   : 1e-12,}
    ),
    'i10' : ('ls.petsc',
             {'method' : 'cg', # ksp_type
              'precond' : 'icc', # pc_type
              'eps_a' : 1e-12, # abstol
              'eps_r' : 1e-12, # rtol
              'i_max' : 1000,} # maxits
    ),
    'i20' : ('ls.scipy_iterative',
             {'method' : 'cg',
              'i_max'   : 1000,
              'eps_r'   : 1e-12,}
    ),
    'i21' : ('ls.scipy_iterative',
             {'method' : 'bicgstab',
              'i_max'   : 1000,
              'eps_r'   : 1e-12,}
    ),
    'i22' : ('ls.scipy_iterative',
             {'method' : 'qmr',
              'i_max'   : 1000,
              'eps_r'   : 1e-12,}
    ),

    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}


options = {
    'nls' : 'newton',
}

from sfepy.base.testing import TestCommon
output_name = 'test_linear_solvers_%s.vtk'

##
# c: 02.05.2008
class Test( TestCommon ):
    can_fail = ['ls.pyamg', 'ls.petsc']

    ##
    # c: 02.05.2008, r: 02.05.2008
    def from_conf( conf, options ):
        from sfepy.discrete import Problem

        problem = Problem.from_conf(conf)
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
        from sfepy.base.base import IndexedStruct
        import os.path as op

        solver_confs = self._list_linear_solvers( self.problem.solver_confs )

        ok = True
        tt = []
        for solver_conf in solver_confs:
            method = solver_conf.get('method', '')
            precond = solver_conf.get('precond', '')
            name = ' '.join( (solver_conf.name, solver_conf.kind,
                              method, precond) ).rstrip()
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
                tt.append( [name, status.time_stats['solve'], status.err] )

                aux = name.replace(' ', '_')
                fname = op.join( self.options.out_dir,
                                op.split( self.conf.output_name )[1] ) % aux
                self.problem.save_state( fname, state )
            else:
                self.report( 'solver failed:' )
                self.report( exc )
                tt.append( [name, 1e10, 1e10] )


        tt.sort( cmp = lambda a, b: cmp( a[1], b[1] ) )
        self.report( 'solution times (rezidual norms):' )
        for row in tt:
            self.report( '%.2f [s]' % row[1], '(%.3e)' % row[2], ':', row[0] )

        return ok
