# 10.07.2007, c
# last revision: 20.02.2008

fileName_meshes = ['database/kostka_medium_tetra.mesh',
                   'database/kostka_medium_tetra.mesh',
                   'database/kostka_medium.mesh']
allYourBases = [{'Omega' : '3_4_P1'},
                {'Omega' : '3_4_P2'},
                {'Omega' : '3_8_Q1'}]

fileName_mesh = None

field_1 = {
    'name' : '3_displacement',
    'dim' : (3,1),
    'flags' : (),
    'domain' : 'Omega',
    'bases' : None
}

material_1 = {
    'name' : 'solid',
    'mode' : 'here',
    'region' : 'Omega',
    'lame' : {'lambda' : 1e1, 'mu' : 1e0},
}

material_2 = {
    'name' : 'spring',
    'mode' : 'here',
    'region' : 'Omega',
    'pars' : {'stiffness' : 1e0, 'projection' : None},
}

variable_1 = {
    'name' : 'u',
    'kind' : 'unknown field',
    'field' : '3_displacement',
    'dofs' : (0, 1, 2),
    'order' : 0,
}
variable_2 = {
    'name' : 'v',
    'kind' : 'test field',
    'field' : '3_displacement',
    'dofs' : (0, 1, 2),
    'dual' : 'u',
}

region_1000 = {
    'name' : 'Omega',
    'select' : 'all',
}

region_1 = {
    'name' : 'Bottom',
    'select' : 'nodes in (z < -0.499)',
}
region_2 = {
    'name' : 'Top',
    'select' : 'nodes in (z > 0.499)',
}

ebc_1 = {
    'name' : 'Load',
    'region' : 'Top',
    'dofs' : (2,),
    'value' : 0.1,
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d3',
}

equations = {
    'balance_of_forces' :
    """dw_sdcc.i1.Omega( solid.lame, v, u )
     = dw_point_lspring.i1.Bottom( spring.pars, v, u )""",
}

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.umfpack',
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'iMax'      : 1,
    'epsA'      : 1e-10,
    'epsR'      : 1.0,
    'macheps'   : 1e-16,
    'linRed'    : 1e-2, # Linear system error < (epsA * linRed).
    'lsRed'     : 0.1,
    'lsRedWarp' : 0.001,
    'lsOn'      : 1.1,
    'lsMin'     : 1e-5,
    'check'     : 0,
    'delta'     : 1e-6,
    'isPlot'    : False,
    'linSolver' : 'umfpack',
    'matrix'    : 'internal', # 'external' or 'internal'
    'problem'   : 'nonlinear', # 'nonlinear' or 'linear' (ignore iMax)
}

##
# FE assembling parameters.
fe = {
    'chunkSize' : 1000
}

from sfe.base.testing import TestCommon

##
# 10.07.2007, c
class Test( TestCommon ):

    ##
    # 10.07.2007, c
    def fromConf( conf, options ):
        return Test( conf = conf, options = options )
    fromConf = staticmethod( fromConf )

    ##
    # c: 10.07.2007, r: 08.02.2008
    def test_linear_terms( self ):
        from sfe.solvers.generic import solveStationary
        from sfe.base.base import IndexedStruct
        import os.path as op

        ok = True
        for ii, bases in enumerate( allYourBases ):
            fname = fileName_meshes[ii]

            self.conf.fileName_mesh = fname
            self.conf.fields['field_1'].bases = bases
            self.report( 'mesh: %s, base: %s' % (fname, bases) )
            status = IndexedStruct()
            problem, vec, data = solveStationary( self.conf,
                                                  nlsStatus = status )
            converged = status.condition == 0
            ok = ok and converged
            self.report( 'converged: %s' % converged )

            name = op.join( self.options.outDir,
                            '_'.join( ('test_elasticity_small_strain',
                                      op.splitext( op.basename( fname ) )[0],
                                      bases.values()[0] ))
                            + '.vtk' )
            problem.saveState( name, vec )

##             trunk = op.join( self.options.outDir,
##                              op.splitext( op.basename( fname ) )[0] )
##             problem.saveFieldMeshes( trunk )
##             problem.saveRegions( trunk )
            
        return ok
