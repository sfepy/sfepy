# c: 19.05.2008, r: 19.05.2008
fileName_mesh = 'database/phono/mesh_circ21.mesh'

is3D = False

if is3D:
    dim, geom = 3, '3_4'
else:
    dim, geom = 2, '2_3'

field_1 = {
    'name' : 'scalar_field',
    'dim' : (1,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '%s_P1' % geom}
}

field_2 = {
    'name' : 'vector_field',
    'dim' : (dim,1),
    'domain' : 'Omega',
    'bases' : {'Omega' : '%s_P1' % geom}
}

variables = {
    'us'  : ('unknown field',   'scalar_field', 0),
    'ts'  : ('test field',      'scalar_field', 'us'),
    'ps1' : ('parameter field', 'scalar_field', 'us'),
    'ps2' : ('parameter field', 'scalar_field', 'us'),
    'uv'  : ('unknown field',   'vector_field', 0),
    'tv'  : ('test field',      'vector_field', 'uv'),
    'pv1' : ('parameter field', 'vector_field', 'uv'),
    'pv2' : ('parameter field', 'vector_field', 'uv'),
}

regions = {
    'Omega' : ('all', {}),
}

integral_1 = {
    'name' : 'i1',
    'kind' : 'v',
    'quadrature' : 'gauss_o2_d2',
}

material_1 = {
    'name' : 'm',
    'mode' : 'function',
    'region' : 'Omega',
    'function' : 'getPars',
}   

fe = {
    'chunkSize' : 100
}

##
# c: 19.05.2008, r: 19.05.2008
def getPars( ts, coor, region, ig, mode = None ):
    nNod, dim = coor.shape
    sym = (dim + 1) * dim / 2

    if mode == 'biot':
        val = nm.zeros( (sym,), dtype = nm.float64 )
        val[:dim] = 0.132
        val[dim:sym] = 0.092
    elif mode == 'biotM':
        val = 1.0 / nm.array( [3.8], dtype = nm.float64 )
    elif mode == 'permeability':
        val = nm.eye( dim, dtype = nm.float64 )
    else:
        raise ValueError

    return {'val' : val}

# ('d' variables, 'dw' variables (test must be paired with unknown, which
# should be last!), mat mode)
testTerms = {
    '%s_biot_div.i1.Omega( m.val, %s, %s )' :
    (('ps1', 'pv1'), ('ts', 'pv1', 'us'), 'biot'),
    '%s_diffusion.i1.Omega( m.val, %s, %s )' :
    (('ps1', 'ps2'), ('ts', 'ps1', 'us'), 'permeability'),
    '%s_volume_wdot.i1.Omega( m.val, %s, %s )' :
    (('ps1', 'ps2'), ('ts', 'ps1', 'us'), 'biotM'),
}

import numpy as nm
from sfepy.base.testing import TestCommon
from sfepy.base.base import debug, pause

##
# c: 19.05.2008
class Test( TestCommon ):

    ##
    # c: 19.05.2008, r: 19.05.2008
    def fromConf( conf, options ):
        from sfepy.fem.problemDef import ProblemDefinition

        problem = ProblemDefinition.fromConf( conf, initVariables = False )
        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    fromConf = staticmethod( fromConf )

    ##
    # c: 19.05.2008, r: 19.05.2008
    def test_consistency_d_dw( self ):
        from sfepy.base.base import selectByNames
        from sfepy.fem.evaluate import evalTermOP

        ok = True
        pb = self.problem
        for termTemplate, (dVars, dwVars, parMode) in testTerms.iteritems():
            print termTemplate, dVars, dwVars, parMode

            varNames = nm.unique1d( dVars + dwVars )
            variables = selectByNames( pb.conf.variables, varNames )
            pb.setVariables( variables )

            term1 = termTemplate % (('d',) + dVars)
            pb.setEquations( {'eq': term1} )

            matArgs = {'m' : {'mode' : parMode}} 
            pb.timeUpdate( extraMatArgs = matArgs )

            vecs = {}
            for varName in dVars:
                var = pb.variables[varName]
                nDof = var.field.nNod * var.field.dim[0]
                vecs[varName] = nm.arange( nDof, dtype = nm.float64 )
                var.dataFromData( vecs[varName] )
            dummy = pb.createStateVector()
            val1 = evalTermOP( dummy, term1, pb )
            self.report( '%s: %s' % (term1, val1) )
            
            term2 = termTemplate % (('dw',) + dwVars[:-1])

            vec = evalTermOP( dummy, term2, pb )
            val2 = nm.dot( vecs[dVars[0]], vec )
            self.report( '%s: %s' % (term2, val2) )

            err = nm.abs( val1 - val2 ) / nm.abs( val1 )
            _ok = err < 1e-12
            self.report( 'relative difference: %e -> %s' % (err, _ok) )

            ok = ok and _ok
            
        return ok
