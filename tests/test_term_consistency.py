# c: 19.05.2008, r: 19.05.2008
from sfepy import top_dir

filename_mesh = top_dir + '/meshes/2d/special/circle_in_square.mesh'

is3d = False

if is3d:
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
    'function' : 'get_pars',
}   

fe = {
    'chunk_size' : 100
}

##
# c: 19.05.2008, r: 19.05.2008
def get_pars( ts, coor, mode=None, region=None, ig=None, term = None ):
    if mode == 'qp':
        n_nod, dim = coor.shape
        sym = (dim + 1) * dim / 2

        if term == 'biot':
            val = nm.zeros( (sym, 1), dtype = nm.float64 )
            val[:dim] = 0.132
            val[dim:sym] = 0.092
        elif term == 'biot_m':
            val = 1.0 / nm.array( [3.8], dtype = nm.float64 )
        elif term == 'permeability':
            val = nm.eye( dim, dtype = nm.float64 )
        else:
            raise ValueError

        return {'val' : nm.tile(val, (coor.shape[0], 1, 1))}

functions = {
    'get_pars' : (get_pars,),
}


# (eval term prefix, parameter corresponding to test variable, 'd' variables,
# 'dw' variables (test must be paired with unknown, which should be last!), mat
# mode)
test_terms = [
    ('%s_biot.i1.Omega( m.val, %s, %s )',
     ('dw', 'ps1', ('pv1', 'ps1'), ('pv1', 'ts', 'us', 'uv'), 'biot')),
    ('%s_biot.i1.Omega( m.val, %s, %s )',
     ('dw', 'pv1', ('pv1', 'ps1'), ('tv', 'ps1', 'uv', 'us'), 'biot')),
    ('%s_diffusion.i1.Omega( m.val, %s, %s )',
     ('dw', 'ps1', ('ps1', 'ps2'), ('ts', 'ps1', 'us'), 'permeability')),
    ('%s_volume_wdot.i1.Omega( m.val, %s, %s )',
     ('dw', 'ps1', ('ps1', 'ps2'), ('ts', 'ps1', 'us'), 'biot_m')),
]

import numpy as nm
from sfepy.base.testing import TestCommon
from sfepy.base.base import debug, pause

##
# c: 19.05.2008
class Test( TestCommon ):

    ##
    # c: 19.05.2008, r: 19.05.2008
    def from_conf( conf, options ):
        from sfepy.fem import ProblemDefinition

        problem = ProblemDefinition.from_conf( conf, init_variables = False )
        test = Test( problem = problem,
                     conf = conf, options = options )
        return test
    from_conf = staticmethod( from_conf )

    ##
    # c: 19.05.2008, r: 19.05.2008
    def test_consistency_d_dw( self ):
        from sfepy.base.base import select_by_names
        from sfepy.fem import eval_term_op, Function

        ok = True
        pb = self.problem
        for aux in test_terms:
            term_template, (prefix, par_name, d_vars, dw_vars, mat_mode) = aux
            print term_template, prefix, par_name, d_vars, dw_vars, mat_mode

            var_names = nm.unique1d( d_vars + dw_vars )
            variables = select_by_names( pb.conf.variables, var_names )
            pb.set_variables( variables )

            vecs = {}
            for var_name in d_vars:
                var = pb.variables[var_name]
                n_dof = var.field.n_nod * var.field.shape[0]
                vecs[var_name] = nm.arange( n_dof, dtype = nm.float64 )
                var.data_from_data( vecs[var_name] )

            term1 = term_template % ((prefix,) + d_vars)
            pb.set_equations( {'eq': term1} )

            pb.materials['m'].function.set_extra_args(term = mat_mode)
            pb.time_update()

            dummy = pb.create_state_vector()
            if prefix == 'd':
                val1 = eval_term_op( dummy, term1, pb )
            else:
                val1 = eval_term_op( dummy, term1, pb, call_mode = 'd_eval' )
                
            self.report( '%s: %s' % (term1, val1) )
            
            term2 = term_template % (('dw',) + dw_vars[:2])

            vec = eval_term_op( dummy, term2, pb )
            pvec = pb.variables.get_state_part_view(vec, dw_vars[2])
            val2 = nm.dot( vecs[par_name], pvec )
            self.report( '%s: %s' % (term2, val2) )

            err = nm.abs( val1 - val2 ) / nm.abs( val1 )
            _ok = err < 1e-12
            self.report( 'relative difference: %e -> %s' % (err, _ok) )

            ok = ok and _ok
            
        return ok
