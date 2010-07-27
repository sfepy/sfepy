from sfepy.base.base import *
import extmods.fem as fem
from sfepy.terms import Term, DataCaches
from region import Region
from equations import Equation, Equations
from integrals import Integrals
from variables import Variables
from fields import setup_dof_conns, setup_extra_data

##
# 02.10.2007, c
class Evaluator( Struct ):
    pass

##
# 02.10.2007, c
class BasicEvaluator( Evaluator ):
    ##
    # 02.10.2007, c
    def __init__( self, problem, mtx = None, **kwargs ):
        Evaluator.__init__( self, problem = problem, data = kwargs )
        if mtx is None:
            self.mtx = problem.mtx_a
        else:
            self.mtx = mtx

    ##
    # c: 11.04.2008, r: 11.04.2008
    def set_term_args( self, **kwargs ):
        self.data = kwargs

    def eval_residual( self, vec, is_full = False ):
        if not is_full:
            vec = self.make_full_vec( vec )
        try:
            pb = self.problem
            vec_r = pb.equations.eval_residuals(vec)

        except StopIteration, exc:
            status = exc.args[0]
            output( 'error %d in term "%s" of equation "%s"!'\
                    % (status, exc.args[1].name, exc.args[2].desc) )
            raise ValueError

        return vec_r
            
    def eval_tangent_matrix( self, vec, mtx = None, is_full = False ):
        if isinstance( vec, str ) and vec == 'linear':
            return get_default( mtx, self.mtx )
        
        if not is_full:
            vec = self.make_full_vec( vec )
        try:
            pb = self.problem
            if mtx is None:
                mtx = self.mtx
            mtx.data[:] = 0.0
            mtx = pb.equations.eval_tangent_matrices(vec, mtx)

        except StopIteration, exc:
            status = exc.args[0]
            output( ('error %d in term "%s" of derivative of equation "%s"'
                     + ' with respect to variable "%s"!')\
                    % (status,
                       exc.args[1].name, exc.args[2].desc, exc.args[3] ) )
            raise ValueError

        return mtx

    ##
    # 02.12.2005, c
    # 09.12.2005
    # 25.07.2006
    # 02.10.2007
    def update_vec( self, vec, delta ):
        self.problem.update_vec( vec, delta )

    def strip_state_vector( self, vec ):
        return self.problem.equations.strip_state_vector(vec,
                                                         follow_epbc=False)

    def make_full_vec( self, vec ):
        return self.problem.equations.make_full_vec(vec)

##
# 04.10.2007, c
class LCBCEvaluator( BasicEvaluator ):

    ##
    # 04.10.2007, c
    def __init__( self, problem, mtx = None, **kwargs ):
        BasicEvaluator.__init__( self, problem, mtx, **kwargs )
        self.op_lcbc = problem.equations.get_lcbc_operator()

    ##
    # 04.10.2007, c
    def eval_residual( self, vec, is_full = False ):
        if not is_full:
            vec = self.make_full_vec( vec )
        vec_r = BasicEvaluator.eval_residual( self, vec, is_full = True )
        vec_rr = self.op_lcbc.T * vec_r
        return vec_rr
            
    ##
    # 04.10.2007, c
    def eval_tangent_matrix( self, vec, mtx = None, is_full = False ):
        if isinstance( vec, str ) and vec == 'linear':
            return get_default( mtx, self.mtx )

        if not is_full:
            vec = self.make_full_vec( vec )
        mtx = BasicEvaluator.eval_tangent_matrix( self, vec, mtx = mtx,
                                                  is_full = True )
        mtx_r = self.op_lcbc.T * mtx * self.op_lcbc
        mtx_r = mtx_r.tocsr()
        mtx_r.sort_indices()
##         import pylab
##         from sfepy.base.plotutils import spy
##         spy( mtx_r )
##         pylab.show()
##         print mtx_r.__repr__()
        return mtx_r

    ##
    # 04.10.2007, c
    def update_vec( self, vec, delta_r ):
        delta = self.op_lcbc * delta_r
        BasicEvaluator.update_vec( self, vec, delta )

    def strip_state_vector( self, vec ):
        vec = BasicEvaluator.strip_state_vector( self, vec )
        return self.op_lcbc.T * vec

def create_evaluable(expression, fields, materials, variables, integrals,
                     ebcs=None, epbcs=None, lcbcs=None, ts=None, functions=None,
                     auto_init=False, mode='eval', verbose=True, kwargs=None):
    """
    Returns
    -------
    equation : Equation instance
        The equation that is ready to be evaluated.
    variables : Variables instance
        The variables used in the equation.
    """
    if kwargs is None:
        kwargs = {}

    domain = fields[fields.keys()[0]].domain
    caches = DataCaches()

    # Create temporary variables.
    aux_vars = Variables(variables)
    equations = Equations.from_conf({'tmp' : expression},
                                    aux_vars, domain.regions,
                                    materials, integrals,
                                    caches=caches, user=kwargs,
                                    verbose=verbose)
    equations.collect_conn_info()

    # The true variables used in the expression.
    variables = equations.variables
    if auto_init:
        for var in variables:
            var.init_data(step=0)

    if mode == 'weak':
        setup_dof_conns(equations.conn_info)
        materials.time_update(ts, domain, equations, verbose=False)
        equations.time_update(ts, domain.regions,
                              ebcs, epbcs, lcbcs, functions)

    else:
        setup_extra_data(equations.conn_info)
        materials.time_update(ts, domain, equations, verbose=False)

    equations.describe_geometry(verbose=False)

    return equations, variables


def eval_equations(equations, variables, clear_caches=True,
                   mode='eval', dw_mode='vector', term_mode=None):
    asm_obj = None

    if mode == 'weak':
        if dw_mode == 'vector':
            asm_obj = equations.create_stripped_state_vector()

        else:
            asm_obj = equations.create_matrix_graph()

    if clear_caches:
        equations.invalidate_term_caches()

    out = equations.evaluate(mode=mode, dw_mode=dw_mode, term_mode=term_mode,
                             asm_obj=asm_obj)

    if variables.has_lcbc and mode == 'weak':
        op_lcbc = variables.op_lcbc
        if dw_mode == 'vector':
            out = op_lcbc.T * out

        elif dw_mode == 'matrix':
            out = op_lcbc.T * out * op_lcbc
            out = out.tocsr()
            out.sort_indices()

    return out

def evaluate(expression, fields, materials, variables, integrals,
             ebcs=None, epbcs=None, lcbcs=None, ts=None, functions=None,
             auto_init=False, mode='eval', dw_mode='vector', term_mode=None,
             ret_variables=False, verbose=True, kwargs=None):
    """
    Convenience function calling create_evaluable() and eval_equations().

    Parameters
    ----------
    mode : one of 'eval', 'el_avg', 'qp', 'weak'
        The evaluation mode.
    """
    aux = create_evaluable(expression, fields, materials, variables, integrals,
                           ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs, ts=ts,
                           functions=functions, auto_init=auto_init,
                           mode=mode, verbose=verbose, kwargs=kwargs)
    equations, variables = aux

    out = eval_equations(equations, variables,
                         mode=mode, dw_mode=dw_mode, term_mode=term_mode)

    if ret_variables:
        out = (out, variables)

    return out

def assemble_by_blocks(conf_equations, problem, ebcs=None, epbcs=None,
                       dw_mode='matrix'):
    """Instead of a global matrix, return its building blocks as defined in
    `conf_equations`. The name and row/column variables of each block have to
    be encoded in the equation's name, as in:

    conf_equations = {
      'A,v,u' : "dw_lin_elastic_iso.i1.Y2( inclusion.lame, v, u )",
    }

    Notes
    -----
    `ebcs`, `epbcs` must be either lists of BC names, or BC configuration
    dictionaries.
    """
    if isinstance( ebcs, list ) and isinstance( epbcs, list ):
        bc_mode = 0
    elif isinstance( ebcs, dict ) and isinstance( epbcs, dict ):
        bc_mode = 1
    else:
        raise TypeError('bad BC!')

    matrices = {}
    for key, mtx_term in conf_equations.iteritems():
        ks = key.split( ',' )
        mtx_name, var_names = ks[0], ks[1:]
        output( mtx_name, var_names )

        problem.set_equations({'eq': mtx_term})
	variables = problem.get_variables()
	indx = variables.get_indx
	dummy = variables.create_state_vector()

        if bc_mode == 0:
            problem.select_bcs( ebc_names = ebcs, epbc_names = epbcs )

        else:
            problem.time_update( conf_ebc = ebcs, conf_epbc = epbcs )

        ir = indx( var_names[0], stripped = True, allow_dual = True )
        ic = indx( var_names[1], stripped = True, allow_dual = True )

        mtx = problem.evaluate(mtx_term, dummy, dw_mode='matrix',
                               copy_materials=False,
			       update_materials=False)
        matrices[mtx_name] = mtx[ir,ic]

    return matrices
