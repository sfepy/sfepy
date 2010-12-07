from sfepy.base.base import *
from sfepy.terms import DataCaches
from equations import Equations
from variables import Variables
from fields import setup_dof_conns, setup_extra_data

##
# 02.10.2007, c
class Evaluator( Struct ):
    pass

##
# 02.10.2007, c
class BasicEvaluator( Evaluator ):

    def __init__(self, problem, matrix_hook=None):
        Evaluator.__init__(self, problem=problem,
                           matrix_hook=matrix_hook)

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
            return get_default(mtx, self.problem.mtx_a)

        if not is_full:
            vec = self.make_full_vec( vec )
        try:
            pb = self.problem
            if mtx is None:
                mtx = pb.mtx_a
            mtx.data[:] = 0.0
            mtx = pb.equations.eval_tangent_matrices(vec, mtx)

        except StopIteration, exc:
            status = exc.args[0]
            output( ('error %d in term "%s" of derivative of equation "%s"'
                     + ' with respect to variable "%s"!')\
                    % (status,
                       exc.args[1].name, exc.args[2].desc, exc.args[3] ) )
            raise ValueError

        if self.matrix_hook is not None:
            mtx = self.matrix_hook(mtx, self.problem, call_mode='basic')

        return mtx

    def make_full_vec( self, vec ):
        return self.problem.equations.make_full_vec(vec)

##
# 04.10.2007, c
class LCBCEvaluator( BasicEvaluator ):

    ##
    # 04.10.2007, c
    def __init__(self, problem, matrix_hook=None):
        BasicEvaluator.__init__(self, problem, matrix_hook=matrix_hook)
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
            return get_default(mtx, self.problem.mtx_a)

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

        if self.matrix_hook is not None:
            mtx_r = self.matrix_hook(mtx_r, self.problem, call_mode='lcbc')

        return mtx_r

def create_evaluable(expression, fields, materials, variables, integrals,
                     update_materials=True,
                     ebcs=None, epbcs=None, lcbcs=None, ts=None, functions=None,
                     auto_init=False, mode='eval', extra_args=None,
                     verbose=True, kwargs=None):
    """
    Create evaluable object (equations and corresponding variables)
    from the `expression` string.

    Parameters
    ----------
    expression : str
        The expression to evaluate.
    fields : dict
        The dictionary of fields used in `variables`.
    materials : Materials instance
        The materials used in the expression.
    variables : Variables instance
        The variables used in the expression.
    integrals : Integrals instance
        The integrals to be used.
    update_materials : bool
        Call time update function of the materials. Safe but can be slow.
    ebcs : Conditions instance, optional
        The essential (Dirichlet) boundary conditions for 'weak'
        mode.
    epbcs : Conditions instance, optional
        The periodic boundary conditions for 'weak'
        mode.
    lcbcs : Conditions instance, optional
        The linear combination boundary conditions for 'weak'
        mode.
    ts : TimeStepper instance, optional
        The time stepper.
    functions : Functions instance, optional
        The user functions for boundary conditions, materials
        etc.
    auto_init : bool
        Set values of all variables to all zeros.
    mode : one of 'eval', 'el_avg', 'qp', 'weak'
        The evaluation mode - 'weak' means the finite element
        assembling, 'qp' requests the values in quadrature points,
        'el_avg' element averages and 'eval' means integration over
        each term region.
    extra_args : dict, optional
        Extra arguments to be passed to terms in the expression.
    verbose : bool
        If False, reduce verbosity.
    kwargs : dict, optional
        The variables (dictionary of (variable name) : (Variable
        instance)) to be used in the expression.

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

    if extra_args is None:
        extra_args = kwargs

    else:
        extra_args = copy(extra_args)
        extra_args.update(kwargs)

    equations = Equations.from_conf({'tmp' : expression},
                                    aux_vars, domain.regions,
                                    materials, integrals,
                                    setup=False,
                                    caches=caches, user=extra_args,
                                    verbose=verbose)
    equations.collect_conn_info()
    equations.assign_geometries()

    # The true variables used in the expression.
    variables = equations.variables
    if auto_init:
        for var in variables:
            var.init_data(step=0)

    # The true materials used in the expression.
    materials = equations.materials

    if mode == 'weak':
        setup_dof_conns(equations.conn_info)
        if update_materials:
            materials.time_update(ts, domain, equations, verbose=False)
        equations.time_update(ts, ebcs, epbcs, lcbcs, functions)

    else:
        setup_extra_data(equations.conn_info)
        if update_materials:
            materials.time_update(ts, domain, equations, verbose=False)

    return equations, variables


def eval_equations(equations, variables, clear_caches=True,
                   mode='eval', dw_mode='vector', term_mode=None):
    """
    Evaluate the equations.

    Parameters
    ----------
    equations : Equations instance
        The equations returned by :func:`create_evaluable()`.
    variables : Variables instance
        The variables returned by :func:`create_evaluable()`.
    clear_caches : bool
        If True, clear term caches.
    mode : one of 'eval', 'el_avg', 'qp', 'weak'
        The evaluation mode - 'weak' means the finite element
        assembling, 'qp' requests the values in quadrature points,
        'el_avg' element averages and 'eval' means integration over
        each term region.
    dw_mode : 'vector' or 'matrix'
        The assembling mode for 'weak' evaluation mode.
    term_mode : str
        The term call mode - some terms support different call modes
        and depending on the call mode different values are
        returned.

    Returns
    -------
    out : array
        The result of the evaluation.
    """
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
             update_materials=True,
             ebcs=None, epbcs=None, lcbcs=None, ts=None, functions=None,
             auto_init=False, mode='eval', dw_mode='vector', term_mode=None,
             ret_variables=False, extra_args=None, verbose=True, kwargs=None):
    """
    Convenience function calling :func:`create_evaluable()` and
    :func:`eval_equations()`.

    Parameters
    ----------
    ... : arguments
        See docstrings of :func:`create_evaluable()` and
        :func:`eval_equations()`.
    """
    aux = create_evaluable(expression, fields, materials, variables, integrals,
                           update_materials=update_materials,
                           ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs, ts=ts,
                           functions=functions, auto_init=auto_init,
                           mode=mode, extra_args=extra_args,
                           verbose=verbose, kwargs=kwargs)
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

        if bc_mode == 0:
            problem.select_bcs( ebc_names = ebcs, epbc_names = epbcs )

        else:
            problem.time_update(ebcs=ebcs, epbcs=epbcs)

        ir = indx( var_names[0], stripped = True, allow_dual = True )
        ic = indx( var_names[1], stripped = True, allow_dual = True )

        mtx = problem.evaluate(mtx_term, auto_init=True,
                               mode='weak', dw_mode='matrix',
                               copy_materials=False)
        matrices[mtx_name] = mtx[ir,ic]

    return matrices
