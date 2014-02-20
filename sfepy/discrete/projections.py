"""
Construct projections between FE spaces.
"""
from sfepy.base.base import output, IndexedStruct
from sfepy.discrete import FieldVariable, Integral, Equation, Equations, Material
from sfepy.discrete import Problem
from sfepy.terms import Term
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton

def create_mass_matrix(field):
    """
    Create scalar mass matrix corresponding to the given field.

    Returns
    -------
    mtx : csr_matrix
        The mass matrix in CSR format.
    """
    u = FieldVariable('u', 'unknown', field)
    v = FieldVariable('v', 'test', field, primary_var_name='u')

    integral = Integral('i', order=field.approx_order * 2)
    term = Term.new('dw_volume_dot(v, u)', integral, field.region, v=v, u=u)
    eq = Equation('aux', term)
    eqs = Equations([eq])
    eqs.time_update(None)

    dummy = eqs.create_state_vector()

    mtx = eqs.create_matrix_graph()
    mtx = eqs.eval_tangent_matrices(dummy, mtx)

    return mtx

def make_l2_projection(target, source):
    """
    Project a scalar `source` field variable to a scalar `target` field
    variable using the :math:`L^2` dot product.
    """
    def eval_variable(ts, coors, mode, **kwargs):
        val = source.evaluate_at(coors)
        val.shape = val.shape + (1,)
        return val

    make_l2_projection_data(target, eval_variable)

def make_l2_projection_data(target, eval_data):
    """
    Project scalar data given by a material-like `eval_data()` function to a
    scalar `target` field variable using the :math:`L^2` dot product.
    """
    order = target.field.approx_order * 2
    integral = Integral('i', order=order)

    un = target.name
    v = FieldVariable('v', 'test', target.field, primary_var_name=un)
    lhs = Term.new('dw_volume_dot(v, %s)' % un, integral,
                   target.field.region, v=v, **{un : target})

    def _eval_data(ts, coors, mode, **kwargs):
        if mode == 'qp':
            val = eval_data(ts, coors, mode, **kwargs)
            return {'val' : val}

    m = Material('m', function=_eval_data)
    rhs = Term.new('dw_volume_lvf(m.val, v)', integral, target.field.region,
                   m=m, v=v)

    eq = Equation('projection', lhs - rhs)
    eqs = Equations([eq])

    ls = ScipyDirect({})

    nls_status = IndexedStruct()
    nls = Newton({}, lin_solver=ls, status=nls_status)

    pb = Problem('aux', equations=eqs, nls=nls, ls=ls)

    pb.time_update()

    # This sets the target variable with the projection solution.
    pb.solve()

    if nls_status.condition != 0:
        output('L2 projection: solver did not converge!')

def make_h1_projection_data(target, eval_data):
    """
    Project scalar data given by a material-like `eval_data()` function to a
    scalar `target` field variable using the :math:`H^1` dot product.
    """
    order = target.field.approx_order * 2
    integral = Integral('i', order=order)

    un = target.name
    v = FieldVariable('v', 'test', target.field, primary_var_name=un)
    lhs1 = Term.new('dw_volume_dot(v, %s)' % un, integral,
                    target.field.region, v=v, **{un : target})
    lhs2 = Term.new('dw_laplace(v, %s)' % un, integral,
                    target.field.region, v=v, **{un : target})

    def _eval_data(ts, coors, mode, **kwargs):
        if mode == 'qp':
            val = eval_data(ts, coors, mode, 'val', **kwargs)
            gval = eval_data(ts, coors, mode, 'grad', **kwargs)
            return {'val' : val, 'gval' : gval}

    m = Material('m', function=_eval_data)
    rhs1 = Term.new('dw_volume_lvf(m.val, v)', integral, target.field.region,
                    m=m, v=v)
    rhs2 = Term.new('dw_diffusion_r(m.gval, v)', integral, target.field.region,
                    m=m, v=v)

    eq = Equation('projection', lhs1 + lhs2 - rhs1 - rhs2)
    eqs = Equations([eq])

    ls = ScipyDirect({})

    nls_status = IndexedStruct()
    nls = Newton({}, lin_solver=ls, status=nls_status)

    pb = Problem('aux', equations=eqs, nls=nls, ls=ls)

    pb.time_update()

    # This sets the target variable with the projection solution.
    pb.solve()

    if nls_status.condition != 0:
        output('H1 projection: solver did not converge!')
