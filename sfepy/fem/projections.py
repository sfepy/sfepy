"""
Construct projections between FE spaces.
"""
from sfepy.base.base import output, IndexedStruct
from sfepy.fem import FieldVariable, Integral, Equation, Equations, Material
from sfepy.fem import ProblemDefinition
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
    u = FieldVariable('u', 'unknown', field, 1)
    v = FieldVariable('v', 'test', field, 1, primary_var_name='u')

    integral = Integral('i', order=field.approx_order**2)
    term = Term.new('dw_mass_scalar(v, u)', integral, field.region, v=v, u=u)
    eq = Equation('aux', term)
    eqs = Equations([eq])
    eqs.time_update(None)

    dummy = eqs.create_state_vector()

    mtx = eqs.create_matrix_graph()
    mtx = eqs.eval_tangent_matrices(dummy, mtx)

    return mtx

def make_l2_projection(target, source):
    """
    Project `source` field variable to `target` field variable using
    :math:`L^2` dot product.
    """
    order = target.field.get_true_order()**2
    integral = Integral('i', order=order)

    un = target.name
    v = FieldVariable('v', 'test', target.field, 1, primary_var_name=un)
    lhs = Term.new('dw_mass_scalar(v, %s)' % un, integral,
                   target.field.region, v=v, **{un : target})

    def eval_variable(ts, coors, mode, **kwargs):
        if mode == 'qp':
            val = source.evaluate_at(coors)
            val.shape = val.shape + (1,)
            out = {'val' : val}
            return out

    m = Material('m', function=eval_variable)
    rhs = Term.new('dw_volume_lvf(m.val, v)', integral, target.field.region,
                   m=m, v=v)

    eq = Equation('projection', lhs - rhs)
    eqs = Equations([eq])

    ls = ScipyDirect({})

    nls_status = IndexedStruct()
    nls = Newton({}, lin_solver=ls, status=nls_status)

    pb = ProblemDefinition('aux', equations=eqs, nls=nls, ls=ls)

    pb.time_update()

    # This sets the target variable with the projection solution.
    pb.solve()

    if nls_status.condition != 0:
        output('L2 projection: solver did not converge!')
