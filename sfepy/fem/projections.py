"""
Construct projections between FE spaces.
"""
from sfepy.fem import FieldVariable, Integral, Equation, Equations
from sfepy.terms import Term

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
