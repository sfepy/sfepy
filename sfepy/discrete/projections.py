"""
Construct projections between FE spaces.
"""
import numpy as nm
import scipy.sparse as sps

from sfepy.base.base import output, IndexedStruct
from sfepy.discrete import (FieldVariable, Integral,
                            Equation, Equations, Material)
from sfepy.discrete import Problem
from sfepy.terms import Term
from sfepy.solvers.ls import ScipyDirect, solve
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
    term = Term.new('dw_dot(v, u)', integral, field.region, v=v, u=u)
    eq = Equation('aux', term)
    eqs = Equations([eq])
    eqs.time_update(None)
    eqs.init_state()

    dummy = eqs.create_vec()

    mtx = eqs.create_matrix_graph()
    mtx = eqs.eval_tangent_matrices(dummy, mtx)

    return mtx

def project_by_component(tensor, tensor_qp, component, order,
                         ls=None, nls_options=None):
    """
    Wrapper around make_l2_projection_data() for non-scalar fields.
    """
    aux = []
    for ic in range(tensor_qp.shape[-2]):
        make_l2_projection_data(component, tensor_qp[..., ic, :].copy(),
                                order=order, ls=ls, nls_options=nls_options)
        aux.append(component())
    tensor.set_data(nm.array(aux).T.ravel())

def make_l2_projection(target, source, ls=None, nls_options=None):
    """
    Project a scalar `source` field variable to a scalar `target` field
    variable using the :math:`L^2` dot product.
    """
    def eval_variable(ts, coors, mode, **kwargs):
        return source.evaluate_at(coors)

    make_l2_projection_data(target, eval_variable,
                            ls=ls, nls_options=nls_options)

def make_l2_projection_data(target, eval_data, order=None,
                            ls=None, nls_options=None):
    """
    Project scalar data to a scalar `target` field variable using the
    :math:`L^2` dot product.

    Parameters
    ----------
    target : FieldVariable instance
        The target variable.
    eval_data : callable or array
        Either a material-like function `eval_data()`, or an array of values in
        quadrature points that has to be reshapable to the shape required by
        `order`.
    order : int, optional
        The quadrature order. If not given, it is set to
        `2 * target.field.approx_order`.
    """
    if order is None:
       order = 2 * target.field.approx_order
    integral = Integral('i', order=order)

    un = FieldVariable('u', 'unknown', target.field)

    v = FieldVariable('v', 'test', un.field, primary_var_name=un.name)
    lhs = Term.new('dw_dot(v, %s)' % un.name, integral,
                   un.field.region, v=v, **{un.name : un})

    def _eval_data(ts, coors, mode, **kwargs):
        if mode == 'qp':
            if callable(eval_data):
                val = eval_data(ts, coors, mode, **kwargs)

            else:
                val = eval_data.reshape((coors.shape[0], 1, 1))

            return {'val' : val}

    m = Material('m', function=_eval_data)
    rhs = Term.new('dw_volume_lvf(m.val, v)', integral, un.field.region,
                   m=m, v=v)

    eq = Equation('projection', lhs - rhs)
    eqs = Equations([eq])

    if ls is None:
        ls = ScipyDirect({})

    if nls_options is None:
        nls_options = {}

    nls_status = IndexedStruct()
    nls = Newton(nls_options, lin_solver=ls, status=nls_status)

    pb = Problem('aux', equations=eqs)
    pb.set_solver(nls)

    # This sets the un variable with the projection solution.
    pb.solve(save_results=False)

    # Copy the projection solution to target.
    target.set_data(un())

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
    lhs1 = Term.new('dw_dot(v, %s)' % un, integral,
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

    pb = Problem('aux', equations=eqs)
    pb.set_solver(nls)

    # This sets the target variable with the projection solution.
    pb.solve(save_results=False)

    if nls_status.condition != 0:
        output('H1 projection: solver did not converge!')

def project_to_facets(region, fun, dpn, field):
    """
    Project a function `fun` to the `field` in facets of the given `region`.
    """
    aux = field.get_dofs_in_region(region)
    nods = nm.unique(aux)
    n_dof = len(nods)

    # Region facet connectivity.
    lconn = field.get_econn('facet', region, local=True)

    # Cell and face(cell) ids for each facet.
    fis = region.get_facet_indices()

    all_qps, all_fbfs, all_dets = field.get_surface_basis(region)

    # DOF values in the physical BQP.
    all_qps = nm.concatenate(all_qps)
    vals = nm.asarray(fun(all_qps))
    if (vals.ndim > 1) and (vals.shape != (len(all_qps), dpn)):
        raise ValueError('The projected function return value should be'
                         ' (n_point, dpn) == %s, instead of %s!'
                         % ((len(all_qps), dpn), vals.shape))
    vals.shape = (len(all_qps), dpn)

    n_qp_face = all_dets[0].shape[0]

    # Assemble l2 projection system.
    rhs = nm.zeros((dpn, n_dof), dtype=nm.float64)
    rows, cols, mvals = [], [], []
    for ii, (ie, ifa) in enumerate(fis):
        # Assembling indices.
        elc = lconn[ii]

        fvals = vals[n_qp_face * ii : n_qp_face * (ii + 1)]

        fbfs = all_fbfs[ii]
        dets = all_dets[ii]

        # Local projection system.
        for idof in range(dpn):
            lrhs = (fbfs * (fvals[:, idof, None] * dets)).sum(0)
            rhs[idof, elc] += lrhs

        lmtx = ((fbfs[..., None] * fbfs[:, None, :])
                * dets[..., None]).sum(0)

        er, ec = nm.meshgrid(elc, elc)
        rows.append(er.ravel())
        cols.append(ec.ravel())
        mvals.append(lmtx.ravel())

    rows = nm.concatenate(rows)
    cols = nm.concatenate(cols)
    mvals = nm.concatenate(mvals)
    mtx = sps.coo_matrix((mvals, (rows, cols)), shape=(n_dof, n_dof)).tocsc()

    vals = nm.zeros((n_dof, dpn), dtype=nm.float64)

    # Solve l2 projection system.
    for idof in range(dpn):
        dofs = solve(mtx, rhs[idof, :])
        vals[:, idof] = dofs

    return vals
