"""
Tests of the contact support based on the contactino package: its binding
(frictional traction, gap operator, configuration correction, nodal fields,
frozen state) and the contact terms with quadratic facets.
"""
import numpy as nm
import pytest

import sfepy.base.testing as tst

# contactino is an optional dependency.
pytest.importorskip('contactino')

# Two conforming 2 x 2 quads: the lower one at z = -0.1 (normal +z, nodes
# 0-3), the upper one at z = 0 (normal -z, nodes 4-7).
quad_pair_coors = nm.array(
    [-1, -1, -0.1,  1, -1, -0.1,  1, 1, -0.1,  -1, 1, -0.1,
     -1, -1,  0.0,  1, -1,  0.0,  1, 1,  0.0,  -1, 1,  0.0],
    dtype=nm.float64)
quad_pair_conn = nm.array([0, 1, 2, 3,  4, 7, 6, 5], dtype=nm.int32)

def init_quad_pair(epsn=1e4, epst=0.0, mu=0.0, algorithm=None):
    import contactino as cc

    iparm = nm.zeros(64, dtype=nm.int32)
    iparm[cc.IPARM_NSD] = 3
    iparm[cc.IPARM_NNP] = 8
    iparm[cc.IPARM_NCE] = 2
    iparm[cc.IPARM_NORMAL_ENFORCEMENT] = \
        cc.NORMAL_CONTACT_ENFORCEMENT_METHOD_PENALTY
    if mu > 0.0:
        iparm[cc.IPARM_CONTACT_TYPE] = cc.CONTACT_TYPE_FRICTIONAL
        iparm[cc.IPARM_TANGENT_ENFORCEMENT] = \
            cc.TANGENT_CONTACT_ENFORCEMENT_METHOD_PENALTY

    else:
        iparm[cc.IPARM_CONTACT_TYPE] = cc.CONTACT_TYPE_FRICTIONLESS
        iparm[cc.IPARM_TANGENT_ENFORCEMENT] = \
            cc.TANGENT_CONTACT_ENFORCEMENT_METHOD_NONE
    iparm[cc.IPARM_NQP_1D] = 2
    if algorithm is not None:
        iparm[cc.IPARM_ALGORITHM] = algorithm

    rparm = nm.zeros(64, dtype=nm.float64)
    rparm[cc.RPARM_EPS_N] = epsn
    if mu > 0.0:
        rparm[cc.RPARM_EPS_T] = epst
        rparm[cc.RPARM_MU] = mu

    types = nm.full(2, cc.CONTACTINO_QUAD, dtype=nm.int32)
    return cc.ContactinoData(iparm, rparm, quad_pair_coors, quad_pair_conn,
                             types)

def quad_pair_displacements(press, slide=0.0):
    """
    The upper quad pressed down by `press` (penetration `press` - 0.1) and
    shifted along x by `slide`.
    """
    uc = nm.zeros(24, dtype=nm.float64)
    uc[12::3] = slide
    uc[14::3] = -press
    return uc

def test_frictional_quad_pair():
    """
    The frictional contact protocol: press, commit, slide. A small slide
    sticks (traction epsT * slip), a large one slips at the Coulomb limit.
    """
    epsn, epst, mu = 1e4, 1e3, 0.3
    cd = init_quad_pair(epsn, epst, mu)

    rc = cd.eval_residual(quad_pair_displacements(0.101))
    fn = -rc[14::3].sum()
    tst.report('normal force %.6e' % fn)
    assert abs(fn - epsn * 1e-3 * 4.0) < 1e-9

    cd.commit_step()
    rc = cd.eval_residual(quad_pair_displacements(0.101, 1e-3))
    ft = -rc[12::3].sum()
    tst.report('stick: tangential force %.6e' % ft)
    assert abs(ft + epst * 1e-3 * 4.0) < 1e-9

    rc, rows, cols, vals = cd.eval_residual_and_tangent(
        quad_pair_displacements(0.101, 1e-2))
    ft = -rc[12::3].sum()
    tst.report('slip: tangential force %.6e, Kuu nnz %d' % (ft, len(vals)))
    assert abs(ft + mu * fn) < 1e-9

def test_gap_operator():
    """
    The gap operator of the conforming pair is symmetric, positive
    semidefinite, differs from the tangent divided by the penalty by the
    geometric term proportional to the penetration, and gives the change of
    the residual for a rigid shift.
    """
    import scipy.sparse as sps

    def get_matrix(rows, cols, vals):
        return sps.coo_matrix((vals, (rows, cols)), shape=(24, 24)).toarray()

    epsn = 1e4
    cd = init_quad_pair(epsn)
    uc = quad_pair_displacements(0.101)

    mtx_g = get_matrix(*cd.eval_gap_operator(uc))
    rc, rows, cols, vals = cd.eval_residual_and_tangent(uc)
    mtx_k = get_matrix(rows, cols, vals)

    asym = nm.abs(mtx_g - mtx_g.T).max() / nm.abs(mtx_g).max()
    tst.report('asymmetry %.3e' % asym)
    assert asym < 1e-13

    eigs = nm.linalg.eigvalsh(0.5 * (mtx_g + mtx_g.T))
    assert eigs.min() > -1e-13 * eigs.max()

    # The geometric term p_c Dn of the tangent is of the order of the
    # penetration over the facet size.
    rel1 = nm.abs(mtx_g - mtx_k / epsn).max() / nm.abs(mtx_g).max()
    uc2 = quad_pair_displacements(0.1005)
    rc2, rows, cols, vals = cd.eval_residual_and_tangent(uc2)
    rel2 = nm.abs(mtx_g - get_matrix(rows, cols, vals) / epsn).max() \
        / nm.abs(mtx_g).max()
    tst.report('|G - K / epsN| / |G| = %.3e (penetration 1e-3),'
               ' %.3e (5e-4)' % (rel1, rel2))
    assert rel1 < 1e-2
    assert abs(rel2 / rel1 - 0.5) < 1e-6

    # The residual change for a rigid shift (constant normal) is epsN G du.
    assert nm.allclose(rc - rc2, epsn * mtx_g @ (uc - uc2), rtol=0.0,
                       atol=1e-9)

def test_configuration_correction():
    """
    The correction of a predicted configuration with equal lumped masses
    splits the penetration in half, conserves momentum and leaves no contact
    force.
    """
    cd = init_quad_pair()
    uc = quad_pair_displacements(0.101)
    dofs = nm.arange(24, dtype=nm.int32)
    duc = cd.correct_configuration(uc, dofs, dofs, nm.ones(24))

    assert nm.allclose(duc[2::3][:4], -0.5e-3, rtol=0.0, atol=1e-15)
    assert nm.allclose(duc[2::3][4:], 0.5e-3, rtol=0.0, atol=1e-15)
    assert nm.allclose(duc[0::3], 0.0, rtol=0.0, atol=1e-15)
    assert abs(duc[2::3].sum()) < 1e-15

    rc = cd.eval_residual(uc + duc)
    tst.report('remaining contact force %.3e' % nm.abs(rc).max())
    assert nm.abs(rc).max() < 1e-9

    # A rigid master (no inverse mass) leaves the whole correction to the
    # slave.
    dofs = nm.arange(12, 24, dtype=nm.int32)
    duc = cd.correct_configuration(uc, dofs, dofs, nm.ones(12))
    assert nm.allclose(duc[2::3][:4], 0.0, rtol=0.0, atol=1e-15)
    assert nm.allclose(duc[2::3][4:], 1e-3, rtol=0.0, atol=1e-15)

def test_nodal_fields():
    """
    The nodal pressure of a flat conforming pair under a rigid shift equals
    the penalty times the penetration, the nodal gap the penetration.
    """
    epsn = 1e4
    cd = init_quad_pair(epsn)
    gap, pressure = cd.eval_nodal_fields(quad_pair_displacements(0.101))
    assert nm.allclose(gap, 1e-3, rtol=1e-12, atol=0.0)
    assert nm.allclose(pressure, epsn * 1e-3, rtol=1e-12, atol=0.0)

    gap, pressure = cd.eval_nodal_fields(quad_pair_displacements(0.099))
    assert nm.allclose(gap, -1e-3, rtol=1e-12, atol=0.0)
    assert nm.allclose(pressure, 0.0, rtol=0.0, atol=1e-15)

def test_one_pass_quad_pair():
    """
    The one-pass algorithm needs the element roles; the slave pass alone
    gives the full contact force on both bodies with exact action and
    reaction.
    """
    import contactino as cc

    epsn = 1e4
    uc = quad_pair_displacements(0.101)
    cd = init_quad_pair(epsn, algorithm=cc.CONTACT_ALGORITHM_ONE_PASS)
    with pytest.raises(RuntimeError) as exc:
        cd.eval_residual(uc)
    tst.report('without roles: %s' % exc.value)
    assert 'role' in str(exc.value).lower()

    cd = init_quad_pair(epsn, algorithm=cc.CONTACT_ALGORITHM_ONE_PASS)
    cd.set_element_roles([cc.CONTACT_ELEMENT_ROLE_MASTER,
                          cc.CONTACT_ELEMENT_ROLE_SLAVE])
    rc = cd.eval_residual(uc)
    f_lower, f_upper = -rc[2:12:3].sum(), -rc[14::3].sum()
    tst.report('one-pass forces: lower %.6e, upper %.6e' % (f_lower, f_upper))
    assert abs(f_upper - epsn * 1e-3 * 4.0) < 1e-9
    assert abs(f_lower + f_upper) < 1e-9

def test_frozen_state():
    """
    A frozen active point stays active when the contact opens (tension), the
    release counts the changed points and the free evaluation is restored.
    """
    epsn = 1e4
    cd = init_quad_pair(epsn)
    closed = quad_pair_displacements(0.101)
    opened = quad_pair_displacements(0.099)

    rc = cd.eval_residual(closed)
    cd.freeze_state()
    rc = cd.eval_residual(opened)
    f_upper = -rc[14::3].sum()
    tst.report('frozen, opened by 1e-3: force on the upper quad %.6e'
               % f_upper)
    assert abs(f_upper + epsn * 1e-3 * 4.0) < 1e-9

    n_changed = cd.release_state(opened)
    tst.report('changed points after the release: %d' % n_changed)
    assert n_changed > 0
    assert nm.abs(cd.eval_residual(opened)).max() == 0.0

    rc = cd.eval_residual(closed)
    cd.freeze_state()
    assert cd.release_state(closed) == 0

def test_error_message():
    """
    The reason of a failed call is in the exception message.
    """
    cd = init_quad_pair()
    uc = quad_pair_displacements(0.101)
    empty = nm.zeros(0, dtype=nm.int32)
    with pytest.raises(RuntimeError) as exc:
        cd.correct_configuration(uc, empty, empty, nm.zeros(0))
    tst.report('immovable constraint: %s' % exc.value)
    assert 'movable' in str(exc.value)

@pytest.mark.parametrize('dim, order, basis', [
    (2, 1, 'lagrange'), (2, 2, 'lagrange'),
    (3, 1, 'lagrange'), (3, 2, 'serendipity'), (3, 2, 'lagrange'),
])
def test_quadratic_facets(dim, order, basis):
    """
    A uniform penetration of two non-conforming blocks gives the contact
    force :math:`\\varepsilon_N \\delta A` on each body for the linear and
    the quadratic contact facets (the facet nodes of the second order fields
    are passed to contactino in its node order).
    """
    from sfepy.discrete.fem import FEDomain, Field
    from sfepy.discrete import FieldVariable, Material, Integral
    from sfepy.terms import Term
    from sfepy.mesh.mesh_generators import gen_block_mesh

    epsn, gap, delta = 1e4, 0.1, 1e-3
    dims = [2.0] * (dim - 1) + [1.0]
    centre = [1.0] * (dim - 1) + [0.5]
    lower = gen_block_mesh(dims, [4] * (dim - 1) + [2], centre,
                           name='lower', verbose=False)
    centre[-1] += 1.0 + gap
    upper = gen_block_mesh(dims, [3] * (dim - 1) + [2], centre,
                           name='upper', verbose=False)
    domain = FEDomain('domain', lower + upper)
    omega = domain.create_region('Omega', 'all')
    ax = 'xyz'[dim - 1]
    gamma = domain.create_region('Gamma',
                                 'vertices in (%s > 0.99) & (%s < %f)'
                                 % (ax, ax, 1.01 + gap), 'facet')

    field = Field.from_args('fu', nm.float64, 'vector', omega,
                            approx_order=order, poly_space_basis=basis)
    u = FieldVariable('u', 'unknown', field)
    v = FieldVariable('v', 'test', field, primary_var_name='u')
    m = Material('m', values={'.algorithm' : 'two_pass', '.epss' : epsn})
    term = Term.new('dw_contact(m.algorithm, m.epss, v, u)',
                    Integral('i', order=2), gamma, m=m, v=v, u=u)
    term.setup()

    is_upper = field.coors[:, -1] > 1.0 + 0.5 * gap
    vec = nm.zeros((field.n_nod, dim))
    vec[is_upper, -1] = -(gap + delta)
    u.set_data(vec.ravel())

    term.evaluate(mode='weak')
    nodes, fc = term.get_nodal_forces()
    f_upper = fc[is_upper[nodes], -1].sum()
    f_lower = fc[~is_upper[nodes], -1].sum()
    f_exact = epsn * delta * 2.0**(dim - 1)
    tst.report('%dD, order %d (%s), %d facet nodes: forces upper %.8e,'
               ' lower %.8e, exact %.8e'
               % (dim, order, basis, term.ci.sd.n_fp, f_upper, f_lower,
                  f_exact))
    assert nm.isclose(f_upper, f_exact, rtol=1e-10)
    assert nm.isclose(f_lower, -f_exact, rtol=1e-10)
