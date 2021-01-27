from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import assert_
from sfepy.linalg import norm_l2_along_axis as norm
from sfepy.linalg import dot_sequences, insert_strided_axis
from sfepy.discrete import PolySpace
from sfepy.discrete.fem.mappings import VolumeMapping
from sfepy.mechanics.tensors import dim2sym
from six.moves import range

def create_transformation_matrix(coors):
    """
    Create a transposed coordinate transformation matrix, that
    transforms 3D coordinates of element face nodes so that the
    transformed nodes are in the `x-y` plane. The rotation is performed
    w.r.t. the first node of each face.

    Parameters
    ----------
    coors : array
        The coordinates of element nodes, shape `(n_el, n_ep, dim)`.

    Returns
    -------
    mtx_t : array
        The transposed transformation matrix :math:`T`, i.e.
        :math:`X_{inplane} = T^T X_{3D}`.

    Notes
    -----
    :math:`T = [t_1, t_2, n]`, where :math:`t_1`, :math:`t_2`, are unit
    in-plane (column) vectors and :math:`n` is the unit normal vector,
    all mutually orthonormal.
    """
    # Local coordinate system.
    t1 = coors[:, 1, :] - coors[:, 0, :]
    t2 = coors[:, -1, :] - coors[:, 0, :]
    n = nm.cross(t1, t2)
    t2 = nm.cross(n, t1)

    t1 = t1 / norm(t1)[:, None]
    t2 = t2 / norm(t2)[:, None]
    n = n / norm(n)[:, None]

    # Coordinate transformation matrix (transposed!).
    mtx_t = nm.concatenate((t1[:, :, None],
                            t2[:, :, None],
                            n[:, :, None]), axis=2)

    return mtx_t

def transform_asm_vectors(out, mtx_t):
    """
    Transform vector assembling contributions to global coordinate system, one
    node at a time.

    Parameters
    ----------
    out : array
        The array of vectors, transformed in-place.
    mtx_t : array
        The transposed transformation matrix :math:`T`, see
        :func:`create_transformation_matrix`.
    """
    n_ep = out.shape[2] // mtx_t.shape[2]
    for iep in range(n_ep):
        ir = slice(iep, None, n_ep)
        fn = out[:, 0, ir, 0]
        fn[:] = dot_sequences(mtx_t, fn, 'AB')

def transform_asm_matrices(out, mtx_t):
    """
    Transform matrix assembling contributions to global coordinate system, one
    node at a time.

    Parameters
    ----------
    out : array
        The array of matrices, transformed in-place.
    mtx_t : array
        The transposed transformation matrix :math:`T`, see
        :func:`create_transformation_matrix`.
    """
    n_ep = out.shape[-1] // mtx_t.shape[-1]
    for iepr in range(n_ep):
        ir = slice(iepr, None, n_ep)
        for iepc in range(n_ep):
            ic = slice(iepc, None, n_ep)
            fn = out[:, 0, ir, ic]
            fn[:] = dot_sequences(dot_sequences(mtx_t, fn, 'AB'), mtx_t, 'ABT')

def create_mapping(coors, gel, order):
    """
    Create mapping from transformed (in `x-y` plane) element faces to
    reference element faces.

    Parameters
    ----------
    coors : array
        The transformed coordinates of element nodes, shape `(n_el,
        n_ep, dim)`. The function verifies that the all `z` components
        are zero.
    gel : GeometryElement instance
        The geometry element corresponding to the faces.
    order : int
        The polynomial order of the mapping.

    Returns
    -------
    mapping : VolumeMapping instance
        The reference element face mapping.
    """
    # Strip 'z' component (should be 0 now...).
    assert_(nm.allclose(coors[:, :, -1], 0.0, rtol=1e-12, atol=1e-12))
    coors = coors[:, :, :-1].copy()

    # Mapping from transformed element to reference element.
    sh = coors.shape
    seq_coors = coors.reshape((sh[0] * sh[1], sh[2]))
    seq_conn = nm.arange(seq_coors.shape[0], dtype=nm.int32)
    seq_conn.shape = sh[:2]

    mapping = VolumeMapping(seq_coors, seq_conn, gel=gel, order=1)

    return mapping

def describe_geometry(field, region, integral):
    """
    Describe membrane geometry in a given region.

    Parameters
    ----------
    field : Field instance
        The field defining the FE approximation.
    region : Region instance
        The surface region to describe.
    integral : Integral instance
        The integral defining the quadrature points.

    Returns
    -------
    mtx_t : array
        The transposed transformation matrix :math:`T`, see
        :func:`create_transformation_matrix`.
    membrane_geo : CMapping instance
        The mapping from transformed elements to a reference elements.
    """
    # Coordinates of element vertices.
    sg, _ = field.get_mapping(region, integral, 'surface')
    sd = field.surface_data[region.name]
    coors = field.coors[sd.econn[:, :sg.n_ep]]

    # Coordinate transformation matrix (transposed!).
    mtx_t = create_transformation_matrix(coors)

    # Transform coordinates to the local coordinate system.
    coors_loc = dot_sequences((coors - coors[:, 0:1, :]), mtx_t)

    # Mapping from transformed elements to reference elements.
    gel = field.gel.surface_facet
    vm = create_mapping(coors_loc, gel, 1)

    qp = integral.get_qp(gel.name)
    ps = PolySpace.any_from_args(None, gel, field.approx_order)
    membrane_geo = vm.get_mapping(qp[0], qp[1], poly_space=ps)
    membrane_geo.bf[:] = ps.eval_base(qp[0])

    return mtx_t, membrane_geo

def describe_deformation(el_disps, bfg):
    """
    Describe deformation of a thin incompressible 2D membrane in 3D
    space, composed of flat finite element faces.

    The coordinate system of each element (face), i.e. the membrane
    mid-surface, should coincide with the `x`, `y` axes of the `x-y`
    plane.

    Parameters
    ----------
    el_disps : array
        The displacements of element nodes, shape `(n_el, n_ep, dim)`.
    bfg : array
        The in-plane base function gradients, shape `(n_el, n_qp, dim-1,
        n_ep)`.

    Returns
    -------
    mtx_c ; array
        The in-plane right Cauchy-Green deformation tensor
        :math:`C_{ij}`, :math:`i, j = 1, 2`.
    c33 : array
        The component :math:`C_{33}` computed from the incompressibility
        condition.
    mtx_b : array
        The discrete Green strain variation operator.
    """
    sh = bfg.shape
    n_ep = sh[3]

    dim = el_disps.shape[2]
    sym2 = dim2sym(dim-1)

    # Repeat el_disps by number of quadrature points.
    el_disps_qp = insert_strided_axis(el_disps, 1, bfg.shape[1])

    # Transformed (in-plane) displacement gradient with
    # shape (n_el, n_qp, 2 (-> a), 3 (-> i)), du_i/dX_a.
    du = dot_sequences(bfg, el_disps_qp)

    # Deformation gradient F w.r.t. in plane coordinates.
    # F_{ia} = dx_i / dX_a,
    # a \in {1, 2} (rows), i \in {1, 2, 3} (columns).
    mtx_f = du + nm.eye(dim - 1, dim, dtype=du.dtype)

    # Right Cauchy-Green deformation tensor C.
    # C_{ab} = F_{ka} F_{kb}, a, b \in {1, 2}.
    mtx_c = dot_sequences(mtx_f, mtx_f, 'ABT')

    # C_33 from incompressibility.
    c33 = 1.0 / (mtx_c[..., 0, 0] * mtx_c[..., 1, 1]
                 - mtx_c[..., 0, 1]**2)

    # Discrete Green strain variation operator.
    mtx_b = nm.empty((sh[0], sh[1], sym2, dim * n_ep), dtype=nm.float64)
    mtx_b[..., 0, 0*n_ep:1*n_ep] = bfg[..., 0, :] * mtx_f[..., 0, 0:1]
    mtx_b[..., 0, 1*n_ep:2*n_ep] = bfg[..., 0, :] * mtx_f[..., 0, 1:2]
    mtx_b[..., 0, 2*n_ep:3*n_ep] = bfg[..., 0, :] * mtx_f[..., 0, 2:3]
    mtx_b[..., 1, 0*n_ep:1*n_ep] = bfg[..., 1, :] * mtx_f[..., 1, 0:1]
    mtx_b[..., 1, 1*n_ep:2*n_ep] = bfg[..., 1, :] * mtx_f[..., 1, 1:2]
    mtx_b[..., 1, 2*n_ep:3*n_ep] = bfg[..., 1, :] * mtx_f[..., 1, 2:3]
    mtx_b[..., 2, 0*n_ep:1*n_ep] = bfg[..., 1, :] * mtx_f[..., 0, 0:1] \
                                   + bfg[..., 0, :] * mtx_f[..., 1, 0:1]
    mtx_b[..., 2, 1*n_ep:2*n_ep] = bfg[..., 0, :] * mtx_f[..., 1, 1:2] \
                                   + bfg[..., 1, :] * mtx_f[..., 0, 1:2]
    mtx_b[..., 2, 2*n_ep:3*n_ep] = bfg[..., 0, :] * mtx_f[..., 1, 2:3] \
                                   + bfg[..., 1, :] * mtx_f[..., 0, 2:3]

    return mtx_c, c33, mtx_b

def get_tangent_stress_matrix(stress, bfg):
    """
    Get the tangent stress matrix of a thin incompressible 2D membrane
    in 3D space, given a stress.

    Parameters
    ----------
    stress : array
        The components `11, 22, 12` of the second Piola-Kirchhoff stress
        tensor, shape `(n_el, n_qp, 3, 1)`.
    bfg : array
        The in-plane base function gradients, shape `(n_el, n_qp, dim-1,
        n_ep)`.

    Returns
    -------
    mtx : array
        The tangent stress matrix, shape `(n_el, n_qp, dim*n_ep, dim*n_ep)`.
    """
    n_el, n_qp, dim, n_ep = bfg.shape
    dim += 1

    mtx = nm.zeros((n_el, n_qp, dim * n_ep, dim * n_ep), dtype=nm.float64)

    g1tg1 = dot_sequences(bfg[..., 0:1, :], bfg[..., 0:1, :], 'ATB')
    g1tg2 = dot_sequences(bfg[..., 0:1, :], bfg[..., 1:2, :], 'ATB')
    g2tg1 = dot_sequences(bfg[..., 1:2, :], bfg[..., 0:1, :], 'ATB')
    g2tg2 = dot_sequences(bfg[..., 1:2, :], bfg[..., 1:2, :], 'ATB')

    aux = stress[..., 0:1, :] * g1tg1 + stress[..., 2:3, :] * g1tg2 \
          + stress[..., 2:3, :] * g2tg1 + stress[..., 1:2, :] * g2tg2

    mtx[..., 0 * n_ep : 1 * n_ep, 0 * n_ep : 1 * n_ep] = aux
    mtx[..., 1 * n_ep : 2 * n_ep, 1 * n_ep : 2 * n_ep] = aux
    mtx[..., 2 * n_ep : 3 * n_ep, 2 * n_ep : 3 * n_ep] = aux

    return mtx

def get_invariants(mtx_c, c33):
    """
    Get the first and second invariants of the right Cauchy-Green
    deformation tensor describing deformation of an incompressible
    membrane.

    Parameters
    ----------
    mtx_c ; array
        The in-plane right Cauchy-Green deformation tensor
        :math:`C_{ij}`, :math:`i, j = 1, 2`, shape `(n_el, n_qp, dim-1,
        dim-1)`.
    c33 : array
        The component :math:`C_{33}` computed from the incompressibility
        condition, shape `(n_el, n_qp)`.

    Returns
    -------
    i1 : array
        The first invariant of :math:`C_{ij}`.
    i2 : array
        The second invariant of :math:`C_{ij}`.
    """
    i1 = mtx_c[..., 0, 0] + mtx_c[..., 1, 1] + c33

    i2 = mtx_c[..., 0, 0] * mtx_c[..., 1,1] \
         + mtx_c[..., 1, 1] * c33 \
         + mtx_c[..., 0, 0] * c33 \
         - mtx_c[..., 0, 1]**2

    return i1, i2

def get_green_strain_sym3d(mtx_c, c33):
    r"""
    Get the 3D Green strain tensor in symmetric storage.

    Parameters
    ----------
    mtx_c ; array
        The in-plane right Cauchy-Green deformation tensor
        :math:`C_{ij}`, :math:`i, j = 1, 2`, shape `(n_el, n_qp, dim-1,
        dim-1)`.
    c33 : array
        The component :math:`C_{33}` computed from the incompressibility
        condition, shape `(n_el, n_qp)`.

    Returns
    -------
    mtx_e : array
        The membrane Green strain :math:`E_{ij} = \frac{1}{2} (C_{ij}) -
        \delta_{ij}`, symmetric storage: items (11, 22, 33, 12, 13, 23),
        shape `(n_el, n_qp, sym, 1)`.
    """
    n_el, n_qp, dm, _ = mtx_c.shape
    dim = dm + 1
    sym = dim2sym(dim)

    mtx_e = nm.empty((n_el, n_qp, sym, 1), dtype=mtx_c.dtype)

    mtx_e[..., 0, 0] = 0.5 * (mtx_c[..., 0, 0] - 1.0)
    mtx_e[..., 1, 0] = 0.5 * (mtx_c[..., 1, 1] - 1.0)
    mtx_e[..., 2, 0] = 0.5 * (c33 - 1.0)
    mtx_e[..., 3, 0] = 0.5 * mtx_c[..., 0, 1]
    mtx_e[..., 4:, 0] = 0.0

    return mtx_e
