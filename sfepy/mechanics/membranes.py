import numpy as nm

from sfepy.linalg import dot_sequences, insert_strided_axis
from sfepy.mechanics.tensors import dim2sym

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
                 + mtx_c[..., 0, 1]**2)

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
