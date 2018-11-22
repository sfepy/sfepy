import numpy as nm

from sfepy.base.base import assert_
from sfepy.terms.extmods import terms

def eval_real(vec, conn, geo, mode, shape, bf=None):
    """
    Evaluate basic derived quantities of a real variable given its DOF
    vector, connectivity and reference mapping.
    """
    n_el, n_qp, dim, n_en, n_comp = shape
    dtype = nm.float64

    if mode == 'val':
        function = terms.dq_state_in_qp

        out = nm.empty((n_el, n_qp, n_comp, 1), dtype=dtype)
        if bf is not None:
            function(out, vec, bf, conn)
        else:
            function(out, vec, geo.bf, conn)

    elif mode == 'grad':
        function = terms.dq_grad

        out = nm.empty((n_el, n_qp, dim, n_comp), dtype=dtype)
        function(out, vec, geo, conn)

    elif mode == 'div':
        assert_(n_comp == dim)
        function = terms.dq_div_vector

        out = nm.empty((n_el, n_qp, 1, 1), dtype=dtype)
        function(out, vec, geo, conn)

    elif mode == 'cauchy_strain':
        assert_(n_comp == dim)
        function = terms.dq_cauchy_strain

        sym = (dim + 1) * dim // 2
        out = nm.empty((n_el, n_qp, sym, 1), dtype=dtype)
        function(out, vec, geo, conn)
    # TODO repair this ugly hotfix!
    # in DG I dont see reason why get values of variable in qp points
    elif mode == 'dg':
        out = nm.empty((geo.n_el, 2), dtype=dtype)
        out[:, 0] = vec[:geo.n_el]
        out[:, 1] = vec[geo.n_el:]


    else:
        raise ValueError('unsupported variable evaluation mode! (%s)'
                         % mode)

    return out

def eval_complex(vec, conn, geo, mode, shape, bf=None):
    """
    Evaluate basic derived quantities of a complex variable given its DOF
    vector, connectivity and reference mapping.
    """
    n_el, n_qp, dim, n_en, n_comp = shape

    if mode == 'val':
        function = terms.dq_state_in_qp

        rout = nm.empty((n_el, n_qp, n_comp, 1), dtype=nm.float64)
        iout = nm.empty((n_el, n_qp, n_comp, 1), dtype=nm.float64)
        if bf is not None:
            function(rout, vec.real.copy(), bf, conn)
            function(iout, vec.imag.copy(), bf, conn)
        else:
            function(rout, vec.real.copy(), geo.bf, conn)
            function(iout, vec.imag.copy(), geo.bf, conn)
        out = rout + 1j * iout

    elif mode == 'grad':
        function = terms.dq_grad

        rout = nm.empty((n_el, n_qp, dim, n_comp), dtype=nm.float64)
        iout = nm.empty((n_el, n_qp, dim, n_comp), dtype=nm.float64)
        function(rout, vec.real.copy(), geo, conn)
        function(iout, vec.imag.copy(), geo, conn)
        out = rout + 1j * iout

    elif mode == 'div':
        assert_(n_comp == dim)
        function = terms.dq_div_vector

        rout = nm.empty((n_el, n_qp, 1, 1), dtype=nm.float64)
        iout = nm.empty((n_el, n_qp, 1, 1), dtype=nm.float64)
        function(rout, vec.real.copy(), geo, conn)
        function(iout, vec.imag.copy(), geo, conn)
        out = rout + 1j * iout

    elif mode == 'cauchy_strain':
        assert_(n_comp == dim)
        function = terms.dq_cauchy_strain

        sym = (dim + 1) * dim // 2
        rout = nm.empty((n_el, n_qp, sym, 1), dtype=nm.float64)
        iout = nm.empty((n_el, n_qp, sym, 1), dtype=nm.float64)
        function(rout, vec.real.copy(), geo, conn)
        function(iout, vec.imag.copy(), geo, conn)
        out = rout + 1j * iout

    else:
        raise ValueError('unsupported variable evaluation mode! (%s)'
                         % mode)

    return out
