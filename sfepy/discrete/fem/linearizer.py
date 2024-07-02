"""
Linearization of higher order solutions for the purposes of visualization.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.linalg import dot_sequences
from sfepy.discrete.fem.refine import refine_reference
from six.moves import range

def get_eval_dofs(dofs, dof_conn, ps, ori=None):
    """
    Get default function for evaluating field DOFs given a list of elements and
    reference element coordinates.
    """
    def _eval(iels, rx):
        edofs = dofs[dof_conn[iels]]

        if ori is not None:
            eori = ori[iels]

        else:
            eori = None

        bf = ps.eval_basis(rx, ori=eori, force_axis=True)[...,0,:]
        rvals = dot_sequences(bf, edofs)

        return rvals

    return _eval

def get_eval_coors(coors, conn, ps):
    """
    Get default function for evaluating physical coordinates given a list of
    elements and reference element coordinates.
    """
    def _eval(iels, rx):
        ecoors = coors[conn[iels]]
        aux = ecoors.transpose((0, 2, 1))

        bf = ps.eval_basis(rx).squeeze()
        phys_coors = nm.dot(aux, bf.T).transpose((0, 2, 1))
        return phys_coors

    return _eval

def create_output(eval_dofs, eval_coors, n_el, ps, min_level=0, max_level=2,
                  eps=1e-4):
    """
    Create mesh with linear elements that approximates DOFs returned by
    `eval_dofs()` corresponding to a higher order approximation with a relative
    precision given by `eps`. The DOFs are evaluated in physical coordinates
    returned by `eval_coors()`.
    """

    def _get_msd(iels, rx, ree):
        rvals = eval_dofs(iels, rx)
        rng = rvals.max() - rvals.min()
        n_components = rvals.shape[-1]

        msd = 0.0
        for ic in range(n_components):
            rval = rvals[..., ic]

            sd = rval[:, ree]
            # ~ max. second derivative.
            msd += nm.abs(sd[..., 0] + sd[..., 2]
                          - 2.0 * sd[..., 1]).max(axis=-1)

        msd /= n_components

        return msd, rng

    rx0 = ps.geometry.coors

    rc0 = ps.geometry.conn[None, :]
    rx, rc, ree = refine_reference(ps.geometry, 1)

    factor = rc.shape[0] / rc0.shape[0]

    iels = nm.arange(n_el)
    msd, rng = _get_msd(iels, rx, ree)
    eps_r = rng * eps
    flag = msd > eps_r

    iels0 = flag0 = None

    coors = []
    conns = []
    vdofs = []

    inod = 0
    for level in range(max_level + 1):
        if level < min_level:
            flag.fill(True) # Force refinement everywhere.

        elif level == max_level:
            # Last level - take everything.
            flag.fill(False)

        # Deal with finished elements.
        if flag0 is not None:
            ii = nm.searchsorted(iels0, iels)
            expand_flag0 = flag0[ii].repeat(factor, axis=1)

        else:
            expand_flag0 = nm.ones_like(flag)

        ie, ir = nm.where((flag == False) & (expand_flag0 == True))
        if len(ie):
            uie, iies = nm.unique(ie, return_inverse=True)

            # Each (sub-)element has own coordinates - no shared vertices.
            xes = eval_coors(iels[uie], rx0)
            des = eval_dofs(iels[uie], rx0)

            # Vectorize (how??) or use cython?
            cc = []
            vd = []
            for ii, iie in enumerate(iies):
                ce = rc0[ir[ii]]

                xe = xes[iie]
                cc.append(xe[ce])

                de = des[iie]
                vd.append(de[ce])

            cc = nm.vstack(cc)
            vd = nm.vstack(vd)

            nc = cc.shape[0]
            np = rc0.shape[1]
            conn = nm.arange(nc, dtype=nm.int32).reshape((nc // np, np))

            coors.append(cc)
            conns.append(conn + inod)
            vdofs.append(vd)

            inod += nc

        if not flag.any():
            break

        iels0 = iels
        flag0 = flag

        # Deal with elements to refine.
        if level < max_level:
            eflag = flag.sum(axis=1, dtype=bool)
            iels = iels[eflag]

            rc0 = rc
            rx0 = rx
            rx, rc, ree = refine_reference(ps.geometry, level + 2)

            msd, rng = _get_msd(iels, rx, ree)
            eps_r = rng * eps
            flag = msd > eps_r

    all_coors = nm.concatenate(coors, axis=0)
    conn = nm.concatenate(conns, axis=0)
    all_vdofs = nm.concatenate(vdofs, axis=0)

    mat_ids = nm.zeros(conn.shape[0], dtype=nm.int32)

    return level, all_coors, conn, all_vdofs, mat_ids
