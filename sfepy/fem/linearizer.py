import numpy as nm

from sfepy.base.base import Struct
from sfepy.fem.geometry_element import GeometryElement
from sfepy.fem.poly_spaces import PolySpace
from sfepy.fem import Mesh, Domain, Field
from sfepy.fem.refine import refine_reference

def create_output(dofs, dof_coors, dof_conn, ps, max_level=3, eps=1e-4):
    """
    Create mesh with linear elements that approximates `dofs`
    corresponding to a higher order approximation with a relative
    precision given by `eps`.
    """

    def _get_msd(iels, rx, ree):
        edofs = dofs[dof_conn[iels]]
        n_components = edofs.shape[-1]

        bf = ps.eval_base(rx).squeeze()

        msd = 0.0
        for ic in range(n_components):
            edof = edofs[..., ic]
            rvals = nm.dot(edof, bf.T)

            sd = rvals[:, ree]
            # ~ max. second derivative.
            msd += nm.abs(sd[..., 0] + sd[..., 2]
                          - 2.0 * sd[..., 1]).max(axis=-1)

        msd /= n_components

        return bf, msd

    eps_r = (dofs.max() - dofs.min()) * eps

    n_el = dof_conn.shape[0]

    bf0 = ps.eval_base(ps.geometry.coors).squeeze()

    rc0 = ps.geometry.conn[None, :]
    rx, rc, ree = refine_reference(ps.geometry, 1)

    iels = nm.arange(n_el)
    bf, msd = _get_msd(iels, rx, ree)
    flag = msd > eps_r

    coors = []
    conns = []
    vdofs = []

    inod = 0

    for level in range(max_level):
        if level == (max_level - 1):
            # Last level - take everything.
            flag.fill(False)

        # Deal with finished elements.
        ie, ir = nm.where(flag == False)
        if len(ie):
            uie, iies = nm.unique(ie, return_inverse=True)

            # Each (sub-)element has own coordinates - no shared vertices.
            ecoors = dof_coors[dof_conn[iels[uie]]]
            aux = ecoors.transpose((0, 2, 1))
            xes = nm.dot(aux, bf0.T).transpose((0, 2, 1))

            edofs = dofs[dof_conn[iels[uie]]]
            aux = edofs.transpose((0, 2, 1))
            des = nm.dot(aux, bf0.T).transpose((0, 2, 1))

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
            conn = nm.arange(nc, dtype=nm.int32).reshape((nc / np, np))

            coors.append(cc)
            conns.append(conn + inod)
            vdofs.append(vd)

            inod += nc

        if not flag.any():
            break

        # Deal with elements to refine.
        if level < (max_level - 1):
            eflag = flag.sum(axis=1, dtype=nm.bool)
            iels = iels[eflag]

            rc0 = rc
            rx, rc, ree = refine_reference(ps.geometry, level + 2)

            bf0 = bf
            bf, msd = _get_msd(iels, rx, ree)

            flag = msd > eps_r

    all_coors = nm.concatenate(coors, axis=0)
    conn = nm.concatenate(conns, axis=0)
    all_vdofs = nm.concatenate(vdofs, axis=0)

    mat_ids = nm.zeros(conn.shape[0], dtype=nm.int32)

    return level, all_coors, conn, all_vdofs, mat_ids

if __name__ == '__main__':

    mesh = Mesh.from_file('meshes/elements/2_3_1.mesh')
    domain = Domain('', mesh)
    domain = domain.refine()
    domain = domain.refine()

    domain.mesh.write('linearizer0.mesh')

    omega = domain.create_region('Omega', 'all')

    field = Field('fu', nm.float64, 'scalar', omega,
                  space='H1', poly_space_base='lagrange', approx_order=3)

    dof_conns = {}
    field.setup_dof_conns(dof_conns, 1, 'volume', omega)

    gel = GeometryElement('2_3')
    ps = PolySpace.any_from_args(None, gel, 3)

    dofs = nm.arange(field.n_nod, dtype=nm.float64)
    dofs = nm.c_[dofs, dofs[::-1]]
    ## dofs = nm.arange(field.n_nod, dtype=nm.float64)
    ## dofs = dofs.reshape((field.n_nod, 1))
    dof_conn = dof_conns.values()[0]

    dof_coors = field.get_coor()

    aux = create_output(dofs, dof_coors, dof_conn, ps, 5, 1e-1)
    level, all_coors, conn, all_vdofs, mat_ids = aux

    mm = mesh.from_data('mm', all_coors, None, [conn], [mat_ids], ['2_3'])

    out = {
        'aa' : Struct(name='output_data', mode='vertex', data=all_vdofs,
                      var_name='aa', dofs=None)
    }

    mm.write('linearizer.mesh')
    mm.write('linearizer.vtk', out=out)

    print level
