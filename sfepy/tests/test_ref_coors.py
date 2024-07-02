import os.path as op

import numpy as nm

import sfepy
from sfepy.discrete.common import Field
import sfepy.discrete.common.global_interp as gi
import sfepy.base.testing as tst

def test_ref_coors_fem():
    from sfepy.discrete.fem import Mesh, FEDomain

    mesh = Mesh.from_file('meshes/3d/special/cross3d.mesh',
                          prefix_dir=sfepy.data_dir)
    domain = FEDomain('domain', mesh)

    omega = domain.create_region('Omega', 'all')

    field = Field.from_args('linear', nm.float64, 'scalar', omega,
                            approx_order=1)


    mcoors = field.domain.get_mesh_coors()
    conn = field.domain.get_conn()

    bbox = field.domain.get_mesh_bounding_box()
    ray = nm.linspace(bbox[0, 0], bbox[1, 0], 7)
    coors = nm.zeros((ray.shape[0], 3), dtype=nm.float64)

    def gen_rays():
        coors[:, 0] = ray
        yield coors

        coors.fill(0.0)
        coors[:, 1] = ray
        yield coors

        coors.fill(0.0)
        coors[:, 2] = ray
        yield coors

    ok = True

    ctx = field.create_basis_context()._geo_ctx

    for ir, coors in enumerate(gen_rays()):
        tst.report('ray %d' % ir)
        ref_coors, cells, status = gi.get_ref_coors(field, coors,
                                                    strategy='general',
                                                    close_limit=0.0,
                                                    verbose=False)
        tst.report(ref_coors)
        tst.report(cells)
        tst.report(status)

        # In the distorted cell 2, the Newton method finds a solution
        # outside of the cell. This will be fixed when box constraints
        # are applied.
        _ok = nm.all((status == 0) | ((cells == 2) & (status == 3)))
        if not _ok:
            tst.report('wrong status %s for ray %d!' % (status, ir))

        ok = ok and _ok

        for ic, cell in enumerate(cells):
            ctx.iel = cell
            bf = ctx.evaluate(ref_coors[ic:ic+1], check_errors=False)

            cell_coors = mcoors[conn[cell]]
            coor = nm.dot(bf, cell_coors).ravel()

            _ok = nm.allclose(coor, coors[ic], atol=1e-14, rtol=0.0)
            if not _ok:
                tst.report('ray %d point %d:' % (ir, ic))
                tst.report(' - wrong reference coordinates %s!'
                           % ref_coors[ic])
                tst.report(' - given point: %s' % coors[ic])
                tst.report(' - found point: %s' % coor)
            ok = ok and _ok

    assert ok

def test_ref_coors_iga():
    from sfepy.discrete.iga.domain import IGDomain

    domain = IGDomain.from_file(op.join(sfepy.data_dir,
                                        'meshes/iga/block2d.iga'))

    omega = domain.create_region('Omega', 'all')

    field = Field.from_args('iga', nm.float64, 'scalar', omega,
                            approx_order='iga', poly_space_basis='iga')

    mcoors = field.nurbs.cps
    conn = field.get_econn('cell', field.region)

    bbox = domain.eval_mesh.get_bounding_box()
    ray = nm.linspace(bbox[0, 0], bbox[1, 0], 11)
    coors = nm.c_[ray, ray]

    ref_coors, cells, status = gi.get_ref_coors(field, coors,
                                                strategy='general',
                                                close_limit=0.0,
                                                verbose=False)
    tst.report(ref_coors)
    tst.report(cells)
    tst.report(status)

    ok = nm.all(status == 0)

    ctx = field.create_basis_context()

    for ic, cell in enumerate(cells):
        ctx.iel = cell
        bf = ctx.evaluate(ref_coors[ic:ic+1])

        cell_coors = mcoors[conn[cell]]
        coor = nm.dot(bf, cell_coors).ravel()

        _ok = nm.allclose(coor, coors[ic], atol=1e-14, rtol=0.0)
        if not _ok:
            tst.report('point %d:' % ic)
            tst.report(' - wrong reference coordinates %s!'
                       % ref_coors[ic])
            tst.report(' - given point: %s' % coors[ic])
            tst.report(' - found point: %s' % coor)
        ok = ok and _ok

    assert ok
