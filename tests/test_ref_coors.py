import numpy as nm

import sfepy
from sfepy.discrete.fem import Mesh, FEDomain, Field
import sfepy.discrete.fem.global_interp as gi

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        mesh = Mesh.from_file('meshes/3d/special/cross3d.mesh',
                              prefix_dir=sfepy.data_dir)
        domain = FEDomain('domain', mesh)

        omega = domain.create_region('Omega', 'all')

        field = Field.from_args('linear', nm.float64, 'scalar', omega,
                                approx_order=1)

        test = Test(conf=conf, options=options, omega=omega, field=field)
        return test

    def test_ref_coors(self):
        field = self.field

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

        for ir, coors in enumerate(gen_rays()):
            self.report('ray %d' % ir)
            ref_coors, cells, status = gi.get_ref_coors(field, coors,
                                                        strategy='general',
                                                        close_limit=0.0,
                                                        verbose=False)
            self.report(ref_coors)
            self.report(cells)
            self.report(status)

            # In the distorted cell 2, the Newton method founds a solution
            # outside of the cell. This will be fixed when box constraints
            # are applied.
            _ok = nm.all((status == 0) | ((cells == 2) & (status == 3)))
            if not _ok:
                self.report('wrong status %s for ray %d!' % (status, ir))

            ok = ok and _ok

            for ic, cell in enumerate(cells):
                ps = field.ap.get_poly_space('v')
                bf = ps.eval_base(ref_coors[ic:ic+1], suppress_errors=True)

                cell_coors = mcoors[conn[cell]]
                coor = nm.dot(bf, cell_coors).ravel()

                _ok = nm.allclose(coor, coors[ic], atol=1e-14, rtol=0.0)
                if not _ok:
                    self.report('ray %d point %d:' % (ir, ic))
                    self.report(' - wrong reference coordinates %s!'
                                % ref_coors[ic])
                    self.report(' - given point: %s' % coors[ic])
                    self.report(' - found point: %s' % coor)
                ok = ok and _ok

        return ok
