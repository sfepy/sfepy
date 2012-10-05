import os
import numpy as nm

from sfepy.base.testing import TestCommon

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        test.join = lambda x: os.path.join(test.options.out_dir, x)
        return test

    def test_linearization(self):
        from sfepy.base.base import Struct
        from sfepy.fem import Mesh, Domain, H1NodalVolumeField
        from sfepy import data_dir

        geometries = ['2_3', '2_4', '3_4', '3_8']
        approx_orders = [1, 2]
        funs = [nm.cos, nm.sin, lambda x: x]

        ok = True
        for geometry in geometries:
            name = os.path.join(data_dir,
                                'meshes/elements/%s_1.mesh' % geometry)
            mesh = Mesh.from_file(name)

            domain = Domain('', mesh)
            domain = domain.refine()

            domain.mesh.write(self.join('linearizer-%s-0.mesh' % geometry))

            omega = domain.create_region('Omega', 'all')

            for approx_order in approx_orders:
                for dpn in [1, mesh.dim]:
                    self.report('geometry: %s, approx. order: %d, dpn: %d' %
                                (geometry, approx_order, dpn))

                    field = H1NodalVolumeField('fu', nm.float64, dpn, omega,
                                               approx_order=approx_order)

                    cc = field.get_coor()
                    dofs = nm.zeros((field.n_nod, dpn), dtype=nm.float64)

                    for ic in range(dpn):
                        dofs[:, ic] = funs[ic](3 * (cc[:, 0] * cc[:, 1]))

                    vmesh, vdofs, levels = field.linearize(dofs,
                                                           min_level=0,
                                                           max_level=3,
                                                           eps=1e-2)
                    level = levels[0]

                    if approx_order == 1:
                        _ok = level == 0

                    else:
                        _ok = level > 0
                    self.report('max. refinement level: %d: %s' % (level, _ok))

                    ok = ok and _ok

                    rdofs = nm.zeros((vmesh.n_nod, dpn), dtype=nm.float64)
                    cc = vmesh.coors
                    for ic in range(dpn):
                        rdofs[:, ic] = funs[ic](3 * (cc[:, 0] * cc[:, 1]))

                    _ok = nm.allclose(rdofs, vdofs, rtol=0.0, atol=0.03)
                    self.report('interpolation: %s' % _ok)
                    ok = ok and _ok

                    out = {
                        'u' : Struct(name='output_data',
                                     mode='vertex', data=vdofs,
                                     var_name='u', dofs=None)
                    }

                    name = self.join('linearizer-%s-%d-%d'
                                     % (geometry, approx_order, dpn))

                    vmesh.write(name + '.mesh')
                    vmesh.write(name + '.vtk', out=out)

        return ok
