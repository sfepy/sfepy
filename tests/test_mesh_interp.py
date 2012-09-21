import os.path as op

import numpy as nm

from sfepy.base.conf import transform_variables
from sfepy.base.testing import TestCommon

variables = {
    'u'       : ('unknown field', 'f', 0),
    'v'       : ('test field',    'f', 'u'),
}

def in_dir(adir):
    return lambda x: op.join(adir, x)

def gen_datas(meshes):
    datas = {}

    for key, mesh in meshes.iteritems():
        bbox = mesh.get_bounding_box()
        nx = bbox[1,0] - bbox[0,0]
        centre = 0.5 * bbox.sum(axis=0)
        mesh.coors -= centre

        data = nm.sin(4.0 * nm.pi * mesh.coors[:,0:1] / nx)
        datas['scalar_' + key] = data

        data = nm.zeros_like(mesh.coors)
        data[:,0] = 0.05 * nx * nm.sin(4.0 * nm.pi * mesh.coors[:,0] / nx)
        data[:,2] = 0.05 * nx * nm.cos(4.0 * nm.pi * mesh.coors[:,0] / nx)
        datas['vector_' + key] = data

    return datas

def do_interpolation(m2, m1, data, field_name, force=False):
    """Interpolate data from m1 to m2. """
    from sfepy.fem import Domain, H1NodalVolumeField, Variables

    fields = {
        'scalar_si' : ((1,1), 'Omega', 2),
        'vector_si' : ((3,1), 'Omega', 2),
        'scalar_tp' : ((1,1), 'Omega', 1),
        'vector_tp' : ((3,1), 'Omega', 1),
    }

    d1 = Domain('d1', m1)

    omega1 = d1.create_region('Omega', 'all')

    f = fields[field_name]

    field1 = H1NodalVolumeField('f', nm.float64, f[0], d1.regions[f[1]],
                                approx_order=f[2])
    ff = {field1.name : field1}

    vv = Variables.from_conf(transform_variables(variables), ff)
    u1 = vv['u']
    u1.set_from_mesh_vertices(data)

    d2 = Domain('d2', m2)
    omega2 = d2.create_region('Omega', 'all')

    field2 = H1NodalVolumeField('f', nm.float64, f[0], d2.regions[f[1]],
                                approx_order=f[2])
    ff2 = {field2.name : field2}

    vv2 = Variables.from_conf(transform_variables(variables), ff2)
    u2 = vv2['u']

    if not force:
        # Performs interpolation, if other field differs from self.field
        # or, in particular, is defined on a different mesh.
        u2.set_from_other(u1, strategy='interpolation', close_limit=0.5)

    else:
        coors = u2.field.get_coor()
        vals = u1.evaluate_at(coors, close_limit=0.5)
        u2.data_from_any(vals)

    return u1, u2

class Test(TestCommon):

    @staticmethod
    def from_conf(conf, options):
        test = Test(conf=conf, options=options)
        return test

    def test_interpolation(self):
        from sfepy import data_dir
        from sfepy.fem import Mesh
        from sfepy.linalg import make_axis_rotation_matrix

        fname = in_dir(self.options.out_dir)

        meshes = {
            'tp' : Mesh('original mesh', data_dir + '/meshes/3d/block.mesh'),
            'si' : Mesh('original mesh', data_dir + '/meshes/3d/cylinder.mesh'),
        }

        datas = gen_datas(meshes)

        for field_name in ['scalar_si', 'vector_si', 'scalar_tp', 'vector_tp']:
            m1 = meshes[field_name[-2:]]

            for ia, angle in enumerate(nm.linspace(0.0, nm.pi, 11)):
                self.report('%s: %d. angle: %f' % (field_name, ia, angle))
                shift = [0.0, 0.0, 0.0]
                mtx = make_axis_rotation_matrix([0, 1, 0], angle)

                m2 = m1.copy('rotated mesh')
                m2.transform_coors(mtx)

                data = datas[field_name]
                u1, u2 = do_interpolation(m2, m1, data, field_name)

                if ia == 0:
                    u1.save_as_mesh(fname('test_mesh_interp_%s_u1.vtk'
                                          % field_name))

                u2.save_as_mesh(fname('test_mesh_interp_%s_u2.%03d.vtk'
                                      % (field_name, ia)))

        return True

    def test_interpolation_two_meshes(self):
        from sfepy import data_dir
        from sfepy.fem import Mesh, Domain, H1NodalVolumeField, Variables

        m1 = Mesh('source mesh', data_dir + '/meshes/3d/block.mesh')

        m2 = Mesh('target mesh', data_dir + '/meshes/3d/cube_medium_tetra.mesh')
        m2.coors *= 2.0

        bbox = m1.get_bounding_box()
        dd = bbox[1,:] - bbox[0,:]
        data = nm.sin(4.0 * nm.pi * m1.coors[:,0:1] / dd[0]) \
               * nm.cos(4.0 * nm.pi * m1.coors[:,1:2] / dd[1])

        variables1 = {
            'u'       : ('unknown field', 'scalar_tp', 0),
            'v'       : ('test field',    'scalar_tp', 'u'),
        }

        variables2 = {
            'u'       : ('unknown field', 'scalar_si', 0),
            'v'       : ('test field',    'scalar_si', 'u'),
        }

        d1 = Domain('d1', m1)
        omega1 = d1.create_region('Omega', 'all')
        field1 = H1NodalVolumeField('scalar_tp', nm.float64, (1,1), omega1,
                                    approx_order=1)
        ff1 = {field1.name : field1}

        d2 = Domain('d2', m2)
        omega2 = d2.create_region('Omega', 'all')
        field2 = H1NodalVolumeField('scalar_si', nm.float64, (1,1), omega2,
                                    approx_order=0)
        ff2 = {field2.name : field2}

        vv1 = Variables.from_conf(transform_variables(variables1), ff1)
        u1 = vv1['u']
        u1.set_from_mesh_vertices(data)

        vv2 = Variables.from_conf(transform_variables(variables2), ff2)
        u2 = vv2['u']

        # Performs interpolation, if other field differs from self.field
        # or, in particular, is defined on a different mesh.
        u2.set_from_other(u1, strategy='interpolation', close_limit=0.1)

        fname = in_dir(self.options.out_dir)
        u1.save_as_mesh(fname('test_mesh_interp_block_scalar.vtk'))
        u2.save_as_mesh(fname('test_mesh_interp_cube_scalar.vtk'))

        return True

    def test_invariance(self):
        from sfepy import data_dir
        from sfepy.fem import Mesh

        meshes = {
            'tp' : Mesh('original mesh', data_dir + '/meshes/3d/block.mesh'),
            'si' : Mesh('original mesh', data_dir + '/meshes/3d/cylinder.mesh'),
        }
        datas = gen_datas(meshes)

        ok = True
        for field_name in ['scalar_si', 'vector_si', 'scalar_tp', 'vector_tp']:
            m1 = meshes[field_name[-2:]]

            data = datas[field_name]
            u1, u2 = do_interpolation(m1, m1, data, field_name, force=True)

            self.report('max. difference:', nm.abs(u1() - u2()).max())
            _ok = nm.allclose(u1(), u2(), rtol=0.0, atol=1e-12)
            self.report('invariance for %s field: %s' % (field_name, _ok))

            ok = ok and _ok

        return ok

    def test_invariance_qp(self):
        from sfepy import data_dir
        from sfepy.fem import (Mesh, Domain, H1NodalVolumeField,
                               Variables, Integral)
        from sfepy.terms import Term
        from sfepy.fem.mappings import get_physical_qps

        mesh = Mesh('source mesh', data_dir + '/meshes/3d/block.mesh')

        bbox = mesh.get_bounding_box()
        dd = bbox[1,:] - bbox[0,:]
        data = nm.sin(4.0 * nm.pi * mesh.coors[:,0:1] / dd[0]) \
               * nm.cos(4.0 * nm.pi * mesh.coors[:,1:2] / dd[1])

        variables = {
            'u'       : ('unknown field', 'scalar_tp', 0),
            'v'       : ('test field',    'scalar_tp', 'u'),
        }

        domain = Domain('domain', mesh)
        omega = domain.create_region('Omega', 'all')
        field = H1NodalVolumeField('scalar_tp', nm.float64, 1, omega,
                                   approx_order=1)
        ff = {field.name : field}

        vv = Variables.from_conf(transform_variables(variables), ff)
        u = vv['u']
        u.set_from_mesh_vertices(data)

        integral = Integral('i', order=2)
        term = Term.new('ev_volume_integrate(u)', integral, omega, u=u)
        term.setup()
        val1, _ = term.evaluate(mode='qp')
        val1 = val1.ravel()

        qps = get_physical_qps(omega, integral)
        coors = qps.get_merged_values()

        val2 = u.evaluate_at(coors).ravel()

        self.report('max. difference:', nm.abs(val1 - val2).max())
        ok = nm.allclose(val1, val2, rtol=0.0, atol=1e-12)
        self.report('invariance in qp: %s' % ok)

        return ok
