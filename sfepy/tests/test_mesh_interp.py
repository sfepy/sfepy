import os.path as op

import numpy as nm

from sfepy.base.conf import transform_variables
import sfepy.base.testing as tst

variables = {
    'u'       : ('unknown field', 'f', 0),
    'v'       : ('test field',    'f', 'u'),
}

def in_dir(adir):
    return lambda x: op.join(adir, x)

def gen_datas(meshes):
    datas = {}

    for key, mesh in meshes.items():
        bbox = mesh.get_bounding_box()
        nx = bbox[1,0] - bbox[0,0]
        centre = 0.5 * bbox.sum(axis=0)
        mesh.coors[:] -= centre

        data = nm.sin(4.0 * nm.pi * mesh.coors[:,0:1] / nx)
        datas['scalar_' + key] = data

        data = nm.zeros_like(mesh.coors)
        data[:,0] = 0.05 * nx * nm.sin(4.0 * nm.pi * mesh.coors[:,0] / nx)
        data[:,2] = 0.05 * nx * nm.cos(4.0 * nm.pi * mesh.coors[:,0] / nx)
        datas['vector_' + key] = data

    return datas

def do_interpolation(m2, m1, data, field_name, force=False):
    """Interpolate data from m1 to m2. """
    from sfepy.discrete import Variables
    from sfepy.discrete.fem import FEDomain, Field

    fields = {
        'scalar_si' : ((1,1), 'Omega', 2),
        'vector_si' : ((3,1), 'Omega', 2),
        'scalar_tp' : ((1,1), 'Omega', 1),
        'vector_tp' : ((3,1), 'Omega', 1),
    }

    d1 = FEDomain('d1', m1)

    d1.create_region('Omega', 'all')

    f = fields[field_name]

    field1 = Field.from_args('f', nm.float64, f[0], d1.regions[f[1]],
                             approx_order=f[2])
    ff = {field1.name : field1}

    vv = Variables.from_conf(transform_variables(variables), ff)
    u1 = vv['u']
    u1.set_from_mesh_vertices(data)

    d2 = FEDomain('d2', m2)
    d2.create_region('Omega', 'all')

    field2 = Field.from_args('f', nm.float64, f[0], d2.regions[f[1]],
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
        u2.set_data(vals)

    return u1, u2

def prepare_variable(filename, n_components):
    from sfepy.discrete import FieldVariable
    from sfepy.discrete.fem import Mesh, FEDomain, Field

    mesh = Mesh.from_file(filename)

    bbox = mesh.get_bounding_box()
    dd = bbox[1,:] - bbox[0,:]
    data = (nm.sin(4.0 * nm.pi * mesh.coors[:,0:1] / dd[0])
            * nm.cos(4.0 * nm.pi * mesh.coors[:,1:2] / dd[1]))

    domain = FEDomain('domain', mesh)
    omega = domain.create_region('Omega', 'all')
    field = Field.from_args('field', nm.float64, n_components, omega,
                            approx_order=2)

    u = FieldVariable('u', 'parameter', field,
                      primary_var_name='(set-to-None)')
    u.set_from_mesh_vertices(data * nm.arange(1, n_components + 1)[None, :])

    return u

def test_interpolation(output_dir):
    from sfepy import data_dir
    from sfepy.discrete.fem import Mesh
    from sfepy.linalg import make_axis_rotation_matrix

    fname = in_dir(output_dir)

    meshes = {
        'tp' : Mesh.from_file(data_dir + '/meshes/3d/block.mesh'),
        'si' : Mesh.from_file(data_dir + '/meshes/3d/cylinder.mesh'),
    }

    datas = gen_datas(meshes)

    for field_name in ['scalar_si', 'vector_si', 'scalar_tp', 'vector_tp']:
        m1 = meshes[field_name[-2:]]
        for ia, angle in enumerate(nm.linspace(0.0, nm.pi, 11)):
            tst.report('%s: %d. angle: %f' % (field_name, ia, angle))
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

def test_interpolation_two_meshes(output_dir):
    from sfepy import data_dir
    from sfepy.discrete import Variables
    from sfepy.discrete.fem import Mesh, FEDomain, Field

    m1 = Mesh.from_file(data_dir + '/meshes/3d/block.mesh')

    m2 = Mesh.from_file(data_dir + '/meshes/3d/cube_medium_tetra.mesh')
    m2.coors[:] *= 2.0

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

    d1 = FEDomain('d1', m1)
    omega1 = d1.create_region('Omega', 'all')
    field1 = Field.from_args('scalar_tp', nm.float64, (1,1), omega1,
                             approx_order=1)
    ff1 = {field1.name : field1}

    d2 = FEDomain('d2', m2)
    omega2 = d2.create_region('Omega', 'all')
    field2 = Field.from_args('scalar_si', nm.float64, (1,1), omega2,
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

    fname = in_dir(output_dir)
    u1.save_as_mesh(fname('test_mesh_interp_block_scalar.vtk'))
    u2.save_as_mesh(fname('test_mesh_interp_cube_scalar.vtk'))

def test_invariance():
    from sfepy import data_dir
    from sfepy.discrete.fem import Mesh

    meshes = {
        'tp' : Mesh.from_file(data_dir + '/meshes/3d/block.mesh'),
        'si' : Mesh.from_file(data_dir + '/meshes/3d/cylinder.mesh'),
    }
    datas = gen_datas(meshes)

    ok = True
    for field_name in ['scalar_si', 'vector_si', 'scalar_tp', 'vector_tp']:
        m1 = meshes[field_name[-2:]]

        data = datas[field_name]
        u1, u2 = do_interpolation(m1, m1, data, field_name, force=True)

        tst.report('max. difference:', nm.abs(u1() - u2()).max())
        _ok = nm.allclose(u1(), u2(), rtol=0.0, atol=1e-12)
        tst.report('invariance for %s field: %s' % (field_name, _ok))

        ok = ok and _ok

    assert ok

def test_invariance_qp():
    from sfepy import data_dir
    from sfepy.discrete import Integral
    from sfepy.terms import Term
    from sfepy.discrete.common.mappings import get_physical_qps

    ok = True
    for name in ['meshes/3d/block.mesh', 'meshes/3d/cylinder.mesh',
                 'meshes/2d/square_quad.mesh',
                 'meshes/2d/square_unit_tri.mesh']:
        tst.report(name)

        u = prepare_variable(op.join(data_dir, name), n_components=3)
        omega = u.field.region

        integral = Integral('i', order=3)
        qps = get_physical_qps(omega, integral)
        coors = qps.values

        term = Term.new('ev_integrate(u)', integral, omega, u=u)
        term.setup()
        val1 = term.evaluate(mode='qp')
        val1 = val1.ravel()

        val2 = u.evaluate_at(coors).ravel()

        tst.report('value: max. difference:', nm.abs(val1 - val2).max())
        ok1 = nm.allclose(val1, val2, rtol=0.0, atol=1e-12)
        tst.report('->', ok1)

        term = Term.new('ev_grad(u)', integral, omega, u=u)
        term.setup()
        val1 = term.evaluate(mode='qp')
        val1 = val1.transpose((0, 1, 3, 2)).ravel()

        val2 = u.evaluate_at(coors, mode='grad').ravel()

        tst.report('gradient: max. difference:', nm.abs(val1 - val2).max())
        ok2 = nm.allclose(val1, val2, rtol=0.0, atol=1e-10)
        tst.report('->', ok2)

        ok = ok and ok1 and ok2

    assert ok

def test_field_gradient():
    from sfepy import data_dir

    ok = True
    for name in ['meshes/3d/block.mesh', 'meshes/3d/cylinder.mesh',
                 'meshes/2d/square_quad.mesh',
                 'meshes/2d/square_unit_tri.mesh']:
        tst.report(name)

        u = prepare_variable(op.join(data_dir, name), n_components=5)

        bbox = u.field.domain.get_mesh_bounding_box()
        coors = nm.c_[tuple([nm.linspace(ii[0] + 1e-3, ii[1] - 1e-3, 100)
                             for ii in bbox.T])]

        grad, cells, status = u.evaluate_at(coors, mode='grad',
                                            close_limit=0.0,
                                            ret_status=True)
        agrad = nm.dot(grad[:, :, :], nm.ones((grad.shape[2], 1)))[..., 0]

        eps = 1e-5
        val0 = u.evaluate_at(coors - eps, close_limit=0.0)[..., 0]
        val1 = u.evaluate_at(coors + eps, close_limit=0.0)[..., 0]

        ngrad = 0.5 * (val1 - val0) / eps

        ii = nm.where(status == 0)

        tst.report('max. difference:', nm.abs(agrad[ii] - ngrad[ii]).max())

        _ok = nm.allclose(agrad[ii], ngrad[ii], rtol=0.0, atol=10 * eps)
        tst.report('->', _ok)

        ok = ok and _ok

        for ic in range(1, u.n_components):
            _ok = nm.allclose((ic + 1) * agrad[ii, 0], agrad[ii, ic],
                              rtol=0.0, atol=1e-12)
            tst.report('component %d / component 0: mean: %.2f'
                        % (ic, (agrad[ii, ic] / agrad[ii, 0]).mean()))
            tst.report('->', _ok)

            ok = ok and _ok

    assert ok

def test_evaluate_at():
    from sfepy import data_dir
    from sfepy.discrete.fem import Mesh
    from sfepy.discrete import Variables
    from sfepy.discrete.fem import FEDomain, Field

    meshes = {
        'tp' : Mesh.from_file(data_dir + '/meshes/3d/block.mesh'),
    }
    datas = gen_datas(meshes)

    fields = {
        'scalar_tp' : ((1,1), 'Omega', 1),
        'vector_tp' : ((3,1), 'Omega', 1),
    }

    ok = True
    for field_name in ['scalar_tp', 'vector_tp']:
        d = FEDomain('d', meshes['tp'])
        d.create_region('Omega', 'all')

        f = fields[field_name]
        field = Field.from_args('f', nm.complex128, f[0],
                                d.regions[f[1]],
                                approx_order=f[2])
        ff = {field.name : field}

        vv = Variables.from_conf(transform_variables(variables), ff)
        u = vv['u']

        bbox = d.get_mesh_bounding_box()
        t = nm.expand_dims(nm.linspace(0, 1, 100), 1)
        coors = nm.expand_dims(bbox[1] - bbox[0], 0) * t + bbox[0]

        data_r = datas[field_name]
        data_i = 2. / (1 + datas[field_name])

        u.set_from_mesh_vertices(data_r)
        vals_r = u.evaluate_at(coors)
        u.set_from_mesh_vertices(data_i)
        vals_i = u.evaluate_at(coors)
        u.set_from_mesh_vertices(data_r + data_i * 1j)
        vals = u.evaluate_at(coors)

        _ok = nm.allclose(vals_r + vals_i * 1j, vals, rtol=0.0, atol=1e-12)
        _ok = _ok and nm.abs(vals).sum() > 1
        tst.report('evaluating complex field %s: %s' % (field_name, _ok))

        ok = ok and _ok

    assert ok
