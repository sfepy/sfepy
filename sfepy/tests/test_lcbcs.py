import os.path as op
import numpy as nm
import pytest

import sfepy.base.testing as tst

def test_stokes_slip_bc(output_dir):
    import scipy.sparse as sp

    from sfepy.base.conf import ProblemConf
    from sfepy.discrete import Problem
    import sfepy.examples.navier_stokes.stokes_slip_bc as ssb

    conf = ProblemConf.from_dict(
        ssb.define(output_dir=output_dir, save_lcbc_vecs=True),
        ssb,
    )
    pb = Problem.from_conf(conf, init_solvers=False)

    pb.time_update()
    variables = pb.get_variables()

    adi = variables.adi
    lcdi = variables.lcdi
    mtx = variables.mtx_lcbc

    ok = adi.var_names == lcdi.var_names
    tst.report('same adi-lcdi ordering:', ok)

    ublock = mtx[adi.indx['u']]
    ir, ic = ublock.nonzero()
    ir += adi.indx['u'].start

    i0, i1 = adi.indx['u'].start, adi.indx['u'].stop
    _ok0 = (i0 <= ir).all() and  (ir < i1).all()
    tst.report('u block rows in [%d %d[: %s' % (i0, i1, _ok0))

    i0, i1 = lcdi.indx['u'].start, lcdi.indx['u'].stop
    _ok1 = (i0 <= ic).all() and  (ic < i1).all()
    tst.report('u block cols in [%d %d[: %s' % (i0, i1, _ok1))

    ok = ok and _ok0 and _ok1

    pblock = mtx[adi.indx['p']]
    ir, ic, iv = sp.find(pblock)
    ir += adi.indx['p'].start

    i0, i1 = adi.indx['p'].start, adi.indx['p'].stop
    _ok0 = (i0 <= ir).all() and  (ir < i1).all()
    tst.report('p block rows in [%d %d[: %s' % (i0, i1, _ok0))

    i0, i1 = lcdi.indx['p'].start, lcdi.indx['p'].stop
    _ok1 = (i0 <= ic).all() and  (ic < i1).all()
    tst.report('p block cols in [%d %d[: %s' % (i0, i1, _ok1))

    ok = ok and _ok0 and _ok1

    _ok0 = (len(ir) == adi.n_dof['p'])
    tst.report('p block size correct:', _ok0)
    _ok1 = ((ir - adi.indx['p'].start) == (ic - lcdi.indx['p'].start)).all()
    tst.report('p block diagonal:', _ok1)
    _ok2 = (iv == 1.0).all()
    tst.report('p block identity:', _ok2)

    ok = ok and _ok0 and _ok1 and _ok2

    assert ok

def test_laplace_shifted_periodic(output_dir):
    from sfepy.mesh.mesh_generators import gen_block_mesh
    from sfepy.discrete.fem import FEDomain
    from sfepy.examples.diffusion.laplace_shifted_periodic import run

    dims = [2.0, 1.0]
    shape = [21, 11]
    centre = [0.0, 0.0]
    mesh = gen_block_mesh(dims, shape, centre, name='block-fem')
    domain = FEDomain('test_laplace_shifted_periodic', mesh)

    pb, state = run(domain, 3, output_dir=output_dir)

    gamma3 = pb.domain.regions['Gamma3']
    gamma4 = pb.domain.regions['Gamma4']

    field = pb.fields['fu']

    # Check that the shift equals to one.
    i3 = field.get_dofs_in_region(gamma3, merge=True)
    i4 = field.get_dofs_in_region(gamma4, merge=True)

    i_corners = nm.array([0, shape[0] - 1])
    ii = nm.setdiff1d(nm.arange(len(i3)), i_corners)

    vals = state()

    shift = vals[i3] - vals[i4]

    ok = (shift[i_corners] == 0.0).all()

    ok = ok and nm.allclose(shift[ii], 1.0, rtol=0.0, atol=1e-14)

    if not ok:
        tst.report('wrong shift:', shift)

    assert ok

@pytest.mark.parametrize('mesh_filename', [
    'meshes/2d/special/circle_in_square.mesh',
    'meshes/3d/special/cube_sphere.mesh',
])
def test_elasticity_rigid(mesh_filename, output_dir):
    from sfepy import data_dir
    from sfepy.base.base import IndexedStruct
    from sfepy.discrete import (FieldVariable, Material, Integral,
                                Equation, Equations, Problem)
    from sfepy.discrete.fem import Mesh, FEDomain, Field
    from sfepy.mechanics.matcoefs import stiffness_from_lame
    from sfepy.terms import Term
    from sfepy.discrete.conditions import (Conditions, EssentialBC,
                                           LinearCombinationBC)
    from sfepy.solvers.ls import ScipyDirect
    from sfepy.solvers.nls import Newton

    filename = op.join(data_dir, mesh_filename)
    mesh = Mesh.from_file(filename)
    domain = FEDomain('domain', mesh)

    min_x, max_x = domain.get_mesh_bounding_box()[:,0]
    eps = 1e-8 * (max_x - min_x)

    axis = {2 : 'y', 3 : 'z'}[mesh.dim]
    order = {2 : 2, 3 : 1}[mesh.dim]

    omega = domain.create_region('Omega', 'all')
    yrig = domain.create_region('Yrig', 'cells of group 2')
    bottom = domain.create_region('Bottom',
                                  f'vertices in {axis} < {min_x + eps}',
                                  'facet')
    top = domain.create_region('Top',
                               f'vertices in {axis} > {max_x - eps}',
                               'facet')
    field = Field.from_args('fu', nm.float64, 'vector', omega,
                            approx_order=order)
    u = FieldVariable('u', 'unknown', field)
    v = FieldVariable('v', 'test', field, primary_var_name='u')

    m = Material('m', D=stiffness_from_lame(dim=mesh.dim, lam=10.0, mu=1.0))

    integral = Integral('i', order=2 * order)
    t1 = Term.new('dw_lin_elastic(m.D, v, u)',
                  integral, omega, m=m, v=v, u=u)
    eq = Equation('balance', t1)
    eqs = Equations([eq])

    fix = EssentialBC('fix', bottom, {'u.all' : 0.0})
    load = EssentialBC('load', top, {'u.all' : 0.3})
    rig = LinearCombinationBC('rig', [yrig, None], {'u.all' : None},
                              None, 'rigid')

    ls = ScipyDirect({})

    nls_status = IndexedStruct()
    nls = Newton({}, lin_solver=ls, status=nls_status)

    pb = Problem('elasticity', equations=eqs)

    trunk = f'test_elasticity_rigid_{mesh.dim}d'
    pb.setup_output(output_filename_trunk=trunk, output_dir=output_dir)
    pb.save_regions_as_groups(op.join(output_dir, trunk + '_regions'))

    pb.set_bcs(ebcs=Conditions([fix, load]), lcbcs=Conditions([rig]))
    pb.set_solver(nls)

    status = IndexedStruct()
    pb.solve(status=status)

    assert nls_status.condition == 0

    strain = u.evaluate('cauchy_strain', region=yrig, integral=integral)

    mstrain = nm.abs(strain).max()
    tst.report('max |strain| in rigid region:', mstrain)

    assert mstrain < 1e-13
