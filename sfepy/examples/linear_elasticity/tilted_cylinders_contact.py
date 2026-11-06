r"""
Impact of two elastic cylinders with crossed axes (the tilted cylinder contact
of [1]_, Sec. 7.5): explicit dynamics with the penalty contact or with the
bipenalty predictor-corrector scheme, frictionless or with the Coulomb
friction, using the contactino library.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} \rho \ul{v} \cdot \ddot{\ul{u}}
    + \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    + \int_{\Gamma_{c}} \left(\varepsilon_N \langle g_N(\ul{u}) \rangle \ul{n}
    + \ul{t}_T(\ul{u})\right) \cdot \ul{v}
    = 0
    \;, \quad \forall \ul{v} \;.

The master cylinder (axis :math:`z`, Young's modulus 1000) is at rest, the
slave cylinder (axis :math:`y`, Young's modulus 500) moves towards it along
:math:`x` with the velocity ``v0`` (and optionally along :math:`z` with
``v_t``, which makes the impact oblique and the friction visible). Both have
the radius 1, the length 2, the Poisson's ratio 0.3 and the density 1. The
meshes are O-grid ("butterfly") hexahedral meshes of solid cylinders, 5632
elements each by default (11520 in [1]_).

With ``method='bipenalty_pc'``, the stepper ``ts.central_difference_pc``
predicts the step without contact and the contactino library corrects the
predicted configuration (``correction='terms'``): the forward-increment
Lagrange multiplier correction, i.e. the exact limit of the stabilized
bipenalty scheme of [1]_, which does not depend on the penalty parameters.
With friction, the correction also removes the slip increment of the points
in contact relative to the committed configuration, inside the Coulomb cone.
With ``method='penalty'``, the stepper ``ts.central_difference`` uses the
penalty contact forces directly.

The time step is ``dt_factor`` times the critical time step :math:`2 /
\omega_{max}` of the discretization without contact, computed from the
largest eigenvalue of the stiffness matrix with respect to the lumped mass
matrix.

Notes on the setting of [1]_: the text gives the radius 2 and the height 4,
but the axes of its Figure 28 span :math:`[-1, 1]` and the contact force
history (impulse about 30 N s) corresponds to the mass :math:`\pi r^2 h \rho
= 2 \pi` of the radius 1 and the length 2, so these are used here. The
minimum kinetic energy of about one half of the initial one in Figure 28
corresponds to one cylinder moving and the other at rest.

.. [1] Y.-J. Kwon, J.-G. Kim, S. S. Cho, J. A. González, A stabilized
       bipenalty formulation using a predictor-corrector scheme for explicit
       contact-impact analysis, Int. J. Numer. Meth. Engng 126 (2025)
       e70119.

Usage examples::

  # The bipenalty predictor-corrector, frictionless (as in [1]_).
  sfepy-run sfepy/examples/linear_elasticity/tilted_cylinders_contact.py

  # With the Coulomb friction and an oblique impact, in the mortar
  # discretization: with the default one (quadrature points), the tangential
  # correction of contactino does not converge.
  sfepy-run sfepy/examples/linear_elasticity/tilted_cylinders_contact.py -d "mu=0.3, v_t=2.0, algorithm='two_pass,discretization=mortar', output_dir='output/tilted_cylinders_mu03'"

  # The penalty method for comparison.
  sfepy-run sfepy/examples/linear_elasticity/tilted_cylinders_contact.py -d "method='penalty', eps_n=1e4, output_dir='output/tilted_cylinders_penalty'"

  sfepy-view output/tilted_cylinders/cylinders.*.vtk -f u:wu:f1:p0 1:vw:p0
"""
import os.path as op
import sys
from functools import partial
import inspect

import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.ioutils import markdown_table_from_dict
import sfepy.mechanics.matcoefs as mc
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_solid_cylinder_mesh

def gen_cylinders_mesh(radius, length, shape, gap0):
    """
    The master cylinder (group 0, axis z, centre at the origin) and the slave
    cylinder (group 1, axis y, centre at x = 2 radius + gap0).
    """
    master = gen_solid_cylinder_mesh(radius, length, shape, (0.0, 0.0, 0.0),
                                     axis='z', mat_id=0, verbose=False)
    slave = gen_solid_cylinder_mesh(radius, length, shape,
                                    (2.0 * radius + gap0, 0.0, 0.0),
                                    axis='y', mat_id=1, verbose=False)
    mesh = master + slave
    mesh.name = 'cylinders'
    return mesh

def get_max_eigenfrequency(mesh, conf_dict):
    """
    The largest eigenfrequency of the discretization without contact: the
    square root of the largest eigenvalue of the stiffness matrix with
    respect to the lumped mass matrix.
    """
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.discrete import Problem
    from sfepy.solvers.ts_solvers import get_max_eigenfrequency as get_omega

    def mesh_hook(mesh_, mode):
        if mode == 'read':
            return mesh

    conf_dict = dict(conf_dict, filename_mesh=UserMeshIO(mesh_hook),
                     equations={'K' : conf_dict['stiffness'] + ' = 0'},
                     ebcs={}, solvers={}, options={})
    required, other = get_standard_keywords()
    conf = ProblemConf.from_dict(conf_dict, sys.modules[__name__],
                                 required=[], other=required + other)
    pb = Problem.from_conf(conf, init_solvers=False)
    pb.time_update()
    pb.update_materials()
    mtx_k = pb.evaluate(conf.stiffness, mode='weak', dw_mode='matrix',
                        verbose=False)
    mtx_m = pb.evaluate(conf.mass, mode='weak', dw_mode='matrix',
                        verbose=False)
    return get_omega(mtx_k, mtx_m)

def define(
        method='bipenalty_pc',
        mu=0.0,
        v0=5.0,
        v_t=0.0,

        radius=1.0,
        length=2.0,
        shape=(17, 23),
        gap0=1e-3,

        young_m=1000.0,
        young_s=500.0,
        poisson=0.3,
        rho=1.0,

        eps_n=1e2,
        eps_t=None,
        algorithm='two_pass',

        dt_factor=0.9,
        t1=0.3,
        mass_lumping='row_sum',
        fast_rmm=False,
        save_times=31,

        output_dir='output/tilted_cylinders',
        verbose=True,
):
    """
    Parameters
    ----------
    method : 'bipenalty_pc' or 'penalty'
        The bipenalty predictor-corrector scheme with the contactino
        correction, or the penalty contact with central differences.
    mu : float
        The Coulomb friction coefficient (0 = frictionless).
    v0, v_t : float
        The initial velocity of the slave cylinder towards the master one
        (along -x) and along the master axis (z).
    radius, length, shape, gap0 :
        The geometry of both cylinders, the mesh shape (nodes along the sides
        of the cross-section square, nodes along the axis) and the initial
        gap.
    young_m, young_s, poisson, rho : float
        The material parameters of the master and slave cylinders.
    eps_n, eps_t : float
        The normal and tangential penalty parameters (eps_t defaults to
        eps_n). The predictor-corrector correction does not depend on them.
    algorithm : 'two_pass' or 'two_half_pass'
        The contact algorithm of contactino; the two-full-pass gives exact
        action and reaction. Append ``',discretization=mortar'`` for the
        mortar discretization (the L2-projected gap), needed with friction.
    dt_factor : float
        The time step in units of the critical time step without contact.
    t1 : float
        The final time.
    mass_lumping : str
        The mass lumping of ``de_mass``.
    fast_rmm : bool
        Use the reciprocal mass matrix solver ``ls.rmm``.
    save_times : int
        The number of saved time steps.
    """
    args = locals()

    inodir = partial(op.join, output_dir)

    output.set_output(filename=inodir('output_log.txt'), quiet=not verbose,
                      combined=True)

    if method not in ('bipenalty_pc', 'penalty'):
        raise ValueError("method must be 'bipenalty_pc' or 'penalty'! (%s)"
                         % method)

    dim = 3
    if eps_t is None:
        eps_t = eps_n
    mesh = gen_cylinders_mesh(radius, length, shape, gap0)
    xc = 2.0 * radius + gap0 # The slave axis position.

    def mesh_hook(mesh_, mode):
        if mode == 'read':
            return mesh

    filename_mesh = UserMeshIO(mesh_hook)

    eps = 1e-6 * radius

    def get_master_side(coors, domain=None):
        """
        The lateral surface of the master cylinder facing the slave one.
        """
        rr = nm.linalg.norm(coors[:, :2], axis=1)
        return nm.where((nm.abs(rr - radius) < eps)
                        & (coors[:, 0] > 0.2 * radius))[0]

    def get_slave_side(coors, domain=None):
        """
        The lateral surface of the slave cylinder facing the master one.
        """
        rr = nm.hypot(coors[:, 0] - xc, coors[:, 2])
        return nm.where((nm.abs(rr - radius) < eps)
                        & (coors[:, 0] < xc - 0.2 * radius))[0]

    regions = {
        'Omega' : 'all',
        'Omega0' : 'cells of group 0',
        'Omega1' : 'cells of group 1',
        'Contact0' : ('vertices by get_master_side', 'facet', 'Omega0'),
        'Contact1' : ('vertices by get_slave_side', 'facet', 'Omega1'),
        'Contact' : ('r.Contact0 +s r.Contact1', 'facet'),
    }

    fields = {
        'displacement' : ('real', 'vector', 'Omega', 1),
    }

    # With auto_transform_equations=True, the solver creates du, ddu, dv, ddv.
    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    materials = {
        'solid' : ({
            'D0' : mc.stiffness_from_youngpoisson(dim, young_m, poisson),
            'D1' : mc.stiffness_from_youngpoisson(dim, young_s, poisson),
            'rho' : rho,
            '.lumping' : mass_lumping,
            '.beta' : 1.0,
        },),
        'contact' : ({
            '.epss' : eps_n,
            '.epst' : eps_t,
            '.mu' : mu,
            '.algorithm' : algorithm,
        },),
    }

    integrals = {
        'i' : 2,
    }

    stiffness = ('dw_lin_elastic.i.Omega0(solid.D0, v, u)'
                 ' + dw_lin_elastic.i.Omega1(solid.D1, v, u)')
    mass = 'de_mass.i.Omega(solid.rho, solid.lumping, solid.beta, v, u)'

    # The critical time step of the discretization without contact.
    omega_max = get_max_eigenfrequency(mesh, dict(
        regions={key : val for key, val in regions.items()
                 if key.startswith('Omega')},
        fields=fields, variables=variables, materials=materials,
        integrals=integrals, stiffness=stiffness, mass=mass,
    ))
    dt_crit = 2.0 / omega_max
    dt = dt_factor * dt_crit

    mass_m = rho * nm.pi * radius**2 * length

    signature = inspect.signature(define)
    arg_names = list(signature.parameters.keys())
    options = {name : args[name] for name in arg_names}
    options.update(omega_max=omega_max, dt_crit=dt_crit, dt=dt,
                   n_step=int(nm.ceil(t1 / dt)), n_el=mesh.n_el)
    with open(inodir('options.md'), 'w') as fd:
        fd.write(markdown_table_from_dict(options))

    output('tilted cylinders: %d cells, omega_max = %.4e, dt_crit = %.4e,'
           ' dt = %.4e, %d steps' % (mesh.n_el, omega_max, dt_crit, dt,
                                     int(nm.ceil(t1 / dt))))

    def get_ic_u(coors, ic):
        return nm.zeros_like(coors)

    def get_ic_du(coors, ic):
        """
        The slave cylinder moves with (-v0, 0, v_t), the master one is at
        rest.
        """
        val = nm.zeros_like(coors)
        slave = coors[:, 0] > radius + 0.5 * gap0
        val[slave, 0] = -v0
        val[slave, 2] = v_t
        return val

    functions = {
        'get_master_side' : (get_master_side,),
        'get_slave_side' : (get_slave_side,),
        'get_ic_u' : (get_ic_u,),
        'get_ic_du' : (get_ic_du,),
    }

    ics = {
        'ic' : ('Omega', {'u.all' : 'get_ic_u', 'du.all' : 'get_ic_du'}),
    }

    ebcs = {}

    # The mass term must come first and have ',v,' in the middle of its
    # arguments: transform_equations_ed() replaces ',v,' by ',ddv,'.
    if mu > 0.0:
        contact_term = ('dw_contact_friction.i.Contact(contact.algorithm,'
                        ' contact.epss, contact.epst, contact.mu, v, u)')

    else:
        contact_term = ('dw_contact.i.Contact(contact.algorithm, contact.epss,'
                        ' v, u)')
    contact_name = contact_term.split('.')[0]

    equations = {
        'balance' : ('de_mass.i.Omega(solid.rho, solid.lumping, solid.beta,'
                     ' v, ddu) + ' + stiffness + ' + ' + contact_term
                     + ' = 0'),
    }

    # History: time, contact force on the master cylinder, momentum, angular
    # velocity of the slave cylinder about its axis, energies, the numbers of
    # contact nodes with a correction force and with a positive gap.
    hist = Struct(name='history', rows=[], mtx_m=None, mtx_k=None,
                  inodes=None, areas=None)

    def get_contact_term(pb):
        return [tt for tt in pb.equations['K'].terms
                if tt.name == contact_name][0]

    def step_hook(pb, ts, variables):
        term = get_contact_term(pb)
        u = variables['u']()
        du = variables['du']()

        if hist.mtx_m is None:
            hist.mtx_m = pb.evaluate(mass.replace(', v, u)', ', ddv, ddu)'),
                                     mode='weak', dw_mode='matrix',
                                     verbose=False).tocsr()
            hist.mtx_k = pb.evaluate(stiffness, mode='weak',
                                     dw_mode='matrix', verbose=False).tocsr()
            hist.master = pb.domain.regions['Omega0'].vertices
            hist.slave = pb.domain.regions['Omega1'].vertices

        if method == 'bipenalty_pc':
            # The total contact force M da of the last correction.
            fall = pb.get_solver().contact_force
            if fall is None:
                fall = nm.zeros_like(u)
            fall = fall.reshape((-1, dim))

        else:
            nodes, fc = term.get_nodal_forces()
            fall = nm.zeros((pb.domain.mesh.n_nod, dim))
            fall[nodes] = fc

        force = fall[hist.master].sum(axis=0)
        n_force = int((nm.abs(fall).sum(axis=1) > 0).sum())

        mv = (hist.mtx_m @ du).reshape((-1, dim))
        momentum = mv.sum(axis=0)

        # The angular velocity of the slave cylinder about its axis (y).
        coors = pb.domain.mesh.coors + u.reshape((-1, dim))
        ii = hist.slave
        mi = hist.mtx_m.diagonal().reshape((-1, dim))[ii, 0]
        com = (mi[:, None] * coors[ii]).sum(axis=0) / mi.sum()
        rel = coors[ii] - com
        vv = du.reshape((-1, dim))[ii]
        ly = (mi * (rel[:, 2] * vv[:, 0] - rel[:, 0] * vv[:, 2])).sum()
        iy = (mi * (rel[:, 0]**2 + rel[:, 2]**2)).sum()

        n_gap = 0
        if term.ci is not None:
            nodes, gap, pressure = term.get_nodal_fields()
            n_gap = int((nm.isfinite(gap) & (gap > 0)).sum())

        ek = 0.5 * du @ (hist.mtx_m @ du)
        es = 0.5 * u @ (hist.mtx_k @ u)
        hist.rows.append((ts.time, *force, *momentum, ly / iy, ek, es,
                          n_force, n_gap))

    def get_nodal_areas(pb):
        """
        The lumped areas of the contact surface nodes (a quarter of the
        adjacent quadrilateral facet areas in the initial configuration).
        """
        if hist.get('areas') is None:
            cmesh = pb.domain.cmesh
            cmesh.setup_connectivity(2, 0)
            conn = cmesh.get_conn(2, 0)
            facets = pb.domain.regions['Contact'].facets
            fareas = cmesh.get_volumes(2)[facets]
            areas = nm.zeros(pb.domain.mesh.n_nod, dtype=nm.float64)
            for facet, area in zip(facets, fareas):
                verts = conn.indices[conn.offsets[facet]:conn.offsets[facet + 1]]
                areas[verts] += area / len(verts)
            hist.areas = areas
        return hist.areas

    def post_process(out, pb, state, extend=False):
        term = get_contact_term(pb)
        n_nod = pb.domain.mesh.n_nod
        stress = nm.zeros((pb.domain.mesh.n_el, 1, 6, 1), dtype=nm.float64)
        for ir, dname in (('Omega0', 'D0'), ('Omega1', 'D1')):
            aux = pb.evaluate('ev_cauchy_stress.i.%s(solid.%s, u)'
                              % (ir, dname), mode='el_avg',
                              copy_materials=False, verbose=False)
            stress[pb.domain.regions[ir].cells] = aux
        out['cauchy_stress'] = Struct(name='output_data', mode='cell',
                                      data=stress)

        # The contact forces acting on the bodies: the correction force M da
        # of the predictor-corrector, or the penalty force. The contact
        # pressure is the force magnitude per lumped nodal area (the nodal
        # fields of contactino refer to the penalty traction, which vanishes
        # after the exact correction).
        if method == 'bipenalty_pc':
            fall = pb.get_solver().contact_force
            fall = (nm.zeros((n_nod, dim)) if fall is None
                    else fall.reshape((-1, dim)))

        else:
            fall = nm.zeros((n_nod, dim))
            if term.ci is not None:
                nodes, fc = term.get_nodal_forces()
                fall[nodes] = fc
        areas = get_nodal_areas(pb)
        pressure = nm.zeros((n_nod, 1), dtype=nm.float64)
        ii = areas > 0
        pressure[ii, 0] = nm.linalg.norm(fall[ii], axis=1) / areas[ii]
        out['f_contact'] = Struct(name='output_data', mode='vertex',
                                  data=fall)
        out['p_contact'] = Struct(name='output_data', mode='vertex',
                                  data=pressure)
        return out

    def post_process_final(pb, state=None):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        hh = nm.array(hist.rows)
        t, fx, fy, fz = hh[:, :4].T
        px, py, pz, omega_y, ek, es = hh[:, 4:10].T
        e0 = ek[0] + es[0]
        fn = -fx # The normal force on the master pushes it along -x.
        ii = nm.argmax(fn)
        impulse = nm.trapezoid(fn, t)

        output('tilted cylinders (%s, mu = %g, v_t = %g):' % (method, mu, v_t))
        output('  peak normal force %.4e at t = %.4e, impulse %.4e'
               ' (rigid-body estimate m v0 = %.4e)'
               % (fn[ii], t[ii], impulse, mass_m * v0))
        output('  kinetic energy: min %.4f, final %.4f of the initial'
               % (ek.min() / e0, ek[-1] / e0))
        output('  total energy: max |E - E0| / E0 = %.3e, final %.3e'
               % (nm.abs(ek + es - e0).max() / e0, (ek[-1] + es[-1]) / e0
                  - 1.0))
        output('  momentum change: max |p - p0| / (m v0) = %.3e'
               % (nm.abs(hh[:, 4:7] - hh[0, 4:7]).max() / (mass_m * v0)))
        output('  final angular velocity of the slave cylinder about y:'
               ' %.4e rad/s' % omega_y[-1])
        if mu > 0.0:
            fmax = nm.abs(fz).max()
            output('  max tangential force %.4e, max |F_z| / (mu F_n) at the'
                   ' peak %.4f' % (fmax, abs(fz[ii]) / max(mu * fn[ii],
                                                           1e-300)))

        fig, axs = plt.subplots(1, 3, figsize=(14, 4))
        ax = axs[0]
        ax.plot(t, (ek + es) / e0, 'k-', label='total')
        ax.plot(t, es / e0, 'r:', label='internal')
        ax.plot(t, ek / e0, 'b--', label='kinetic')
        ax.set_ylim(0, 1.5)
        ax.set_ylabel('energy / initial energy')
        ax = axs[1]
        ax.plot(t, fn, 'k-', label='normal $-F_x$')
        if mu > 0.0:
            ax.plot(t, fz, 'g-', label='tangential $F_z$')
            ax.plot(t, mu * fn, 'k:', label=r'$\mu F_n$')
        ax.set_ylabel('contact force on the master cylinder')
        ax = axs[2]
        ax.plot(t, omega_y, 'm-')
        ax.set_ylabel('angular velocity of the slave cylinder')
        for ax in axs:
            ax.set_xlabel('time')
            ax.grid(True)
            if ax.get_legend_handles_labels()[0]:
                ax.legend()
        fig.suptitle('%s, mu = %g, v0 = %g, v_t = %g, %s'
                     % (method, mu, v0, v_t, algorithm))
        fig.tight_layout()
        fig.savefig(inodir('history.png'), dpi=150)
        plt.close(fig)

        nm.savetxt(inodir('history.txt'), hh,
                   header='t F_x F_y F_z (on the master) p_x p_y p_z'
                   ' omega_y(slave) e_kin e_strain n_force_nodes'
                   ' n_positive_gap_nodes')

    if fast_rmm:
        ls_name = 'lsrmm'

    else:
        ls_name = 'lsd'

    solvers = {
        'lsd' : ('ls.scipy_direct', {
            'use_presolve' : True,
        }),
        'lsrmm' : ('ls.rmm', {
            'rmm_term' : """de_mass.i.Omega(solid.rho, solid.lumping,
                                            solid.beta, ddv, ddu)""",
        }),
        'newton' : ('nls.newton', {
            'i_max' : 1,
            'eps_a' : 1e-6,
            'eps_r' : 1.0,
        }),
        'tscd' : ('ts.central_difference_pc' if method == 'bipenalty_pc'
                  else 'ts.central_difference', {
            'contact_terms' : [contact_name],
            'contact_mass_terms' : [],
            'correction' : 'terms',
            't0' : 0.0,
            't1' : t1,
            'dt' : dt,
            'n_step' : None,
            'is_linear' : False,
            'var_names' : {'u' : 'u', 'du' : 'du', 'ddu' : 'ddu'},
            'verbose' : 1,
        }),
    }

    options = {
        'ts' : 'tscd',
        'nls' : 'newton',
        'ls' : ls_name,
        'save_times' : save_times,
        'auto_transform_equations' : True,
        'active_only' : False,
        'output_format' : 'vtk',
        'output_dir' : output_dir,
        'step_hook' : 'step_hook',
        'post_process_hook' : 'post_process',
        'post_process_hook_final' : 'post_process_final',
    }

    return locals()
