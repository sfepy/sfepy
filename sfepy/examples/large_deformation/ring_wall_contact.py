r"""
Impact of a hyperelastic ring on a rigid wall with the Coulomb friction (the
ring contact of [1]_, Sec. 7.2): explicit dynamics in the total Lagrangian
formulation with the bipenalty predictor-corrector scheme and the contactino
library, in 2D plane stress.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} \rho \ul{v} \cdot \ddot{\ul{u}}
    + \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u}; \ul{v})
    + \int_{\Gamma_{c}} \left(\varepsilon_N \langle g_N(\ul{u}) \rangle \ul{n}
    + \ul{t}_T(\ul{u})\right) \cdot \ul{v}
    = 0
    \;, \quad \forall \ul{v} \;,

where :math:`S_{ij} = \mu (\delta_{ij} - C^{-1}_{ij}) + \lambda \ln J
C^{-1}_{ij}` is the second Piola-Kirchhoff stress of the compressible
neo-Hookean material in plane stress (``dw_tl_he_neohook_compressible``).

The ring (radii 8 and 10 mm, 600 four-node elements, :math:`\lambda = 171.41`
MPa, :math:`\mu = 428.51` MPa) starts just above a rigid wall with the
velocity :math:`(20, -7)` mm/s and bounces off it while sliding along it. With
friction (``mu=0.9`` in [1]_), the ring starts to rotate. The stepper
``ts.central_difference_pc`` predicts every step without contact and
corrects it by the contact: with ``correction='host'`` by the stiffness and
inertia penalty forces of the bipenalty method with Criterion 1 of [1]_ (the
method of [1]_, frictionless here), with ``correction='terms'`` by the
correction of the predicted configuration computed by the contactino library
(the forward-increment Lagrange multiplier correction, the exact limit of the
scheme, with the Coulomb friction). The wall is a layer of elements with all
nodes fixed (displacements, velocities and accelerations), i.e. immovable in
the correction.

The density is not given in [1]_. The response depends on it only through the
ratio of the kinetic energy to the stiffness; ``rho`` is chosen so that the
frictionless energy history agrees with Figure 16 of [1]_ (the internal energy
reaches about 0.09 of the initial one around t = 1.3 s).

.. [1] Y.-J. Kwon, J.-G. Kim, S. S. Cho, J. A. González, A stabilized
       bipenalty formulation using a predictor-corrector scheme for explicit
       contact-impact analysis, Int. J. Numer. Meth. Engng 126 (2025)
       e70119.

Usage examples::

  # Frictionless.
  sfepy-run sfepy/examples/large_deformation/ring_wall_contact.py

  # With friction.
  sfepy-run sfepy/examples/large_deformation/ring_wall_contact.py -d "mu=0.9, output_dir='output/ring_wall_mu09'"

  sfepy-view output/ring_wall/ring.*.vtk -2 -f u:wu:f1:p0 1:vw:p0
"""
import os.path as op
import sys
from functools import partial
import inspect

import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.ioutils import markdown_table_from_dict
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mechanics.matcoefs import stiffness_from_lame
from sfepy.mesh.mesh_generators import gen_block_mesh, gen_cylinder_mesh

def gen_ring_wall_mesh(r1, r2, shape, wall_dims, wall_shape, gap0):
    """
    The ring (group 0, centre at the origin) and the wall (group 1) with the
    top at y = -r2 - gap0.
    """
    ring = gen_cylinder_mesh([r1, r1, r2, r2, 0.0], [shape[0], shape[1], 1],
                             [0.0, 0.0, 0.0], axis='z', make_2d=True,
                             name='ring', verbose=False)
    ring.cmesh.cell_groups[:] = 0
    x0, x1, height = wall_dims
    top = -r2 - gap0
    wall = gen_block_mesh((x1 - x0, height), wall_shape,
                          (0.5 * (x0 + x1), top - 0.5 * height), mat_id=1,
                          name='wall', verbose=False)
    mesh = ring + wall
    mesh.name = 'ring'
    return mesh

def get_ring_matrices(mesh, conf_dict):
    """
    The tangent stiffness matrix of the ring in the reference configuration
    and the lumped mass matrix, for the critical time step.
    """
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.discrete import Problem

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
    # The tangent of the hyperelastic term at the zero displacement.
    pb.get_initial_state()
    mtx_k = pb.evaluate(conf.stiffness, mode='weak', dw_mode='matrix',
                        verbose=False)
    mtx_m = pb.evaluate(conf.mass, mode='weak', dw_mode='matrix',
                        verbose=False)
    return mtx_k, mtx_m

def define(
        mu=0.0,
        correction='terms',
        v0=(20.0, -7.0),

        r1=8.0,
        r2=10.0,
        shape=(6, 120),
        gap0=1e-3,
        wall_dims=(-15.0, 90.0, 2.0),
        wall_shape=(106, 2),

        lam=171.41,
        mu_e=428.51,
        rho=0.08,

        eps_n=1e10,
        eps_t=None,
        algorithm='one_pass:Contact0',

        dt_factor=0.9,
        t1=3.0,
        save_times=61,

        output_dir='output/ring_wall',
        verbose=True,
):
    """
    Parameters
    ----------
    mu : float
        The Coulomb friction coefficient (0 = frictionless).
    correction : 'host' or 'terms'
        The corrector of the predictor-corrector scheme: assembled by the
        host from the penalty contact force and the inertia penalty mass
        with Criterion 1 of [1]_ (eps_M = eps_N dt**2 - m, m the surface
        nodal mass of the ring per unit length), frictionless only; or the
        correction of the predicted configuration by the contactino library
        (exact, with friction).
    v0 : tuple of 2 floats
        The initial velocity of the ring.
    r1, r2, shape : float, float, tuple of 2 ints
        The radii of the ring and its mesh shape (nodes across the ring,
        nodes along the circumference = elements along it).
    gap0 : float
        The initial gap between the ring and the wall.
    wall_dims, wall_shape : tuples
        The extent (x0, x1, height) and the mesh shape of the wall.
    lam, mu_e, rho : float
        The Lame parameters and the density of the ring.
    eps_n, eps_t : float
        The normal and tangential penalty parameters (eps_t defaults to
        eps_n); the predictor-corrector correction does not depend on them.
    algorithm : str
        The contact algorithm of contactino, see ``dw_contact``; by default
        the one-pass algorithm with the ring as the slave.
    dt_factor : float
        The time step in units of the critical time step of the ring in the
        reference configuration.
    t1 : float
        The final time.
    save_times : int
        The number of saved time steps.
    """
    args = locals()

    inodir = partial(op.join, output_dir)

    output.set_output(filename=inodir('output_log.txt'), quiet=not verbose,
                      combined=True)

    dim = 2
    if eps_t is None:
        eps_t = eps_n
    v0 = nm.asarray(v0, dtype=nm.float64)
    mesh = gen_ring_wall_mesh(r1, r2, shape, wall_dims, wall_shape, gap0)
    wall_top = -r2 - gap0

    def mesh_hook(mesh_, mode):
        if mode == 'read':
            return mesh

    filename_mesh = UserMeshIO(mesh_hook)

    eps = 1e-6 * r2

    def get_ring_outer(coors, domain=None):
        rr = nm.linalg.norm(coors, axis=1)
        return nm.where(nm.abs(rr - r2) < eps)[0]

    regions = {
        'Omega' : 'all',
        'Omega0' : 'cells of group 0',
        'Omega1' : 'cells of group 1',
        'Contact0' : ('vertices by get_ring_outer', 'facet', 'Omega0'),
        'Contact1' : ('vertices in (y > %.12e) *v r.Omega1'
                      % (wall_top - eps), 'facet', 'Omega1'),
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
        'ring' : ({
            'lm' : nm.array([[lam, mu_e]]),
            'rho' : rho,
            '.lumping' : 'row_sum',
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

    stiffness = 'dw_tl_he_neohook_compressible.i.Omega0(ring.lm, v, u)'
    mass = 'de_mass.i.Omega(ring.rho, ring.lumping, ring.beta, v, u)'

    # The critical time step of the ring in the reference configuration.
    from sfepy.solvers.ts_solvers import get_max_eigenfrequency
    mtx_k, mtx_m = get_ring_matrices(mesh, dict(
        regions={'Omega' : 'all', 'Omega0' : 'cells of group 0'},
        fields=fields, variables=variables,
        materials={'ring' : materials['ring']}, integrals=integrals,
        stiffness=stiffness, mass=mass,
    ))
    omega_max = get_max_eigenfrequency(mtx_k, mtx_m)
    dt_crit = 2.0 / omega_max
    dt = dt_factor * dt_crit

    if (mu > 0.0) and (correction != 'terms'):
        raise ValueError("friction needs correction='terms'!")

    # Criterion 1 of [1]_ (flexible-rigid contact): the surface nodal mass of
    # the ring per unit length of its outer surface.
    rr = nm.linalg.norm(mesh.coors, axis=1)
    outer = nm.where(nm.abs(rr - r2) < 1e-6 * r2)[0]
    m_node = mtx_m.diagonal().reshape((-1, dim))[outer, 0].max()
    m_surf = m_node / (2.0 * nm.pi * r2 / shape[1])
    eps_m = eps_n * dt**2 - m_surf
    if eps_m < 0.0:
        raise ValueError('eps_n too small for Criterion 1: eps_n >= %.3e'
                         ' is needed!' % (m_surf / dt**2))
    materials['contact'][0]['.epsm'] = eps_m
    area = nm.pi * (r2**2 - r1**2)
    mass_ring = rho * area

    signature = inspect.signature(define)
    arg_names = list(signature.parameters.keys())
    options = {name : args[name] for name in arg_names}
    options.update(omega_max=omega_max, dt_crit=dt_crit, dt=dt,
                   n_step=int(nm.ceil(t1 / dt)), m_surf=m_surf, eps_m=eps_m)
    with open(inodir('options.md'), 'w') as fd:
        fd.write(markdown_table_from_dict(options))

    output('ring on a wall: omega_max = %.4e, dt_crit = %.4e, dt = %.4e,'
           ' %d steps' % (omega_max, dt_crit, dt, int(nm.ceil(t1 / dt))))

    def get_ic_u(coors, ic):
        return nm.zeros_like(coors)

    def get_ic_du(coors, ic):
        val = nm.zeros_like(coors)
        ring = coors[:, 1] > wall_top + 0.5 * gap0
        val[ring] = v0
        return val

    functions = {
        'get_ring_outer' : (get_ring_outer,),
        'get_ic_u' : (get_ic_u,),
        'get_ic_du' : (get_ic_du,),
    }

    ics = {
        'ic' : ('Omega', {'u.all' : 'get_ic_u', 'du.all' : 'get_ic_du'}),
    }

    ebcs = {
        # The fixed wall: the velocities and accelerations as well, so that
        # the wall is immovable in the contact correction.
        'fix_wall' : ('Omega1', {'u.all' : 0.0, 'du.all' : 0.0,
                                 'ddu.all' : 0.0}),
    }

    if mu > 0.0:
        contact_term = ('dw_contact_friction.i.Contact(contact.algorithm,'
                        ' contact.epss, contact.epst, contact.mu, v, u)')

    else:
        contact_term = ('dw_contact.i.Contact(contact.algorithm, contact.epss,'
                        ' v, u)')
    contact_name = contact_term.split('.')[0]
    if correction == 'host':
        contact_term += (' + dw_contact_mass.i.Contact(contact.algorithm,'
                         ' contact.epsm, v, ddu)')

    # The mass term must come first and have ',v,' in the middle of its
    # arguments: transform_equations_ed() replaces ',v,' by ',ddv,'.
    equations = {
        'balance' : ('de_mass.i.Omega(ring.rho, ring.lumping, ring.beta,'
                     ' v, ddu) + ' + stiffness + ' + ' + contact_term
                     + ' = 0'),
    }

    # History: time, contact force on the wall, kinetic energy, internal
    # energy (the work of the internal forces), rotation of the ring, centre
    # of mass.
    hist = Struct(name='history', rows=[], mtx_m=None, ring=None, ei=0.0,
                  u_last=None, f_last=None)

    def get_internal_force(pb):
        return pb.evaluate(stiffness, mode='weak', dw_mode='vector',
                           active_only=False, verbose=False)

    def step_hook(pb, ts, variables):
        u = variables['u']().copy()
        du = variables['du']()

        if hist.mtx_m is None:
            hist.mtx_m = pb.evaluate(mass.replace(', v, u)', ', ddv, ddu)'),
                                     mode='weak', dw_mode='matrix',
                                     active_only=False,
                                     verbose=False).tocsr()
            hist.ring = pb.domain.regions['Omega0'].vertices
            hist.wall = pb.domain.regions['Omega1'].vertices
            hist.mi = hist.mtx_m.diagonal().reshape((-1, dim))[hist.ring, 0]
            x0 = pb.domain.mesh.coors[hist.ring]
            hist.x0 = x0 - (hist.mi[:, None] * x0).sum(0) / hist.mi.sum()

        fint = get_internal_force(pb)
        if hist.u_last is not None:
            hist.ei += 0.5 * (fint + hist.f_last) @ (u - hist.u_last)
        hist.u_last, hist.f_last = u, fint

        fall = pb.get_solver().contact_force
        fall = (nm.zeros((pb.domain.mesh.n_nod, dim)) if fall is None
                else fall.reshape((-1, dim)))
        # The wall is immovable in the correction: its reaction is the opposite
        # of the contact force on the ring.
        force = -fall[hist.ring].sum(axis=0)

        # The rotation of the ring: the angle of the best rigid rotation of
        # the nodes about the centre of mass.
        xx = pb.domain.mesh.coors[hist.ring] + u.reshape((-1, dim))[hist.ring]
        com = (hist.mi[:, None] * xx).sum(0) / hist.mi.sum()
        rel = xx - com
        cross = (hist.mi * (hist.x0[:, 0] * rel[:, 1]
                            - hist.x0[:, 1] * rel[:, 0])).sum()
        dot = (hist.mi * (hist.x0 * rel).sum(1)).sum()
        angle = nm.arctan2(cross, dot)

        ek = 0.5 * du @ (hist.mtx_m @ du)
        hist.rows.append((ts.time, *force, ek, hist.ei, angle, *com))

    def post_process(out, pb, state, extend=False):
        fall = pb.get_solver().contact_force
        n_nod = pb.domain.mesh.n_nod
        fall = (nm.zeros((n_nod, dim)) if fall is None
                else fall.reshape((-1, dim)))
        out['f_contact'] = Struct(name='output_data', mode='vertex',
                                  data=fall)
        stress = pb.evaluate('dw_tl_he_neohook_compressible.i.Omega0(ring.lm,'
                             ' v, u)', mode='el_avg', term_mode='stress',
                             verbose=False)
        aux = nm.zeros((pb.domain.mesh.n_el,) + stress.shape[1:],
                       dtype=nm.float64)
        aux[pb.domain.regions['Omega0'].cells] = stress
        out['pk2_stress'] = Struct(name='output_data', mode='cell', data=aux)
        return out

    def post_process_final(pb, state=None):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        hh = nm.array(hist.rows)
        t, fx, fy, ek, ei, angle, cx, cy = hh.T
        e0 = ek[0] + ei[0]
        output('ring on a wall (mu = %g, rho = %g):' % (mu, rho))
        output('  kinetic energy: min %.4f at t = %.3f, final %.4f of the'
               ' initial' % (ek.min() / e0, t[nm.argmin(ek)], ek[-1] / e0))
        output('  internal energy: max %.4f at t = %.3f, final %.4f'
               % (ei.max() / e0, t[nm.argmax(ei)], ei[-1] / e0))
        output('  total energy: final %.4f, max increase %.3e'
               % ((ek[-1] + ei[-1]) / e0,
                  max(0.0, ((ek + ei) / e0 - 1.0).max())))
        output('  final rotation %.2f deg, centre of mass (%.3f, %.3f)'
               % (nm.rad2deg(angle[-1]), cx[-1], cy[-1]))
        contact = nm.where(nm.abs(fy) > 0)[0]
        if len(contact):
            output('  contact from t = %.3f to t = %.3f, peak normal force'
                   ' %.4e' % (t[contact[0]], t[contact[-1]],
                              nm.abs(fy).max()))

        fig, axs = plt.subplots(1, 3, figsize=(14, 4))
        ax = axs[0]
        ax.plot(t, (ek + ei) / e0, 'k-', label='total')
        ax.plot(t, ei / e0, 'b:', label='internal')
        ax.plot(t, ek / e0, 'r--', label='kinetic')
        ax.set_ylim(0, 2)
        ax.set_ylabel('energy / initial energy')
        ax.legend()
        ax = axs[1]
        ax.plot(t, -fy, 'k-', label='normal')
        ax.plot(t, -fx, 'g-', label='tangential')
        ax.set_ylabel('contact force on the ring')
        ax.legend()
        ax = axs[2]
        ax.plot(t, nm.rad2deg(angle), 'm-')
        ax.set_ylabel('rotation of the ring [deg]')
        for ax in axs:
            ax.set_xlabel('time')
            ax.grid(True)
        fig.suptitle('ring on a wall, mu = %g' % mu)
        fig.tight_layout()
        fig.savefig(inodir('history.png'), dpi=150)
        plt.close(fig)

        nm.savetxt(inodir('history.txt'), hh,
                   header='t F_x F_y (on the wall) e_kin e_int rotation'
                   ' com_x com_y')

    solvers = {
        'lsd' : ('ls.scipy_direct', {
            'use_presolve' : True,
        }),
        'newton' : ('nls.newton', {
            'i_max' : 1,
            'eps_a' : 1e-6,
            'eps_r' : 1.0,
        }),
        'tscd' : ('ts.central_difference_pc', {
            'contact_terms' : [contact_name],
            'contact_mass_terms' : ['dw_contact_mass'],
            'correction' : correction,
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
        'ls' : 'lsd',
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
