r"""
Impact of two hyperelastic rings inside a third, fixed one, with friction and
self contact (the two-dimensional impact problem of [1]_, Sec. 6.1): explicit
dynamics in the total Lagrangian formulation with the bipenalty
predictor-corrector scheme and the contactino library, in 2D plane stress.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} \rho \ul{v} \cdot \ddot{\ul{u}}
    + \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u}; \ul{v})
    + \int_{\Gamma_{c}} \left(\varepsilon_N \langle g_N(\ul{u}) \rangle \ul{n}
    + \ul{t}_T(\ul{u})\right) \cdot \ul{v}
    = 0
    \;, \quad \forall \ul{v} \;,

where :math:`S_{ij}` is the second Piola-Kirchhoff stress of the compressible
neo-Hookean material (``dw_tl_he_neohook_compressible``).

The largest ring (diameters 26 and 30, :math:`E = 288000`, :math:`\rho = 1`)
is fixed at its outer boundary. The smallest ring (diameters 8 and 10,
:math:`E = 10000`, :math:`\rho = 0.1`) moves with the initial velocity
:math:`(30, -30)` towards the medium sized ring (diameters 10 and 12,
:math:`E = 2250`, :math:`\rho = 0.01`); both move together towards the inner
surface of the largest ring and the medium sized ring folds and contacts
itself. Poisson's ratio is 0.125, the friction coefficient 0.3.

All the ring surfaces but the fixed outer one form a single contact surface:
the two-half-pass algorithm of contactino has no master / slave labels, so
contact between the rings and self contact of a ring are treated alike, with
the ``max_depth`` option below the ring thickness (otherwise a point of the
inner surface of a thin ring would see the outer surface of the same ring as
a penetration of the depth of the thickness). The stepper
``ts.central_difference_pc`` predicts every step without contact and the
contactino library corrects the predicted configuration (the exact limit of
the bipenalty predictor-corrector scheme, with the Coulomb friction).

Notes on the setting of [1]_: the centres given in the text, :math:`[-7.9,
8.5]` for the medium ring and :math:`[7.9, -8.5]` for the smallest one, do not
fit into the largest ring with the given diameters, and Figure 15 of [1]_ shows
the smallest ring in the upper left and the medium one in the lower right at
half these distances. The rings are placed as in the figure, i.e. the
smallest ring at :math:`(-3.95, 4.25)` and the medium one at :math:`(3.95,
-4.25)`, and the initial velocity :math:`(30, -30)` points from the smallest
ring to the medium one, as the description and Figure 16 require. The paper
uses the implicit Newmark method with :math:`\Delta t = 0.001`.

.. [1] B. Yang, T. A. Laursen, A large deformation mortar formulation of self
       contact with finite sliding, Comput. Methods Appl. Mech. Engrg. 197
       (2008) 756-772.

Usage examples::

  sfepy-run sfepy/examples/large_deformation/three_rings_self_contact.py
  sfepy-view output/three_rings/rings.*.vtk -2 -f u:wu:f1:p0 1:vw:p0
"""
import os.path as op
import sys
from functools import partial
import inspect

import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.ioutils import markdown_table_from_dict
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mechanics.matcoefs import lame_from_youngpoisson
from sfepy.mesh.mesh_generators import gen_cylinder_mesh

def gen_ring_mesh(r1, r2, centre, n_r, n_phi, mat_id):
    """
    An annulus with the radii `r1`, `r2` around `centre`, with `n_r` and
    `n_phi` elements across and along it.
    """
    mesh = gen_cylinder_mesh([r1, r1, r2, r2, 0.0], [n_r + 1, n_phi, 1],
                             [0.0, 0.0, 0.0], axis='z', make_2d=True,
                             name='ring', verbose=False)
    mesh.coors[:] += centre
    mesh.cmesh.cell_groups[:] = mat_id
    return mesh

def get_matrices(mesh, conf_dict):
    """
    The tangent stiffness matrix in the reference configuration and the
    lumped mass matrix, for the critical time step.
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
    pb.get_initial_state()
    mtx_k = pb.evaluate(conf.stiffness, mode='weak', dw_mode='matrix',
                        verbose=False)
    mtx_m = pb.evaluate(conf.mass, mode='weak', dw_mode='matrix',
                        verbose=False)
    return mtx_k, mtx_m

def define(
        mu=0.3,
        v0=(30.0, -30.0),

        rings=((13.0, 15.0, (0.0, 0.0), 288000.0, 1.0, (4, 120)),
               (4.0, 5.0, (-3.95, 4.25), 10000.0, 0.1, (4, 112)),
               (5.0, 6.0, (3.95, -4.25), 2250.0, 0.01, (4, 140))),
        poisson=0.125,

        eps_n=1e6,
        eps_t=None,
        algorithm='two_half_pass,max_depth=0.5',

        dt_factor=0.5,
        t1=0.3,
        save_times=151,

        output_dir='output/three_rings',
        verbose=True,
):
    """
    Parameters
    ----------
    mu : float
        The Coulomb friction coefficient.
    v0 : tuple of 2 floats
        The initial velocity of the smallest ring.
    rings : tuple
        For the largest (fixed), the smallest (moving) and the medium sized
        ring: the inner and outer radii, the centre, Young's modulus, the
        density and the numbers of elements across and along the ring.
    poisson : float
        Poisson's ratio of all rings.
    eps_n, eps_t : float
        The normal and tangential penalty parameters (eps_t defaults to
        eps_n); the predictor-corrector correction does not depend on them.
    algorithm : str
        The contact algorithm and options of contactino, see ``dw_contact``.
    dt_factor : float
        The time step in units of the critical time step in the reference
        configuration. The neo-Hookean material stiffens in compression, so
        the critical time step decreases in the folds of the medium sized
        ring; with 0.9, an element there inverts at t = 0.21.
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

    mesh = None
    for ir, (r1, r2, centre, young, rho, (n_r, n_phi)) in enumerate(rings):
        ring = gen_ring_mesh(r1, r2, centre, n_r, n_phi, ir)
        mesh = ring if mesh is None else mesh + ring
    mesh.name = 'rings'

    def mesh_hook(mesh_, mode):
        if mode == 'read':
            return mesh

    filename_mesh = UserMeshIO(mesh_hook)

    r_fixed = rings[0][1]
    def get_fixed(coors, domain=None):
        rr = nm.linalg.norm(coors - rings[0][2], axis=1)
        return nm.where(nm.abs(rr - r_fixed) < 1e-6 * r_fixed)[0]

    regions = {
        'Omega' : 'all',
        'Surface' : ('vertices of surface', 'facet'),
        'Fixed' : ('vertices by get_fixed', 'facet'),
        'Contact' : ('r.Surface -v r.Fixed', 'facet'),
    }
    for ir in range(len(rings)):
        regions['Omega%d' % ir] = 'cells of group %d' % ir

    fields = {
        'displacement' : ('real', 'vector', 'Omega', 1),
    }

    # With auto_transform_equations=True, the solver creates du, ddu, dv, ddv.
    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    solid = {'.lumping' : 'row_sum', '.beta' : 1.0}
    for ir, (r1, r2, centre, young, rho, _) in enumerate(rings):
        solid['lm%d' % ir] = nm.array([lame_from_youngpoisson(young, poisson)])
        solid['rho%d' % ir] = rho

    materials = {
        'solid' : (solid,),
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

    stiffness = ' + '.join('dw_tl_he_neohook_compressible.i.Omega%d'
                           '(solid.lm%d, v, u)' % (ir, ir)
                           for ir in range(len(rings)))
    mass = ' + '.join('de_mass.i.Omega%d(solid.rho%d, solid.lumping,'
                      ' solid.beta, v, u)' % (ir, ir)
                      for ir in range(len(rings)))

    # The critical time step in the reference configuration.
    from sfepy.solvers.ts_solvers import get_max_eigenfrequency
    mtx_k, mtx_m = get_matrices(mesh, dict(
        regions={key : val for key, val in regions.items()
                 if key.startswith('Omega')},
        fields=fields, variables=variables, materials={'solid' : (solid,)},
        integrals=integrals, stiffness=stiffness, mass=mass,
    ))
    omega_max = get_max_eigenfrequency(mtx_k, mtx_m)
    dt_crit = 2.0 / omega_max
    dt = dt_factor * dt_crit

    signature = inspect.signature(define)
    arg_names = list(signature.parameters.keys())
    options = {name : args[name] for name in arg_names}
    options.update(omega_max=omega_max, dt_crit=dt_crit, dt=dt,
                   n_step=int(nm.ceil(t1 / dt)), n_el=mesh.n_el)
    with open(inodir('options.md'), 'w') as fd:
        fd.write(markdown_table_from_dict(options))

    output('three rings: %d cells, omega_max = %.4e, dt_crit = %.4e,'
           ' dt = %.4e, %d steps' % (mesh.n_el, omega_max, dt_crit, dt,
                                     int(nm.ceil(t1 / dt))))

    def get_ic_u(coors, ic):
        return nm.zeros_like(coors)

    def get_ic_du(coors, ic):
        """
        The smallest ring moves with v0, the others are at rest.
        """
        val = nm.zeros_like(coors)
        r1, r2, centre = rings[1][:3]
        rr = nm.linalg.norm(coors - centre, axis=1)
        val[rr < r2 + 1e-6] = v0
        return val

    functions = {
        'get_fixed' : (get_fixed,),
        'get_ic_u' : (get_ic_u,),
        'get_ic_du' : (get_ic_du,),
    }

    ics = {
        'ic' : ('Omega', {'u.all' : 'get_ic_u', 'du.all' : 'get_ic_du'}),
    }

    # The fixed outer boundary: the velocities and accelerations as well, so
    # that it is immovable in the contact correction.
    ebcs = {
        'fix' : ('Fixed', {'u.all' : 0.0, 'du.all' : 0.0, 'ddu.all' : 0.0}),
    }

    contact_term = ('dw_contact_friction.i.Contact(contact.algorithm,'
                    ' contact.epss, contact.epst, contact.mu, v, u)')

    # The mass terms must come first and have ',v,' in the middle of their
    # arguments: transform_equations_ed() replaces ',v,' by ',ddv,'.
    equations = {
        'balance' : (mass.replace(', v, u)', ', v, ddu)') + ' + '
                     + stiffness + ' + ' + contact_term + ' = 0'),
    }

    # History: time, kinetic and internal energies (the work of the internal
    # forces), the centres of mass of the moving rings.
    hist = Struct(name='history', rows=[], mtx_m=None, ei=0.0, u_last=None,
                  f_last=None)

    def step_hook(pb, ts, variables):
        u = variables['u']().copy()
        du = variables['du']()

        if hist.mtx_m is None:
            hist.mtx_m = pb.evaluate(mass.replace(', v, u)', ', ddv, ddu)'),
                                     mode='weak', dw_mode='matrix',
                                     active_only=False,
                                     verbose=False).tocsr()
            hist.groups = [pb.domain.regions['Omega%d' % ir].vertices
                           for ir in range(len(rings))]
            hist.mi = hist.mtx_m.diagonal().reshape((-1, dim))[:, 0]

        fint = pb.evaluate(stiffness, mode='weak', dw_mode='vector',
                           active_only=False, verbose=False)
        if hist.u_last is not None:
            hist.ei += 0.5 * (fint + hist.f_last) @ (u - hist.u_last)
        hist.u_last, hist.f_last = u, fint

        xx = pb.domain.mesh.coors + u.reshape((-1, dim))
        coms = []
        for ii in hist.groups[1:]:
            mi = hist.mi[ii]
            coms.extend((mi[:, None] * xx[ii]).sum(0) / mi.sum())

        ek = 0.5 * du @ (hist.mtx_m @ du)
        hist.rows.append((ts.time, ek, hist.ei, *coms))

    def post_process(out, pb, state, extend=False):
        fall = pb.get_solver().contact_force
        n_nod = pb.domain.mesh.n_nod
        fall = (nm.zeros((n_nod, dim)) if fall is None
                else fall.reshape((-1, dim)))
        out['f_contact'] = Struct(name='output_data', mode='vertex',
                                  data=fall)
        return out

    def post_process_final(pb, state=None):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        hh = nm.array(hist.rows)
        t, ek, ei = hh[:, :3].T
        e0 = ek[0] + ei[0]
        output('three rings (mu = %g, %s):' % (mu, algorithm))
        output('  kinetic energy: min %.4f, final %.4f of the initial'
               % (ek.min() / e0, ek[-1] / e0))
        output('  total energy: final %.4f, max increase %.3e'
               % ((ek[-1] + ei[-1]) / e0,
                  max(0.0, ((ek + ei) / e0 - 1.0).max())))

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(t, (ek + ei) / e0, 'k-', label='total')
        ax.plot(t, ei / e0, 'b:', label='internal')
        ax.plot(t, ek / e0, 'r--', label='kinetic')
        ax.set_xlabel('time')
        ax.set_ylabel('energy / initial energy')
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        fig.savefig(inodir('history.png'), dpi=150)
        plt.close(fig)

        nm.savetxt(inodir('history.txt'), hh,
                   header='t e_kin e_int com_x com_y (smallest ring)'
                   ' com_x com_y (medium ring)')

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
            'contact_terms' : ['dw_contact_friction'],
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
