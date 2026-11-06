r"""
Impact of two elastic bars: explicit dynamics with the penalty or the
bipenalty contact of the contactino library.

Find :math:`\ul{u}` such that:

.. math::
    \int_{\Omega} \rho\ \ul{v} \cdot \ddot{\ul{u}}
    + \int_{\Omega} D_{ijkl}\ e_{ij}(\ul{v}) e_{kl}(\ul{u})
    + \int_{\Gamma_c} \varepsilon_N \langle g_N(\ul{u}) \rangle \ul{n} \ul{v}
    + \int_{\Gamma_c} \varepsilon_M \ddot{g}_N(\ul{u}) \ul{n} \ul{v}
    = 0 \;, \quad \forall \ul{v} \;,

where the last term (the inertia penalty of the bipenalty method [1]_, [2]_)
is present only for ``method='bipenalty'``. The stiffness penalty alone
raises the maximum eigenfrequency of the discretization and hence lowers the
critical time step of the central difference scheme roughly as
:math:`1 / \sqrt{\varepsilon_N}`. With the inertia penalty on the same
constraint, the maximum eigenfrequency is bounded by :math:`\max(\omega_{max},
\sqrt{\varepsilon_N / \varepsilon_M})`, where :math:`\omega_{max}` is the
maximum eigenfrequency of the unpenalized discretization, so the penalty
ratio :math:`\varepsilon_N / \varepsilon_M \leq \omega_{max}^2` keeps the
critical time step while :math:`\varepsilon_N` can be arbitrarily large.

Two identical bars with the velocities :math:`\pm v_0` impact at :math:`t =
0`. For :math:`\nu = 0` (default) the bars behave as one-dimensional and the
Saint-Venant solution holds: the contact stress is :math:`\sigma_c = \rho c
v_0` with :math:`c = \sqrt{E / \rho}` for :math:`0 < t < 2L/c`, then the bars
separate with the velocities reversed. The example records the contact
force, the penetration, the momentum and the energies in every time step and
compares with this solution.

Notes:

- The problem is linear, so the results in the units of :math:`\rho c v_0`
  and of the element size do not depend on the scale of the geometry.
- The stiffness penalty :math:`\varepsilon_N = E / h` (``eps_s=1``) already
  reduces the critical time step below the unpenalized one: the default
  ``dt_factor=0.5`` is stable, 0.9 is not. The bipenalty is stable with 0.9
  for any ``eps_s``, but the contact stress oscillation grows with
  ``eps_s`` (the inertia penalty at the critical ratio adds the mass
  :math:`\varepsilon_N / \omega_{max}^2` to the interface), so moderate
  values are preferable.

Usage examples::

  # Stiffness penalty (default) with dt = 0.5 of the unpenalized critical
  # time step.
  sfepy-run sfepy/examples/linear_elasticity/explicit_dynamic_contact.py

  # Bipenalty with the same penalty and dt = 0.9 of the critical time step,
  # unstable for the stiffness penalty alone.
  sfepy-run sfepy/examples/linear_elasticity/explicit_dynamic_contact.py -d "method='bipenalty', dt_factor=0.9"
  sfepy-run sfepy/examples/linear_elasticity/explicit_dynamic_contact.py -d "method='penalty', dt_factor=0.9"

  # Bipenalty with a 100 times stiffer penalty: stable, but the contact
  # stress oscillates strongly.
  sfepy-run sfepy/examples/linear_elasticity/explicit_dynamic_contact.py -d "method='bipenalty', eps_s=100.0"

  # Bipenalty with the predictor-corrector scheme and the inertia penalty
  # from the criterion eps_M = eps_N dt^2 - m/2: the predicted penetration
  # is removed in the correction for any (large enough) eps_N.
  sfepy-run sfepy/examples/linear_elasticity/explicit_dynamic_contact.py -d "method='bipenalty_pc', eps_s=100.0, dt_factor=0.9"
  sfepy-run sfepy/examples/linear_elasticity/explicit_dynamic_contact.py -d "method='bipenalty_pc', eps_s=100.0, dt_factor=0.9, fast_rmm=True"

  # The same with the correction phase solved by the contactino library
  # (exact removal of the predicted penetration, inverse mass from the host).
  sfepy-run sfepy/examples/linear_elasticity/explicit_dynamic_contact.py -d "method='bipenalty_pc', eps_s=100.0, dt_factor=0.9, fast_rmm=True, correction='terms'"

  sfepy-view output/bars/two_bars.*.vtk -f u:wu:f1e2:p0 1:vw:p0

.. [1] H. Askes, M. Caramés-Saddler, A. Rodríguez-Ferran, Bipenalty method
       for time domain computational dynamics, Proc. R. Soc. A 466 (2010).
.. [2] J. Kopačka, A. Tkachuk, D. Gabriel, R. Kolman, M. Bischoff, J. Plešek,
       On stability and reflection-transmission analysis of the bipenalty
       method in contact-impact problems: A one-dimensional, homogeneous case
       study, Int. J. Numer. Meth. Engng 113 (2018).
"""
import os.path as op
from functools import partial
import inspect

import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.base.ioutils import markdown_table_from_dict
import sfepy.mechanics.matcoefs as mc
from sfepy.discrete.fem.meshio import UserMeshIO

def gen_two_bars(L, a, shape, gap0):
    """
    Two bars along x, bar 1 in [0, L], bar 2 in [L + gap0, 2L + gap0].
    """
    from sfepy.discrete.fem import Mesh
    from sfepy.mesh.mesh_generators import gen_block_mesh

    dims = nm.array([L, a, a])
    c1 = nm.array([L / 2, 0, 0])
    c2 = nm.array([L + gap0 + L / 2, 0, 0])

    m1 = gen_block_mesh(dims, shape, c1, verbose=False)
    m2 = gen_block_mesh(dims, shape, c2, verbose=False)

    desc = m1.descs[0]
    coors = nm.concatenate((m1.coors, m2.coors), axis=0)
    conn = nm.concatenate((m1.get_conn(desc),
                           m2.get_conn(desc) + m1.n_nod), axis=0)
    ngroups = nm.zeros(coors.shape[0], dtype=nm.int32)
    ngroups[m1.n_nod:] = 1
    mat_id = nm.zeros(conn.shape[0], dtype=nm.int32)
    mat_id[m1.n_el:] = 1

    return Mesh.from_data('two_bars', coors, ngroups, [conn], [mat_id],
                          m1.descs)

def define(
        method='penalty',
        eps_s=1.0,
        penalty_ratio=1.0,
        criterion='kwon2',
        correction='host',
        algorithm='two_half_pass',
        omega_max=None,

        young=210e9,
        poisson=0.0,
        rho=7800.0,

        L=0.2,
        a=0.02,
        shape=(41, 3, 3),
        gap0=1e-6,
        v0=5.0,

        dt_factor=0.5,
        dt=None,
        t_end_factor=1.25,
        mass_lumping='row_sum',
        mass_beta=1.0,
        fast_rmm=False,
        save_times=20,

        output_dir='output/bars',
        verbose=True,
):
    """
    Parameters
    ----------
    method : 'penalty', 'bipenalty' or 'bipenalty_pc'
        The stiffness penalty, the traditional bipenalty with the central
        difference scheme, or the bipenalty with the predictor-corrector
        scheme (ts.central_difference_pc).
    eps_s : float
        The stiffness penalty in units of E / h, h = the element size along
        the bars.
    penalty_ratio : float
        The bipenalty ratio eps_N / eps_M in units of omega_max**2
        (``criterion='ratio'`` or the traditional bipenalty).
    criterion : 'kwon2', 'kwon1' or 'ratio'
        The inertia penalty of the predictor-corrector scheme:
        eps_M = eps_N dt**2 - m_eff, where the effective surface mass m_eff is
        half of the surface node mass per unit area for two deformable
        bodies ('kwon2') or the whole one for a rigid counterpart ('kwon1'),
        or the ratio criterion as for the traditional bipenalty.
    correction : 'host' or 'terms'
        The corrector of the predictor-corrector scheme: assembled by the
        host from the contact force and the inertia penalty mass with the
        scalar criterion above, or provided by the contactino library as the
        exact removal of the predicted penetration distributed by the inverse
        mass of the host (the forward-increment Lagrange multiplier
        correction, independent of eps_N and of the time step).
    algorithm : 'two_half_pass' or 'two_pass'
        The contact algorithm of contactino, see ``dw_contact``.
    omega_max : float, optional
        The maximum eigenfrequency of the unpenalized discretization. If not
        given, the one-dimensional estimate 2 c / h for the lumped mass
        matrix (beta = 1) is used (a lower bound for the 3D hexahedra).
    dt_factor : float
        The time step as a fraction of the unpenalized critical time step
        2 / omega_max. Ignored if `dt` is given.
    """
    args = locals()

    inodir = partial(op.join, output_dir)

    output.set_output(filename=inodir('output_log.txt'), quiet=not verbose,
                      combined=True)

    if method not in ('penalty', 'bipenalty', 'bipenalty_pc'):
        raise ValueError("method must be 'penalty', 'bipenalty' or"
                         " 'bipenalty_pc'! (%s)" % method)

    dim = 3
    lam, mu = mc.lame_from_youngpoisson(young, poisson)
    c_L = nm.sqrt((lam + 2.0 * mu) / rho)
    c_0 = nm.sqrt(young / rho)
    h = L / (shape[0] - 1)
    area = a**2

    if omega_max is None:
        omega_max = 2.0 * c_L / h
    dt_crit = 2.0 / omega_max
    if dt is None:
        dt = dt_factor * dt_crit
    t1 = t_end_factor * 2.0 * L / c_0

    eps_n = eps_s * young / h
    # Surface node mass per unit area of the lumped hexahedra.
    m_surf = rho * h / 2.0
    if (method == 'bipenalty_pc') and (criterion != 'ratio'):
        m_eff = m_surf if criterion == 'kwon1' else 0.5 * m_surf
        eps_m = eps_n * dt**2 - m_eff
        if eps_m < 0.0:
            raise ValueError('the stiffness penalty is too small for the'
                             ' criterion %s: eps_s >= %.3e is needed!'
                             % (criterion, m_eff / dt**2 * h / young))

    else:
        eps_m = eps_n / (penalty_ratio * omega_max**2)
    sigma_sv = rho * c_0 * v0
    args.update(omega_max=omega_max, dt=dt)

    signature = inspect.signature(define)
    arg_names = list(signature.parameters.keys())
    options = {name : args[name] for name in arg_names}
    options.update(eps_N=eps_n, eps_M=eps_m, dt_crit=dt_crit, t1=t1,
                   h=h, c_L=c_L, c_0=c_0, sigma_saint_venant=sigma_sv)
    with open(inodir('options.md'), 'w') as fd:
        fd.write(markdown_table_from_dict(options))

    output('bars: h = %.4e, c_L = %.4e, c_0 = %.4e, omega_max = %.4e,'
           ' dt_crit = %.4e, dt = %.4e (%.3f dt_crit), %d steps'
           % (h, c_L, c_0, omega_max, dt_crit, dt, dt / dt_crit,
              int(nm.ceil(t1 / dt))))
    output('penalties: eps_N = %.4e (%.2f E/h), eps_M = %.4e'
           ' (%.3f surface node masses), eps_N/eps_M = %.4e'
           ' = %.3f omega_max^2'
           % (eps_n, eps_s, eps_m, eps_m / m_surf,
              eps_n / max(eps_m, 1e-300),
              eps_n / max(eps_m, 1e-300) / omega_max**2))

    def mesh_hook(mesh, mode):
        if mode == 'read':
            return gen_two_bars(L, a, shape, gap0)

        elif mode == 'write':
            pass

    filename_mesh = UserMeshIO(mesh_hook)

    eps = 1e-8 * L
    regions = {
        'Omega' : 'all',
        'Bar1' : 'cells of group 0',
        'Bar2' : 'cells of group 1',
        'Contact1' : ('vertices in (x > %.12e) & (x < %.12e)'
                      % (L - eps, L + eps), 'facet', 'Bar1'),
        'Contact2' : ('vertices in (x > %.12e) & (x < %.12e)'
                      % (L + gap0 - eps, L + gap0 + eps), 'facet', 'Bar2'),
        'Contact' : ('r.Contact1 +s r.Contact2', 'facet'),
    }

    materials = {
        'solid' : ({
            'D' : mc.stiffness_from_youngpoisson(dim, young, poisson),
            'rho' : rho,
            '.lumping' : mass_lumping,
            '.beta' : mass_beta,
        },),
        'contact' : ({
            '.epss' : eps_n,
            '.epsm' : eps_m,
            '.algorithm' : algorithm,
        },),
    }

    fields = {
        'displacement' : ('real', 'vector', 'Omega', 1),
    }

    # With auto_transform_equations=True, the solver creates du, ddu, dv, ddv.
    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    def get_ic_u(coors, ic):
        return nm.zeros_like(coors)

    def get_ic_du(coors, ic):
        """Bar 1 moves with +v0, bar 2 with -v0 along x."""
        val = nm.zeros_like(coors)
        mid = L + gap0 / 2
        val[coors[:, 0] < mid, 0] = v0
        val[coors[:, 0] >= mid, 0] = -v0
        return val

    functions = {
        'get_ic_u' : (get_ic_u,),
        'get_ic_du' : (get_ic_du,),
    }

    ics = {
        'ic' : ('Omega', {'u.all' : 'get_ic_u', 'du.all' : 'get_ic_du'}),
    }

    ebcs = {}

    integrals = {
        'i' : 2,
    }

    # The mass term must come first and have ',v,' in the middle of its
    # arguments: transform_equations_ed() replaces ',v,' by ',ddv,'.
    mass_term = 'de_mass.i.Omega(solid.rho, solid.lumping, solid.beta, v, ddu)'
    contact_terms = '+ dw_contact.i.Contact(contact.algorithm, contact.epss, v, u)'
    if method in ('bipenalty', 'bipenalty_pc'):
        contact_terms += ('\n + dw_contact_mass.i.Contact(contact.algorithm,'
                          ' contact.epsm, v, ddu)')

    equations = {
        'balance' : mass_term + """
         + dw_lin_elastic.i.Omega(solid.D, v, u)
        """ + contact_terms + ' = 0',
    }

    # History of the contact force, penetration, momentum and energies.
    hist = Struct(name='history', t=[], f1=[], f2=[], gap=[], p=[], ek=[],
                  es=[], mtx_m=None, mtx_k=None, pairs=None)

    def get_contact_term(pb):
        return [tt for tt in pb.equations['K'].terms
                if tt.name == 'dw_contact'][0]

    def step_hook(pb, ts, variables):
        import scipy.sparse as sps

        term = get_contact_term(pb)
        ci = term.ci
        uvar, duvar = variables['u'], variables['du']
        u = uvar()
        du = duvar()

        if hist.mtx_m is None:
            # Constant matrices for the momentum and the energies.
            hist.mtx_m = pb.evaluate(
                'de_mass.i.Omega(solid.rho, solid.lumping, solid.beta,'
                ' ddv, ddu)', mode='weak', dw_mode='matrix',
                verbose=False).tocsr()
            hist.mtx_k = pb.evaluate(
                'dw_lin_elastic.i.Omega(solid.D, v, u)', mode='weak',
                dw_mode='matrix', verbose=False).tocsr()
            # Pairs of facing nodes of the conforming contact surfaces.
            coors = pb.domain.mesh.coors
            r1 = pb.domain.regions['Contact1'].vertices
            r2 = pb.domain.regions['Contact2'].vertices
            i1 = r1[nm.lexsort((coors[r1, 2], coors[r1, 1]))]
            i2 = r2[nm.lexsort((coors[r2, 2], coors[r2, 1]))]
            hist.pairs = (i1, i2)

        # Contact forces acting on the bars (x components).
        body = nm.isin(ci.sd.nodes, pb.domain.regions['Bar2'].vertices)
        if method == 'bipenalty_pc':
            # The total contact force M da of the last correction.
            tss = pb.get_solver()
            fall = tss.contact_force
            if fall is None:
                fall = nm.zeros(u.shape[0])
            fc = fall.reshape((-1, dim))[ci.sd.nodes]

        else:
            # The residual contribution is the internal force.
            uc = ci.get_surface_displacements(uvar)
            fc = -ci.cd.eval_residual(uc).reshape((-1, dim))

        f1 = fc[~body, 0].sum()
        f2 = fc[body, 0].sum()

        uu = u.reshape((-1, dim))
        i1, i2 = hist.pairs
        gap = ((coors_x := pb.domain.mesh.coors[:, 0])[i2] + uu[i2, 0]
               - coors_x[i1] - uu[i1, 0]).min()

        mv = hist.mtx_m @ du
        hist.t.append(ts.time)
        hist.f1.append(f1)
        hist.f2.append(f2)
        hist.gap.append(gap)
        hist.p.append(mv.reshape((-1, dim))[:, 0].sum())
        hist.ek.append(0.5 * du @ mv)
        hist.es.append(0.5 * u @ (hist.mtx_k @ u))

    def post_process(out, pb, state, extend=False):
        ev = pb.evaluate
        stress = ev('ev_cauchy_stress.i.Omega(solid.D, u)', mode='el_avg',
                    copy_materials=False, verbose=False)
        out['cauchy_stress'] = Struct(name='output_data', mode='cell',
                                      data=stress)
        return out

    def post_process_final(pb, state=None):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        t = nm.array(hist.t)
        f1, f2 = nm.array(hist.f1), nm.array(hist.f2)
        gap = nm.array(hist.gap)
        p, ek, es = nm.array(hist.p), nm.array(hist.ek), nm.array(hist.es)

        tc = 2.0 * L / c_0
        sig1, sig2 = -f1 / area, f2 / area
        ii = (t > 0.1 * tc) & (t < 0.9 * tc)
        if not ii.any():
            ii = slice(None)
        e0 = ek[0]
        output('two bars impact (%s): eps_N = %.3e, eps_M = %.3e, dt/dt_crit'
               ' = %.3f' % (method, eps_n, eps_m, dt / dt_crit))
        output('  Saint-Venant contact stress rho c v0 = %.6e' % sigma_sv)
        output('  mean contact stress in (0.1, 0.9) 2L/c: bar 1 %.6e (%+.2f %%),'
               ' bar 2 %.6e (%+.2f %%)'
               % (sig1[ii].mean(), 100 * (sig1[ii].mean() / sigma_sv - 1),
                  sig2[ii].mean(), 100 * (sig2[ii].mean() / sigma_sv - 1)))
        output('  max |oscillation| of the contact stress: %.2f %%'
               % (100 * nm.abs(sig1[ii] - sigma_sv).max() / sigma_sv))
        output('  max penetration: %.4e (%.3f %% of h, eps_N estimate %.4e)'
               % (-gap.min(), -100 * gap.min() / h, sigma_sv / eps_n))
        output('  momentum x: initial %.6e, max |change| %.3e'
               % (p[0], nm.abs(p - p[0]).max()))
        output('  energy: initial %.6e, max |(ek + es) - e0| / e0 = %.3e,'
               ' final %.3e'
               % (e0, nm.abs(ek + es - e0).max() / e0, (ek[-1] + es[-1]) / e0 - 1))

        fig, axs = plt.subplots(2, 2, figsize=(10, 7))
        ax = axs[0, 0]
        ax.plot(t / tc, sig1 / sigma_sv, label='bar 1')
        ax.plot(t / tc, sig2 / sigma_sv, label='bar 2', ls='--')
        ax.axhline(1.0, color='k', lw=0.8, label=r'$\rho c v_0$')
        ax.set_xlabel('t / (2L/c)')
        ax.set_ylabel(r'contact stress / $\rho c v_0$')
        ax.legend()
        ax = axs[0, 1]
        ax.plot(t / tc, -gap / h)
        ax.set_xlabel('t / (2L/c)')
        ax.set_ylabel('penetration / h')
        ax = axs[1, 0]
        ax.plot(t / tc, (p - p[0]) / (rho * area * L * v0))
        ax.set_xlabel('t / (2L/c)')
        ax.set_ylabel(r'$\Delta$ momentum / ($\rho A L v_0$)')
        ax = axs[1, 1]
        ax.plot(t / tc, ek / e0, label='kinetic')
        ax.plot(t / tc, es / e0, label='strain')
        ax.plot(t / tc, (ek + es) / e0, label='total', color='k')
        ax.set_xlabel('t / (2L/c)')
        ax.set_ylabel('energy / E_0')
        ax.legend()
        for ax in axs.ravel():
            ax.grid(True)
        if (method == 'bipenalty_pc') and (correction == 'terms'):
            mass_info = 'correction by contactino'
        elif (method == 'bipenalty_pc') and (criterion != 'ratio'):
            mass_info = 'eps_M = eps_N dt^2 - m_eff (%s)' % criterion
        else:
            mass_info = 'eps_N/eps_M = %.2g omega_max^2' % penalty_ratio
        fig.suptitle('%s: eps_N = %.2g E/h, %s, dt = %.2f dt_crit'
                     % (method, eps_s, mass_info, dt / dt_crit))
        fig.tight_layout()
        fig.savefig(inodir('history.png'), dpi=150)
        plt.close(fig)

        nm.savetxt(inodir('history.txt'), nm.c_[t, sig1, sig2, gap, p, ek, es],
                   header='t sigma_c1 sigma_c2 gap momentum_x e_kin e_strain')

    if fast_rmm:
        if method == 'bipenalty':
            raise ValueError('the reciprocal mass matrix solver does not'
                             ' support the traditional bipenalty mass yet,'
                             " use method='bipenalty_pc'!")
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
            'contact_terms' : ['dw_contact'],
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
