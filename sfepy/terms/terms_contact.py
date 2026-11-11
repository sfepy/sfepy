import numpy as nm

from sfepy.base.base import output, Struct
from sfepy.terms.terms import Term
from sfepy.terms.extmods import terms

try:
    import contactino as cc

except ImportError:
    # The optional contactino package is not installed: the contact terms
    # cannot be evaluated.
    cc = None

class ContactInfo(Struct):
    """
    Contact-related data of contact terms: the contact surface description
    handed over to contactino and the mapping of the contact surface DOFs to
    the DOFs of the state variable.
    """
    # (space dimension, nodes per facet) -> contactino element type. The
    # facet nodes of the second order fields (the vertices, then the edge
    # midpoints in the order of the facet edges, then the centre) are ordered
    # as contactino expects. The constants are looked up in the contactino
    # module only when needed, so that this module imports without it.
    element_types = {
        (2, 2) : 'CONTACTINO_LINE',
        (2, 3) : 'CONTACTINO_QUADRATIC_EDGE',
        (3, 3) : 'CONTACTINO_TRIANGLE',
        (3, 6) : 'CONTACTINO_QUADRATIC_TRIANGLE',
        (3, 4) : 'CONTACTINO_QUAD',
        (3, 8) : 'CONTACTINO_QUADRATIC_QUAD',
        (3, 9) : 'CONTACTINO_BIQUADRATIC_QUAD',
    }
    # The contact algorithm: the two-half-pass (each pass assembles the rows of
    # its slave only, unbiased, not symmetric) or the two-full-pass (each pass
    # assembles both the slave and the master rows with the weight 1/2, exact
    # action and reaction, symmetric frictionless tangent).
    algorithms = {
        'two_half_pass' : 'CONTACT_ALGORITHM_TWO_HALF_PASS',
        'two_pass' : 'CONTACT_ALGORITHM_TWO_PASS',
        'one_pass' : 'CONTACT_ALGORITHM_ONE_PASS',
    }
    # The contact discretization: the penalty on the gap of every quadrature
    # point, or on the gap L2-projected onto the slave nodes (mortar, linear
    # facets only).
    discretizations = {
        'qpts' : 'CONTACT_DISCRETIZATION_QPTS',
        'mortar' : 'CONTACT_DISCRETIZATION_MORTAR',
    }

    def __init__(self, region, integral, geo, state, epss=1.0, epst=0.0,
                 mu=0.0, algorithm='two_half_pass', verbosity=0):
        if cc is None:
            raise ModuleNotFoundError('the contact terms need the contactino'
                                      ' package (pip install contactino)!',
                                      name='contactino')

        # Uses field connectivity (higher order nodes).
        sd = state.field.extra_data[f'sd_{region.name}']

        nsd = region.dim
        key = (nsd, sd.n_fp)
        if key not in self.element_types:
            raise ValueError('unsupported contact facet: %d nodes per facet'
                             ' in %dD!' % (sd.n_fp, nsd))

        # Number of 1D Gauss points integrating polynomials of integral.order
        # exactly (contactino accepts 1 to 6).
        nqp_1d = min(max(1, (integral.order + 2) // 2), 6)

        iparm = nm.zeros(64, dtype=nm.int32)
        iparm[cc.IPARM_NSD] = nsd
        iparm[cc.IPARM_NNP] = len(sd.nodes)
        iparm[cc.IPARM_NCE] = sd.n_fa
        # mu = 0 gives frictionless contact, epst is then ignored.
        frictional = mu > 0.0
        if frictional and not (epst > 0.0):
            raise ValueError('frictional contact needs a positive tangential'
                             ' penalty epst!')
        if frictional:
            iparm[cc.IPARM_CONTACT_TYPE] = cc.CONTACT_TYPE_FRICTIONAL
            iparm[cc.IPARM_TANGENT_ENFORCEMENT] = \
                cc.TANGENT_CONTACT_ENFORCEMENT_METHOD_PENALTY

        else:
            iparm[cc.IPARM_CONTACT_TYPE] = cc.CONTACT_TYPE_FRICTIONLESS
            iparm[cc.IPARM_TANGENT_ENFORCEMENT] = \
                cc.TANGENT_CONTACT_ENFORCEMENT_METHOD_NONE
        iparm[cc.IPARM_NORMAL_ENFORCEMENT] = \
            cc.NORMAL_CONTACT_ENFORCEMENT_METHOD_PENALTY
        # '<algorithm>[:<slave region>][,<option>=<value>...]': e.g.
        # 'one_pass:<region>' is the one-pass algorithm with the facets of the
        # given subregion of the contact region as the slave ones, the option
        # max_depth bounds the depth of the penetration searched for (e.g.
        # below the thickness of thin bodies in self contact), the option
        # discretization selects the quadrature-point or the mortar
        # discretization.
        spec, *opts = [token.strip() for token in algorithm.split(',')]
        options = {'max_depth' : 0.0, 'discretization' : 'qpts'}
        for opt in opts:
            opt_name, _, opt_val = opt.partition('=')
            opt_name, opt_val = opt_name.strip(), opt_val.strip()
            if opt_name == 'max_depth':
                options[opt_name] = float(opt_val)

            elif opt_name == 'discretization':
                if opt_val not in self.discretizations:
                    raise ValueError('unknown contact discretization %s!'
                                     ' (one of %s)'
                                     % (opt_val,
                                        sorted(self.discretizations.keys())))
                options[opt_name] = opt_val

            else:
                raise ValueError('unknown contact option %s! (%s)'
                                 % (opt_name, algorithm))
        alg_name, _, slave_name = spec.partition(':')
        if alg_name not in self.algorithms:
            raise ValueError('unknown contact algorithm %s! (one of %s)'
                             % (algorithm, sorted(self.algorithms.keys())))
        if (alg_name == 'one_pass') != bool(slave_name):
            raise ValueError("the slave region has to be given as"
                             " 'one_pass:<region name>', and only for the"
                             " one-pass algorithm! (%s)" % algorithm)
        iparm[cc.IPARM_ALGORITHM] = getattr(cc, self.algorithms[alg_name])
        iparm[cc.IPARM_DISCRETIZATION] = \
            getattr(cc, self.discretizations[options['discretization']])
        iparm[cc.IPARM_NQP_1D] = nqp_1d
        iparm[cc.IPARM_VERBOSITY] = {
            0 : cc.VERBOSITY_NONE, 1 : cc.VERBOSITY_INFO,
        }.get(verbosity, cc.VERBOSITY_DEBUG)

        # The contact parameters are those of contactino: the frictional
        # traction depends on the normal pressure, so no scaling of the
        # output is possible. A change of the parameters needs a new instance
        # (see ContactTerm.get_contact_info()).
        rparm = nm.zeros(64, dtype=nm.float64)
        rparm[cc.RPARM_MAX_PENETRATION] = options['max_depth']
        rparm[cc.RPARM_EPS_N] = epss
        if frictional:
            # Contactino rejects a tangential penalty without friction.
            rparm[cc.RPARM_EPS_T] = epst
            rparm[cc.RPARM_MU] = mu

        # sd.nodes: field node indices of the contact surface nodes.
        # sd.leconn: facet connectivity in the local (surface) node numbering,
        # oriented as the reference element facets, i.e. with outward
        # normals, as contactino requires.
        coords = nm.ascontiguousarray(state.field.coors[sd.nodes].ravel(),
                                      dtype=nm.float64)
        connectivity = nm.ascontiguousarray(sd.leconn.ravel(), dtype=nm.int32)
        element_types = nm.full(sd.n_fa, getattr(cc, self.element_types[key]),
                                dtype=nm.int32)

        cd = cc.ContactinoData(iparm, rparm, coords, connectivity,
                               element_types)
        roles = None
        if alg_name == 'one_pass':
            slave = region.domain.regions[slave_name]
            is_slave = nm.isin(sd.nodes[sd.leconn], slave.vertices).all(axis=1)
            if is_slave.all() or not is_slave.any():
                raise ValueError('the slave region %s has to contain a part'
                                 ' of the contact facets!' % slave_name)
            roles = nm.where(is_slave, cc.CONTACT_ELEMENT_ROLE_SLAVE,
                             cc.CONTACT_ELEMENT_ROLE_MASTER).astype(nm.int32)
            cd.set_element_roles(roles)

        # Local contact DOF j * nsd + ic -> variable DOF sd.nodes[j] * nsd + ic.
        global_dofs = (nsd * nm.repeat(sd.nodes, nsd)
                       + nm.tile(nm.arange(nsd), len(sd.nodes)))

        Struct.__init__(self, name='contact_info_%s' % region.name,
                        sd=sd, nsd=nsd, neq=cd.neq, iparm=iparm, rparm=rparm,
                        pars=(epss, epst, mu, algorithm),
                        frictional=frictional, roles=roles,
                        coords=coords, connectivity=connectivity,
                        element_types=element_types,
                        global_dofs=global_dofs.astype(nm.int32), cd=cd)

    def get_surface_displacements(self, state):
        """
        Displacements of the contact surface nodes in the contactino
        (node-major, local surface node) ordering.
        """
        u = state().reshape((-1, self.nsd))
        return nm.ascontiguousarray(u[self.sd.nodes].ravel(), dtype=nm.float64)

class ContactTerm(Term):
    r"""
    Frictionless contact term with a penalty function, evaluated by the
    contactino library.

    The penalty function is defined as :math:`\varepsilon_N \langle g_N(\ul{u})
    \rangle`, where :math:`\varepsilon_N` is the normal penalty parameter and
    :math:`\langle g_N(\ul{u}) \rangle` are the Macaulay's brackets of the gap
    function :math:`g_N(\ul{u})`. See :class:`ContactFrictionTerm` for the
    contact with the Coulomb friction.

    This term has a dynamic connectivity of DOFs in its region.

    :Definition:

    .. math::
        \int_{\Gamma_{c}} \varepsilon_N \langle g_N(\ul{u}) \rangle \ul{n}
        \cdot \ul{v}

    :Arguments:
        - opt_material : the contact algorithm, ``'two_half_pass'`` (default),
          ``'two_pass'`` (exact action and reaction, symmetric frictionless
          tangent) or ``'one_pass:<region>'`` (the facets of the given
          subregion of the contact region are the slaves, the others the
          masters; e.g. the concave side against a convex one), optionally
          followed by ``',max_depth=<value>'``, the maximum depth of the
          penetration searched for (required for the self contact of thin
          bodies: below their thickness), and by
          ``',discretization=<qpts|mortar>'``, the penalty on the gap of
          every quadrature point (default) or on the gap L2-projected onto
          the slave nodes (linear facets only)
        - material     : :math:`\varepsilon_N`
        - virtual      : :math:`\ul{v}`
        - state        : :math:`\ul{u}`
    """
    name = 'dw_contact'
    arg_types = ('opt_material', 'material', 'virtual', 'state')
    arg_shapes = [{'opt_material' : '.: str', 'material' : '.: 1',
                   'virtual' : ('D', 'state'), 'state' : 'D'},
                  {'opt_material' : None}]
    integration = 'facet'

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.ci = None
        self._gap_warned = False

    def get_contact_info(self, geo, state, epss=1.0, epst=0.0, mu=0.0,
                         algorithm=None):
        """
        The contact info with the given contact parameters. A change of the
        parameters re-creates it, which is an error in frictional contact
        (the frictional history would be lost).
        """
        if algorithm is None:
            algorithm = 'two_half_pass'
        pars = (float(epss), float(epst), float(mu), algorithm)
        if (self.ci is not None) and (self.ci.pars != pars):
            if self.ci.frictional:
                raise ValueError('%s: the contact parameters cannot change'
                                 ' in frictional contact!' % self.name)
            self.ci = None

        if self.ci is None:
            self.ci = ContactInfo(self.region, self.integral, geo, state,
                                  *pars, verbosity=self.verbosity)

        return self.ci

    def get_nodal_forces(self, state=None):
        """
        The contact forces acting on the bodies at the contact surface nodes
        (contactino ordering) for the state variable (the term's one by
        default), after an evaluation of the term. The residual contribution
        is the internal force, i.e. minus the force on the body.

        Returns the field node indices of the contact surface nodes and the
        forces, shape ``(n_nod, dim)``.
        """
        ci = self._get_evaluated_ci()
        if state is None:
            state = self.get_args(['state'])[0]
        uc = ci.get_surface_displacements(state)
        forces = -ci.cd.eval_residual(uc).reshape((-1, ci.nsd))
        return ci.sd.nodes, forces

    def get_nodal_fields(self, state=None):
        """
        The normal gap (penetration-positive, -inf where no master was found)
        and the contact pressure at the contact surface nodes (contactino
        ordering) for the state variable (the term's one by default), after
        an evaluation of the term, see
        :meth:`ContactinoData.eval_nodal_fields()
        <contactino.ContactinoData.eval_nodal_fields>`.

        Returns the field node indices of the contact surface nodes, the gaps
        and the pressures.
        """
        ci = self._get_evaluated_ci()
        if state is None:
            state = self.get_args(['state'])[0]
        uc = ci.get_surface_displacements(state)
        gap, pressure = ci.cd.eval_nodal_fields(uc)
        return ci.sd.nodes, gap, pressure

    def _get_evaluated_ci(self):
        if self.ci is None:
            raise ValueError('%s: the term has not been evaluated yet!'
                             % self.name)
        return self.ci

    def is_frictional(self):
        """
        True for a frictional contact, after an evaluation of the term.
        """
        return (self.ci is not None) and self.ci.frictional

    def freeze_state(self):
        """
        Freeze the contact state of the most recent evaluation of the term:
        the active set, the master element of every active point and its
        stick or slip branch; the following evaluations are smooth in the
        displacements, e.g. for an active set strategy of a nonlinear solver.
        """
        self._get_evaluated_ci().cd.freeze_state()

    def release_state(self, state=None):
        """
        Release the frozen contact state and return the number of quadrature
        points whose free state at the current values of the state variable
        (the term's one by default) differs from the frozen one.
        """
        ci = self._get_evaluated_ci()
        if state is None:
            state = self.get_args(['state'])[0]
        uc = ci.get_surface_displacements(state)
        return ci.cd.release_state(uc)

    def commit_correction(self):
        """
        Commit the state of the most recent configuration correction (see
        :meth:`get_configuration_correction()`) as the reference state of the
        frictional history of the next correction. Frictionless contact does
        not need it.
        """
        if self.is_frictional():
            self.ci.cd.commit_correction()

    def commit_state(self):
        """
        Commit the contact state of the most recent evaluation of the term as
        the reference state of the frictional history: the master points and
        the frictional tractions of the active quadrature points. Frictionless
        contact does not need it.
        """
        if self.is_frictional():
            self.ci.cd.commit_step()

    def advance(self, ts):
        """
        Commit the contact state of the most recent evaluation as the
        converged state of the time step (the frictional history). The most
        recent evaluation has to be at the converged state, which holds for
        the implicit time stepping solvers. Solvers committing the state
        themselves (e.g. :class:`CentralDifferencePCTS
        <sfepy.solvers.ts_solvers.CentralDifferencePCTS>`) set the
        `commit_on_advance` attribute to False.
        """
        if getattr(self, 'commit_on_advance', True):
            self.commit_state()

    def get_configuration_correction(self, get_inverse_mass):
        """
        The correction of the predicted configuration of an explicit
        predictor-corrector scheme by the contactino library, see
        :class:`CentralDifferencePCTS
        <sfepy.solvers.ts_solvers.CentralDifferencePCTS>` and
        :meth:`ContactinoData.correct_configuration()
        <contactino.ContactinoData.correct_configuration>`.

        The state variable of the term has to hold the predicted
        displacements. ``get_inverse_mass(dofs)`` returns the inverse mass of
        the host on the variable DOFs `dofs` as COO triplets with indices
        local to `dofs`; DOFs without an entry are immovable.

        Returns the variable DOFs of the contact surface and the displacement
        correction in them.
        """
        virtual, state = self.get_args(['virtual', 'state'])
        ci = self._get_evaluated_ci()
        uc = ci.get_surface_displacements(state)
        rows, cols, vals = get_inverse_mass(ci.global_dofs)
        duc = ci.cd.correct_configuration(uc, rows, cols, vals)
        return ci.global_dofs, duc

    def call_function(self, out, fargs):
        try:
            out, status = self.function(out, *fargs)

        except (RuntimeError, ValueError):
            terms.errclear()
            raise

        if status:
            terms.errclear()
            raise ValueError('term evaluation failed! (%s)' % self.name)

        return out, status

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):
        if mode == 'weak':
            out, status = self.call_function(None, fargs)

        else:
            out = nm.empty(shape, dtype=nm.float64)
            status = self.call_function(out, fargs)

        return out, status

    @staticmethod
    def function_weak(out, out_cc):
        return out_cc, 0

    @staticmethod
    def integrate(out, val_qp, geo, fmode):
        if fmode == 2:
            out[:] = val_qp
            status = 0

        else:
            status = geo.integrate(out, val_qp, fmode)

        return out, status

    @staticmethod
    def function(out, fun, *args):
        return fun(out, *args)

    def get_fargs(self, algorithm, epss, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return self.eval_contact(algorithm, epss, 0.0, 0.0, virtual, state,
                                 mode, term_mode, diff_var)

    def eval_contact(self, algorithm, epss, epst, mu, virtual, state,
                     mode, term_mode, diff_var):
        """
        The common implementation of :meth:`get_fargs()` of the frictionless
        and frictional contact terms.
        """
        geo, _ = self.get_mapping(virtual)

        ci = self.get_contact_info(geo, state, epss=epss, epst=epst, mu=mu,
                                   algorithm=algorithm)

        if mode == 'weak':
            uc = ci.get_surface_displacements(state)

            # A failed contact evaluation (e.g. a degenerate contact facet in
            # a trial state of a nonlinear solver) is reported as ValueError,
            # so that the line search of the Newton solver can reduce the
            # step, as for the "warp violation" of the hyperelastic terms.
            try:
                if diff_var is None:
                    rc = ci.cd.eval_residual(uc)
                    out_cc = (rc, ci.global_dofs, state)

                else:
                    rc, rows, cols, vals = ci.cd.eval_residual_and_tangent(uc)
                    # contactino returns the full (non-symmetric)
                    # two-half-pass tangent as COO triplets with repeated
                    # (row, col) pairs that the COO -> CSR conversion sums up.
                    # Do not mirror.
                    out_cc = (vals, ci.global_dofs[rows],
                              ci.global_dofs[cols], state, state)

            except RuntimeError as exc:
                raise ValueError('%s: contact evaluation failed! (%s)'
                                 % (self.name, exc))

            return self.function_weak, out_cc

        elif mode in ('el_avg', 'qp'):
            fmode = {'el_avg' : 1, 'qp' : 2}[mode]

            if term_mode == 'gap':
                # The contactino C API does not expose the gaps (yet).
                if not self._gap_warned:
                    output('WARNING: %s: gap evaluation is not available'
                           ' with contactino, returning zeros!' % self.name)
                    self._gap_warned = True

                gap = nm.zeros((ci.sd.n_fa, geo.n_qp, 1, 1), dtype=nm.float64)
                return self.integrate, gap, geo, fmode

            else:
                raise ValueError('unsupported term mode in %s! (%s)'
                                 % (self.name, term_mode))

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, algorithm, epss, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, 1, 1), state.dtype

class ContactFrictionTerm(ContactTerm):
    r"""
    Contact term with a penalty function and the penalty-regularized Coulomb
    friction, evaluated by the contactino library.

    The normal penalty traction is :math:`\varepsilon_N \langle g_N(\ul{u})
    \rangle` as in :class:`ContactTerm`. The frictional traction
    :math:`\ul{t}_T` is the Coulomb friction regularized by the tangential
    penalty :math:`\varepsilon_T` (elastic slip stiffness): :math:`|\ul{t}_T|
    \leq \mu \varepsilon_N \langle g_N(\ul{u}) \rangle`, with the slip
    increment measured from the state committed at the end of the previous
    time step (backward Euler, see :meth:`ContactTerm.advance()`); it acts
    from the second time step in contact. The contact parameters cannot
    change during the simulation.

    This term has a dynamic connectivity of DOFs in its region.

    :Definition:

    .. math::
        \int_{\Gamma_{c}} \left(\varepsilon_N \langle g_N(\ul{u}) \rangle
        \ul{n} + \ul{t}_T(\ul{u})\right) \cdot \ul{v}

    :Arguments:
        - opt_material : the contact algorithm, see :class:`ContactTerm`
        - material_n   : :math:`\varepsilon_N`
        - material_t   : :math:`\varepsilon_T`
        - material_mu  : :math:`\mu`
        - virtual      : :math:`\ul{v}`
        - state        : :math:`\ul{u}`
    """
    name = 'dw_contact_friction'
    arg_types = ('opt_material', 'material_n', 'material_t', 'material_mu',
                 'virtual', 'state')
    arg_shapes = [{'opt_material' : '.: str', 'material_n' : '.: 1',
                   'material_t' : '.: 1', 'material_mu' : '.: 1',
                   'virtual' : ('D', 'state'), 'state' : 'D'},
                  {'opt_material' : None}]

    def get_fargs(self, algorithm, epss, epst, mu, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return self.eval_contact(algorithm, epss, epst, mu, virtual, state,
                                 mode, term_mode, diff_var)

    def get_eval_shape(self, algorithm, epss, epst, mu, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        return ContactTerm.get_eval_shape(self, algorithm, epss, virtual,
                                          state, mode=mode,
                                          term_mode=term_mode,
                                          diff_var=diff_var, **kwargs)

class ContactMassTerm(ContactTerm):
    r"""
    Inertia (mass) penalty of the bipenalty method for the contact constraint
    enforced by :class:`ContactTerm` (``dw_contact``), evaluated by the
    contactino library.

    The bipenalty method [1]_, [2]_ adds to the stiffness penalty
    :math:`\varepsilon_N \langle g_N \rangle` the inertia penalty
    :math:`\varepsilon_M \ddot{g}_N` on the same constraint, so that the
    stiffness penalty can be arbitrarily large without reducing the critical
    time step of explicit time integration, provided :math:`\varepsilon_N /
    \varepsilon_M \leq \omega_{max}^2` of the unpenalized discretization. The
    term assembles the penalty mass matrix :math:`\varepsilon_M G(\ul{u})`,
    where :math:`G = \sum_q w_q N_q N_q^T` is the gap operator of the active
    contact constraints in the configuration given by the displacements
    :math:`\ul{u}`, and its product with the state variable
    :math:`\ddot{\ul{u}}`.

    The term belongs to the mass block of the elastodynamics solvers: use it
    as ``dw_contact_mass.i.Contact(contact.epsm, v, ddu)`` with
    ``auto_transform_equations``. The state variable is then the second time
    derivative ``dd<u>`` of the displacements and the displacement variable
    ``<u>`` is looked up by its name. If the state variable name does not
    start with ``dd``, the state itself is taken as the displacements.

    .. [1] H. Askes, M. Caramés-Saddler, A. Rodríguez-Ferran, Bipenalty method
           for time domain computational dynamics, Proc. R. Soc. A 466 (2010).
    .. [2] J. Kopačka, A. Tkachuk, D. Gabriel, R. Kolman, M. Bischoff,
           J. Plešek, On stability and reflection-transmission analysis of the
           bipenalty method in contact-impact problems: A one-dimensional,
           homogeneous case study, Int. J. Numer. Meth. Engng 113 (2018).

    :Definition:

    .. math::
        \int_{\Gamma_{c}} \varepsilon_M \ddot{g}_N(\ul{u}) \ul{n} \ul{v}
        \approx \ul{v}^T \varepsilon_M G(\ul{u}) \ddot{\ul{u}}

    :Arguments:
        - opt_material : the contact algorithm, see :class:`ContactTerm`
        - material : :math:`\varepsilon_M`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ddot{\ul{u}}`
    """
    name = 'dw_contact_mass'
    arg_types = ('opt_material', 'material', 'virtual', 'state')
    arg_shapes = [{'opt_material' : '.: str', 'material' : '.: 1',
                   'virtual' : ('D', 'state'), 'state' : 'D'},
                  {'opt_material' : None}]
    integration = 'facet'

    def get_displacement_variable(self, state):
        """
        The displacement variable ``<u>`` for the state ``dd<u>``.
        """
        name = state.name
        if name.startswith('dd') and (state._variables is not None):
            try:
                return state._variables[name[2:]]

            except (KeyError, IndexError):
                pass

        return state

    def get_fargs(self, algorithm, epsm, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        import scipy.sparse as sps

        if mode != 'weak':
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

        geo, _ = self.get_mapping(virtual)

        # The contact state is given by the displacements. The gap operator
        # G = sum_q w_q j_q N_q N_q^T of the active constraints carries neither
        # the penalty parameter nor the geometric terms.
        uvar = self.get_displacement_variable(state)
        ci = self.get_contact_info(geo, uvar, algorithm=algorithm)
        uc = ci.get_surface_displacements(uvar)
        rows, cols, vals = ci.cd.eval_gap_operator(uc)

        if diff_var is None:
            ac = ci.get_surface_displacements(state)
            mtx_g = sps.coo_matrix((vals, (rows, cols)),
                                   shape=(ci.neq, ci.neq)).tocsr()
            out_cc = (epsm * (mtx_g @ ac), ci.global_dofs, state)

        else:
            out_cc = (epsm * vals, ci.global_dofs[rows], ci.global_dofs[cols],
                      state, state)

        return self.function_weak, out_cc

    def get_eval_shape(self, epsm, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        return (n_el, 1, 1, 1), state.dtype

class ContactIPCTerm(Term):
    r"""
    Contact term based on IPC Toolkit.

    This term has a dynamic connectivity of DOFs in its region.

    :Definition:

    .. math::
        \int_{\Gamma_{c}} \sum_{k \in C} \nabla b(d_k, \ul{u})

    :Arguments:
        - material_m : :math:`m`
        - material_k : :math:`k`
        - material_d : :math:`\hat{d}`
        - material_p : PSDProjectionMethod
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_contact_ipc'
    arg_types = ('material_m', 'material_k', 'material_d', 'material_p',
                 'virtual', 'state')
    arg_shapes = {'material_m' : '.: 1', 'material_k' : '.: 1',
                  'material_d' : '.: 1', 'material_p' : '.: str',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    integration = 'facet'
    geometries = ['3_4', '3_8']

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        import ipctk

        self.ipc = ipctk
        self.ci = None

    def call_function(self, out, fargs):
        try:
            out, status = self.function(out, *fargs)

        except (RuntimeError, ValueError):
            terms.errclear()
            raise

        if status:
            terms.errclear()
            raise ValueError('term evaluation failed! (%s)' % self.name)

        return out, status

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):
        if mode == 'weak':
            out, status = self.call_function(None, fargs)

        else:
            out = nm.empty(shape, dtype=nm.float64)
            status = self.call_function(out, fargs)

        return out, status

    @staticmethod
    def function_weak(out, out_cc):
        return out_cc, 0

    @staticmethod
    def function(out, fun, *args):
        return fun(out, *args)

    def get_contact_info(self, state):
        from sfepy.discrete.fem.mesh import Mesh
        from sfepy.mesh.mesh_tools import triangulate
        from sfepy.discrete.fem.geometry_element import create_geometry_elements
        from sfepy.discrete.variables import create_adof_conn

        if self.ci is not None:
            return self.ci

        tdim = self.region.kind_tdim
        if tdim not in (1, 2):
            raise ValueError(
                'contact region topological dimension must be 1 or 2!'
                f' (is {tdim})'
            )

        smesh = Mesh.from_region(self.region, self.region.domain.mesh,
                                 localize=True, is_surface=True)

        if '2_4' in smesh.descs:
            smesh = triangulate(smesh)

        gels = create_geometry_elements()

        cmesh = smesh.cmesh
        cmesh.set_local_entities(gels)
        cmesh.setup_entities()

        edges = cmesh.get_conn(1, 0).indices.astype(nm.int32)
        edges = edges.reshape((-1, 2))
        if tdim == 2:
            faces = cmesh.get_conn(2, 0).indices.astype(nm.int32)
            faces = faces.reshape((-1, 3))
            collision_mesh = self.ipc.CollisionMesh(smesh.coors, edges, faces)

        else:
            faces = None
            collision_mesh = self.ipc.CollisionMesh(smesh.coors, edges)

        nods = state.field.get_dofs_in_region(self.region, merge=True)

        econn = state.field.get_econn('facet', self.region)
        dofs = nm.unique(
            create_adof_conn(nm.arange(state.n_dof), econn, smesh.dim, 0)
        )

        self.ci = Struct(smesh=smesh, edges=edges, faces=faces, nods=nods,
                         econn=econn, dofs=dofs, dim=smesh.dim,
                         collision_mesh=collision_mesh,
                         prev_min_distance=None,
                         min_distance=None,
                         it=0,
                         ls_it=0,
                         barrier_potential=None,
                         barrier_stiffness=None,
                         max_barrier_stiffness=None,
                         e_grad_full=None,
                         bp_grad_norm=None,
                         e_grad_norm=None)
        return self.ci

    def get_fargs(self, avg_mass, stiffness, dhat, spd_projection,
                  virtual, state,
                  mode=None, term_mode=None, diff_var=None, ci=None,
                  **kwargs):
        if diff_var is not None:
            ci = self.ci

        if ci is None:
            output('WARNING: ci argument was not passed to dw_contact_ipc')
            ci = self.get_contact_info(state)

        collision_mesh = ci.collision_mesh

        uvec = state().reshape((-1, ci.dim))[ci.nods]

        vertices = collision_mesh.rest_positions + uvec
        try:
            collisions = self.ipc.NormalCollisions()

        except AttributeError:
            collisions = self.ipc.Collisions()

        collisions.build(collision_mesh, vertices, dhat)

        ci.min_distance = collisions.compute_minimum_distance(
            collision_mesh, vertices,
        )
        ci.bbox_diagonal = self.ipc.world_bbox_diagonal_length(vertices)

        B = self.ipc.BarrierPotential(dhat)

        actual_stiffness = stiffness
        if actual_stiffness == 0.0: # Adaptive barrier stiffness
            if ci.it == 0:
                # This should be done at the start of each time step.
                bp_grad = B.gradient(collisions, collision_mesh, vertices)

                if ci.e_grad_full is not None:
                    e_grad = ci.e_grad_full[ci.dofs]

                else:
                    e_grad = nm.zeros_like(bp_grad)

                (ci.barrier_stiffness,
                 ci.max_barrier_stiffness) = self.ipc.initial_barrier_stiffness(
                     ci.bbox_diagonal,
                     B.barrier,
                     dhat, avg_mass,
                     e_grad, bp_grad,
                 )
                actual_stiffness = ci.barrier_stiffness
                ci.prev_min_distance = ci.min_distance
                ci.e_grad_norm = nm.linalg.norm(e_grad)
                ci.bp_grad_norm = nm.linalg.norm(bp_grad)

            elif (ci.ls_it == 0) and (diff_var is None):
                # Do not update during line-search!
                actual_stiffness = self.ipc.update_barrier_stiffness(
                    ci.prev_min_distance, ci.min_distance,
                    ci.max_barrier_stiffness, ci.barrier_stiffness,
                    ci.bbox_diagonal,
                )
                ci.barrier_stiffness = actual_stiffness
                ci.prev_min_distance = ci.min_distance

            else:
                actual_stiffness = ci.barrier_stiffness

        else:
            ci.barrier_stiffness = actual_stiffness
            ci.e_grad_norm = 0.0
            ci.bp_grad_norm = 0.0

        ci.barrier_potential = B(collisions, collision_mesh, vertices)

        if self.verbosity > 0:
            output('min. distance::', ci.min_distance)
            output('barrier potential:', ci.barrier_potential)
            output('barrier stiffness:', actual_stiffness)

        if diff_var is None:
            bp_grad = B.gradient(collisions, collision_mesh, vertices)

            out_cc = (actual_stiffness * bp_grad,
                      ci.dofs, state)

        else:
            proj = getattr(self.ipc.PSDProjectionMethod, spd_projection,
                           self.ipc.PSDProjectionMethod.NONE)
            bp_hess = B.hessian(
                collisions, collision_mesh, vertices,
                project_hessian_to_psd=proj,
            )

            ir, ic = bp_hess.nonzero()
            if len(ir):
                vals = bp_hess[ir, ic].A1

            else:
                vals = nm.array([], dtype=bp_hess.dtype)

            out_cc = (actual_stiffness * vals, ci.dofs[ir],
                      ci.dofs[ic], state, state)

        return self.function_weak, out_cc
