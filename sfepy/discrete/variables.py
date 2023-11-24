"""
Classes of variables for equations/terms.
"""
from __future__ import print_function
from __future__ import absolute_import

import numpy as nm

from sfepy.base.base import (real_types, complex_types, assert_, get_default,
                             output, OneTypeList, Container, Struct,
                             iter_dict_of_lists)
from sfepy.base.timing import Timer
import sfepy.linalg as la
from sfepy.discrete.functions import Function
from sfepy.discrete.conditions import get_condition_value
from sfepy.discrete.integrals import Integral
from sfepy.discrete.common.dof_info import (DofInfo, EquationMap,
                                            expand_nodes_to_equations,
                                            is_active_bc)
from sfepy.discrete.fem.lcbc_operators import LCBCOperators
from sfepy.discrete.common.mappings import get_physical_qps
from sfepy.discrete.evaluate_variable import eval_real, eval_complex
import six
from six.moves import range

is_state = 0
is_virtual = 1
is_parameter = 2
is_field = 10

def create_adof_conns(conn_info, var_indx=None, active_only=True, verbose=True):
    """
    Create active DOF connectivities for all variables referenced in
    `conn_info`.

    If a variable has not the equation mapping, a trivial mapping is assumed
    and connectivity with all DOFs active is created.

    DOF connectivity key is a tuple ``(primary variable name, region name,
    type, trace_region)``.

    Notes
    -----
    If `active_only` is False, the DOF connectivities contain all DOFs, with
    the E(P)BC-constrained ones stored as `-1 - <DOF number>`, so that the full
    connectivities can be reconstructed for the matrix graph creation.
    """
    var_indx = get_default(var_indx, {})

    def _create(var, econn):
        offset = var_indx.get(var.name, slice(0, 0)).start
        if var.eq_map is None:
            eq = nm.arange(var.n_dof, dtype=nm.int32)

        else:
            if isinstance(var, DGFieldVariable):
                eq = nm.arange(var.n_dof, dtype=nm.int32)
            else:
                if active_only:
                    eq = var.eq_map.eq
                else:
                    eq = nm.arange(var.n_dof, dtype=nm.int32)
                    eq[var.eq_map.eq_ebc] = -1 - (var.eq_map.eq_ebc + offset)
                    eq[var.eq_map.master] = eq[var.eq_map.slave]

        adc = create_adof_conn(eq, econn, var.n_components, offset)

        return adc

    def _assign(adof_conns, dof_conn_type, region, var, field, trace_region):
        key = (var.name, region.name, dof_conn_type, trace_region)
        if not key in adof_conns:
            econn = field.get_econn(dof_conn_type, region, trace_region)
            if econn is None: return

            adof_conns[key] = _create(var, econn)

        if trace_region is not None:
            key = (var.name, region.name, dof_conn_type, None)
            if not key in adof_conns:
                econn = field.get_econn(dof_conn_type, region,
                                        trace_region=None)

                adof_conns[key] = _create(var, econn)

    if verbose:
        output('setting up dof connectivities...')
        timer = Timer(start=True)

    adof_conns = {}

    for key, ii, info in iter_dict_of_lists(conn_info, return_keys=True):
        if info.primary is not None:
            var = info.primary
            field = var.get_field()

            region = info.get_region()
            field.setup_extra_data(info)

            mreg_name = info.get_region_name(can_trace=False)
            mreg_name = None if region.name == mreg_name else mreg_name
            dct = info.dof_conn_types[var.name]
            _assign(adof_conns, dct, region, var, field, mreg_name)

        if info.has_virtual and info.trace_region is None:
            var = info.virtual
            field = var.get_field()
            field.setup_extra_data(info)

            aux = var.get_primary()
            var = aux if aux is not None else var

            region = info.get_region(can_trace=False)
            dct = info.dof_conn_types[var.name]
            _assign(adof_conns, dct, region, var, field, None)

    if verbose:
        output('...done in %.2f s' % timer.stop())

    return adof_conns

def create_adof_conn(eq, conn, dpn, offset):
    """
    Given a node connectivity, number of DOFs per node and equation mapping,
    create the active dof connectivity.

    Locally (in a connectivity row), the DOFs are stored DOF-by-DOF (u_0 in all
    local nodes, u_1 in all local nodes, ...).

    Globally (in a state vector), the DOFs are stored node-by-node (u_0, u_1,
    ..., u_X in node 0, u_0, u_1, ..., u_X in node 1, ...).
    """
    if dpn == 1:
        aux = nm.take(eq, conn)
        adc = aux + nm.asarray(offset * (aux >= 0), dtype=nm.int32)

    else:
        n_el, n_ep = conn.shape
        adc = nm.empty((n_el, n_ep * dpn), dtype=conn.dtype)
        ii = 0
        for idof in range(dpn):
            aux = nm.take(eq, dpn * conn + idof)
            adc[:, ii : ii + n_ep] = aux + nm.asarray(offset * (aux >= 0),
                                                      dtype=nm.int32)
            ii += n_ep

    return adc

def expand_basis(basis, dpn):
    """
    Expand basis for variables with several components (DOFs per node), in a
    way compatible with :func:`create_adof_conn()`, according to `dpn`
    (DOF-per-node count).
    """
    n_c, n_bf = basis.shape[-2:]
    ebasis = nm.zeros(basis.shape[:2] + (dpn, n_bf * dpn), dtype=nm.float64)
    for ic in range(n_c):
        for ir in range(dpn):
            ebasis[..., n_c*ir+ic, ir*n_bf:(ir+1)*n_bf] = basis[..., ic, :]

    return ebasis

class Variables(Container):
    """
    Container holding instances of Variable.
    """

    @staticmethod
    def from_conf(conf, fields):
        obj = Variables()
        for key, val in six.iteritems(conf):
            var = Variable.from_conf(key, val, fields)

            obj[var.name] = var

        obj.setup_dtype()

        return obj

    def __init__(self, variables=None):
        Container.__init__(self, OneTypeList(Variable),
                           vec=None,
                           r_vec=None,
                           state=set(),
                           virtual=set(),
                           parameter=set(),
                           has_virtual_dcs=False,
                           has_lcbc=False,
                           has_lcbc_rhs=False,
                           has_eq_map=False,
                           ordered_state=[],
                           ordered_virtual=[])

        if variables is not None:
            for var in variables:
                self[var.name] = var

        self.setup_dtype()

        self.adof_conns = {}

    def __setitem__(self, ii, var):
        Container.__setitem__(self, ii, var)

        if var.is_state():
            self.state.add(var.name)

        elif var.is_virtual():
            self.virtual.add(var.name)

        elif var.is_parameter():
            self.parameter.add(var.name)

        var._variables = self

        self.setup_ordering()
        self.setup_dof_info()

    def setup_dtype(self):
        """
        Setup data types of state variables - all have to be of the same
        data type, one of nm.float64 or nm.complex128.
        """
        dtypes = {nm.complex128 : 0, nm.float64 : 0}
        for var in self.iter_state(ordered=False):
            dtypes[var.dtype] += 1

        if dtypes[nm.float64] and dtypes[nm.complex128]:
            raise ValueError("All variables must have the same dtype!")

        elif dtypes[nm.float64]:
            self.dtype = nm.float64

        elif dtypes[nm.complex128]:
            self.dtype = nm.complex128

        else:
            self.dtype = None

    def link_duals(self):
        """
        Link state variables with corresponding virtual variables,
        and assign link to self to each variable instance.

        Usually, when solving a PDE in the weak form, each state
        variable has a corresponding virtual variable.
        """
        for ii in self.state:
            self[ii].dual_var_name = None

        for ii in self.virtual:
            vvar = self[ii]
            try:
                self[vvar.primary_var_name].dual_var_name = vvar.name
            except IndexError:
                pass

    def get_dual_names(self):
        """
        Get names of pairs of dual variables.

        Returns
        -------
        duals : dict
           The dual names as virtual name : state name pairs.
        """
        duals = {}
        for name in self.virtual:
            duals[name] = self[name].primary_var_name

        return duals

    def setup_ordering(self):
        """
        Setup ordering of variables.
        """
        self.link_duals()

        is_given = [self[name]._order is not None for name in self.state]
        if any(is_given) and not all(is_given):
            raise ValueError('either all or none of state variables have to'
                             ' be created with a given order!')

        if all(is_given):
            aux = [self[name]._order for name in self.state]
            orders = [(self[name]._order, name) for name in self.state]
        else:
            orders = [ii for ii in enumerate(self.state)]

        if len(orders):
            self.ordered_state = [name for order, name in sorted(orders)]
            self.ordered_virtual = [var.dual_var_name
                                    for var in self.iter_state(ordered=True)
                                    if var.dual_var_name is not None]

    def has_virtuals(self):
        return len(self.virtual) > 0

    def _create_dof_info(self, ordered_vars, flag, active=False):
        orders = {}
        di = DofInfo(f'{flag}_dof_info')
        for var_name in ordered_vars:
            var = self[var_name]
            order = var._order
            if order in orders:
                di.append_variable(var, active=active, shared=orders[order])
            else:
                di.append_variable(var, active=active)
                if order is not None:
                    orders[order] = var_name

        return di

    def setup_dof_info(self, make_virtual=False):
        """
        Setup global DOF information.
        """
        self.di = self._create_dof_info(self.ordered_state, 'state')

        if make_virtual:
            self.vdi = self._create_dof_info(self.ordered_virtual, 'virtual')
        else:
            self.vdi = self.di

    def setup_lcbc_operators(self, lcbcs, ts=None, functions=None):
        """
        Prepare linear combination BC operator matrix and right-hand side
        vector.
        """
        from sfepy.discrete.common.region import are_disjoint
        if lcbcs is None:
            self.lcdi = self.adi
            return

        self.lcbcs = lcbcs

        if (ts is None) or ((ts is not None) and (ts.step == 0)):
            regs = []
            var_names = []
            for bcs in self.lcbcs:
                for bc in bcs.iter_single():
                    vns = bc.get_var_names()

                    regs.append(bc.regions[0])
                    var_names.append(vns[0])
                    if bc.regions[1] is not None:
                        regs.append(bc.regions[1])
                        var_names.append(vns[1])

            for i0 in range(len(regs) - 1):
                for i1 in range(i0 + 1, len(regs)):
                    if ((var_names[i0] == var_names[i1])
                        and not are_disjoint(regs[i0], regs[i1])):
                        raise ValueError('regions %s and %s are not disjoint!'
                                         % (regs[i0].name, regs[i1].name))

        ops = LCBCOperators('lcbcs', self, functions=functions)

        for bcs in self.lcbcs:
            for bc in bcs.iter_single():
                vns = bc.get_var_names()
                dofs = [self[vn].dofs for vn in vns if vn is not None]
                bc.canonize_dof_names(*dofs)

                if not is_active_bc(bc, ts=ts, functions=functions):
                    continue

                output('lcbc:', bc.name)

                ops.add_from_bc(bc, ts)

        aux = ops.make_global_operator(self.adi)
        self.mtx_lcbc, self.vec_lcbc, self.lcdi = aux

        self.has_lcbc = self.mtx_lcbc is not None
        self.has_lcbc_rhs = self.vec_lcbc is not None

    def get_lcbc_operator(self):
        if self.has_lcbc:
            return self.mtx_lcbc

        else:
            raise ValueError('no LCBC defined!')

    def equation_mapping(self, ebcs, epbcs, ts, functions, problem=None,
                         active_only=True):
        """
        Create the mapping of active DOFs from/to all DOFs for all state
        variables.

        Parameters
        ----------
        ebcs : Conditions instance
            The essential (Dirichlet) boundary conditions.
        epbcs : Conditions instance
            The periodic boundary conditions.
        ts : TimeStepper instance
            The time stepper.
        functions : Functions instance
            The user functions for boundary conditions.
        problem : Problem instance, optional
            The problem that can be passed to user functions as a context.
        active_only : bool
            If True, the active DOF info ``self.adi`` uses the reduced (active
            DOFs only) numbering. Otherwise it is the same as ``self.di``.

        Returns
        -------
        active_bcs : set
            The set of boundary conditions active in the current time.
        """
        self.ebcs = ebcs
        self.epbcs = epbcs

        ##
        # Assing EBC, PBC to variables and regions.
        if ebcs is not None:
            self.bc_of_vars = self.ebcs.group_by_variables()

        else:
            self.bc_of_vars = {}

        if epbcs is not None:
            self.bc_of_vars = self.epbcs.group_by_variables(self.bc_of_vars)

        ##
        # List EBC nodes/dofs for each variable.
        active_bcs = set()
        for var_name in self.di.var_names:
            var = self[var_name]
            bcs = self.bc_of_vars.get(var.name, None)

            var_di = self.di.get_info(var_name)
            active = var.equation_mapping(bcs, var_di, ts, functions,
                                          problem=problem)
            active_bcs.update(active)

            if self.has_virtual_dcs:
                vvar = self[var.dual_var_name]
                vvar_di = self.vdi.get_info(var_name)
                active = vvar.equation_mapping(bcs, vvar_di, ts, functions,
                                               problem=problem)
                active_bcs.update(active)

        self.adi = self._create_dof_info(self.ordered_state, 'active_state',
                                         active=active_only)

        if self.has_virtual_dcs:
            self.avdi = self._create_dof_info(self.ordered_virtual,
                                              'active_virtual',
                                              active=active_only)
        else:
            self.avdi = self.adi

        self.has_eq_map = True

        return active_bcs

    def get_indx(self, var_name, reduced=False, allow_dual=False):
        var = self[var_name]

        if not var.is_state():
            if allow_dual and var.is_virtual():
                var_name = var.primary_var_name
            else:
                msg = '%s is not a state part' % var_name
                raise IndexError(msg)

        if reduced:
            return self.adi.indx[var_name]
        else:
            return self.di.indx[var_name]

    def get_matrix_shape(self):
        if not self.has_eq_map:
            raise ValueError('call equation_mapping() first!')

        return (self.avdi.n_dof_total, self.adi.n_dof_total)

    def setup_initial_conditions(self, ics, functions):
        self.ics = ics
        self.ic_of_vars = self.ics.group_by_variables()

        for var_name in self.di.var_names:
            var = self[var_name]

            ics = self.ic_of_vars.get(var.name, None)
            if ics is None: continue

            var.setup_initial_conditions(ics, self.di, functions)

        for var_name in self.parameter:
            var = self[var_name]
            if hasattr(var, 'special') and ('ic' in var.special):
                setter, sargs, skwargs = var._get_setter('ic', functions)

                var.set_data(setter(*sargs, **skwargs))
                output('IC data of %s set by %s()' % (var.name, setter.name))

    def set_adof_conns(self, adof_conns):
        """
        Set all active DOF connectivities to `self` as well as relevant
        sub-dicts to the individual variables.
        """
        self.adof_conns = adof_conns

        for var in self:
            var.adof_conns = {}

        for key, val in six.iteritems(adof_conns):
            if key[0] in self.names:
                var = self[key[0]]
                var.adof_conns[key] = val

                var = var.get_dual()
                if var is not None:
                    var.adof_conns[key] = val

    def create_vec(self):
        vec = nm.zeros((self.di.n_dof_total,), dtype=self.dtype)
        return vec

    def create_reduced_vec(self):
        vec = nm.zeros((self.adi.n_dof_total,), dtype=self.dtype)
        return vec

    def check_vec_size(self, vec, reduced=False):
        """
        Check whether the shape of the DOF vector corresponds to the
        total number of DOFs of the state variables.

        Parameters
        ----------
        vec : array
            The vector of DOF values.
        reduced : bool
            If True, the size of the DOF vector should be reduced,
            i.e. without DOFs fixed by boundary conditions.
        """
        if not reduced:
            n_dof = self.di.get_n_dof_total()

            if vec.size != n_dof:
                msg = 'incompatible data size!' \
                      ' (%d (variables) == %d (DOF vector))' \
                      % (n_dof, vec.size)
                raise ValueError(msg)

        else:
            if self.has_lcbc:
                n_dof = self.lcdi.get_n_dof_total()

            else:
                n_dof = self.adi.get_n_dof_total()

            if vec.size != n_dof:
                msg = 'incompatible data size!' \
                      ' (%d (active variables) == %d (reduced DOF vector))' \
                      % (n_dof, vec.size)
                raise ValueError(msg)

    def reduce_vec(self, vec, follow_epbc=False, svec=None):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.

        Notes
        -----
        If 'follow_epbc' is True, values of EPBC master dofs are not simply
        thrown away, but added to the corresponding slave dofs, just like when
        assembling. For vectors with state (unknown) variables it should be set
        to False, for assembled vectors it should be set to True.
        """
        if svec is None:
            svec = nm.empty((self.adi.n_dof_total,), dtype=self.dtype)
        for var in self.iter_state():
            aindx = self.adi.indx[var.name]
            svec[aindx] = var.get_reduced(vec, self.di.indx[var.name].start,
                                          follow_epbc)
        return svec

    def make_full_vec(self, svec, force_value=None, vec=None):
        """
        Make a full DOF vector satisfying E(P)BCs from a reduced DOF
        vector.

        Parameters
        ----------
        svec : array
            The reduced DOF vector.
        force_value : float, optional
            Passing a `force_value` overrides the EBC values.
        vec : array, optional
            If given, the buffer for storing the result (zeroed).

        Returns
        -------
        vec : array
            The full DOF vector.
        """
        self.check_vec_size(svec, reduced=True)

        if self.has_lcbc:
            if self.has_lcbc_rhs:
                svec = self.mtx_lcbc * svec + self.vec_lcbc

            else:
                svec = self.mtx_lcbc * svec

        if vec is None:
            vec = self.create_vec()
        for var in self.iter_state():
            indx = self.di.indx[var.name]
            aindx = self.adi.indx[var.name]
            var.get_full(svec, aindx.start, force_value, vec, indx.start)

        return vec

    def set_vec_part(self, vec, var_name, part, reduced=False):
        self.check_vec_size(vec, reduced=reduced)
        vec[self.get_indx(var_name, reduced)] = part

    def get_vec_part(self, vec, var_name, reduced=False):
        self.check_vec_size(vec, reduced=reduced)
        return vec[self.get_indx(var_name, reduced)]

    def invalidate_evaluate_caches(self, step=0):
        for var in self.iter_state():
            var.invalidate_evaluate_cache(step=step)

    def init_state(self, vec=None):
        self.init_history()

        if vec is None:
            vec = self.create_vec()

        for var in self.iter_state():
            var.locked = False

        self.set_data(vec)

        for var in self.iter_state():
            var.locked = True

        self.vec = vec

    def fill_state(self, value):
        """
        Fill the DOF vector with given value.
        """
        if self.r_vec is not None:
            self.r_vec.fill(value)

        self.vec.fill(value)
        self.invalidate_evaluate_caches(step=0)

    def apply_ebc(self, vec=None, force_values=None):
        """
        Apply essential (Dirichlet) and periodic boundary conditions
        to state all variables or the given vector `vec`.
        """
        if vec is None:
            vec = self.vec
            self.invalidate_evaluate_caches(step=0)

        for var in self.iter_state():
            var.apply_ebc(vec, self.di.indx[var.name].start, force_values)

    def apply_ic(self, vec=None, force_values=None):
        """
        Apply initial conditions to all state variables or the given
        vector `vec`.
        """
        if vec is None:
            vec = self.vec
            self.invalidate_evaluate_caches(step=0)

        for var in self.iter_state():
            var.apply_ic(vec, self.di.indx[var.name].start, force_values)

    def has_ebc(self, vec=None, force_values=None, verbose=False):
        if vec is None:
            vec = self.vec

        ok = True
        for var in self.iter_state():
            _ok = self[var.name].has_ebc(vec=vec[self.di.indx[var.name]],
                                         force_values=force_values)
            ok = ok and _ok

            if verbose:
                output(f'variable {var.name} has E(P)BC:', _ok)

        return ok

    def set_data(self, data, step=0, ignore_unknown=False,
                 preserve_caches=False):
        """
        Set data (vectors of DOF values) of variables.

        Parameters
        ----------
        data : array
            The state vector or dictionary of {variable_name : data vector}.
        step : int, optional
            The time history step, 0 (default) = current.
        ignore_unknown : bool, optional
            Ignore unknown variable names if `data` is a dict.
        preserve_caches : bool
            If True, do not invalidate evaluate caches of variables.
        """
        if data is None: return

        if isinstance(data, dict):

            for key, val in six.iteritems(data):
                try:
                    var = self[key]

                except (ValueError, IndexError):
                    if ignore_unknown:
                        pass

                    else:
                        raise KeyError('unknown variable! (%s)' % key)

                else:
                    var.set_data(val, step=step,
                                 preserve_caches=preserve_caches)

        elif isinstance(data, nm.ndarray):
            self.check_vec_size(data)

            for ii in self.state:
                var = self[ii]
                var.set_data(data, self.di.indx[var.name], step=step,
                             preserve_caches=preserve_caches)

        else:
            raise ValueError('unknown data class! (%s)' % data.__class__)

    def set_reduced_state(self, r_vec, preserve_caches=False):
        """
        Set the reduced DOF vector, with EBC and PBC DOFs removed.

        Parameters
        ----------
        r_vec : array
            The reduced DOF vector corresponding to the variables.
        preserve_caches : bool
            If True, do not invalidate evaluate caches of variables.
        """
        self.vec[:] = self.make_full_vec(r_vec)

        if self.has_lcbc:
            self.r_vec = r_vec

        if not preserve_caches:
            self.invalidate_evaluate_caches(step=0)

    def get_reduced_state(self, follow_epbc=False, force=False):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.
        """
        if self.has_lcbc:
            if self.r_vec is None:
                if force:
                    r_vec = self.reduce_vec(self.vec, follow_epbc=follow_epbc)
                    # This just sets the correct vector size (wrong values)!
                    r_vec = self.mtx_lcbc.T * r_vec

                else:
                    raise ValueError('Reduced state DOFs are not available!')

            else:
                r_vec = self.r_vec

        else:
            r_vec = self.reduce_vec(self.vec, follow_epbc=follow_epbc)

        return r_vec

    def set_full_state(self, vec, force=False, preserve_caches=False):
        """
        Set the full DOF vector (including EBC and PBC DOFs). If
        `var_name` is given, set only the DOF sub-vector corresponding
        to the given variable. If `force` is True, setting variables
        with LCBC DOFs is allowed.
        """
        if self.has_lcbc:
            if not force:
                raise ValueError('cannot set full DOF vector with LCBCs!')

            else:
                self.r_vec = None

        self.vec[:] = vec
        if not preserve_caches:
            self.invalidate_evaluate_caches(step=0)

    def __call__(self, var_name=None):
        """
        Get the full DOF vector (including EBC and PBC DOFs). If
        `var_name` is given, return only the DOF vector corresponding to
        the given variable.
        """
        if var_name is None:
            out = self.vec

        else:
            out = self.vec[self.di.indx[var_name]]

        return out

    def set_state(self, vec, reduced=False, force=False, preserve_caches=False,
                  apply_ebc=False):
        preserve_caches = preserve_caches and (not apply_ebc)
        if reduced:
            self.set_reduced_state(vec, preserve_caches=preserve_caches)
            if apply_ebc:
                self.apply_ebc()

        else:
            if apply_ebc:
                self.apply_ebc(vec)
            self.set_full_state(vec, force=force,
                                preserve_caches=preserve_caches)

    def get_state(self, reduced=False, follow_epbc=False, force=False):
        if reduced:
            vec = self.get_reduced_state(follow_epbc=follow_epbc, force=force)

        else:
            vec = self()

        return vec

    def set_state_parts(self, parts, vec=None, force=False):
        """
        Set parts of the DOF vector corresponding to individual state
        variables.

        Parameters
        ----------
        parts : dict
            The dictionary of the DOF vector parts.
        force : bool
            If True, proceed even with LCBCs present.
        """
        if self.has_lcbc and not force:
            raise ValueError('cannot set full DOF vector with LCBCs!')

        if vec is None:
            vec = self.vec
            self.invalidate_evaluate_caches(step=0)

        else:
            self.check_vec_size(vec, reduced=False)

        for key, part in parts.items():
            vec[self.di.indx[key]] = part

    def get_state_parts(self, vec=None):
        """
        Return parts of a state vector corresponding to individual state
        variables.

        Parameters
        ----------
        vec : array, optional
            The state vector. If not given, then the data stored in the
            variables are returned instead.

        Returns
        -------
        out : dict
            The dictionary of the state parts.
        """
        if vec is None:
            vec = self.vec

        else:
            self.check_vec_size(vec, reduced=False)

        out = {}
        for var in self.iter_state():
            out[var.name] = vec[self.di.indx[var.name]]

        return out

    def create_output(self, vec=None, fill_value=None, var_info=None,
                      extend=True, linearization=None):
        """
        Creates an output dictionary with state variables data, that can be
        passed as 'out' kwarg to :func:`Mesh.write()`.

        Then the dictionary entries are formed by components of the
        state vector corresponding to unknown variables according to
        kind of linearization given by `linearization`.
        """
        if vec is None:
            vec = self.vec

        di = self.di

        if var_info is None:
            self.check_vec_size(vec)

            var_info = {}
            for name in di.var_names:
                var_info[name] = (False, name)

        out = {}
        for key, indx in six.iteritems(di.indx):
            var = self[key]

            if key not in list(var_info.keys()): continue
            is_part, name = var_info[key]

            if is_part:
                aux = vec

            else:
                aux = vec[indx]

            out.update(var.create_output(aux, key=name, extend=extend,
                                         fill_value=fill_value,
                                         linearization=linearization))

        return out

    def iter_state(self, ordered=True):

        if ordered:
            for ii in self.ordered_state:
                yield self[ii]

        else:
            for ii in self.state:
                yield self[ii]

    def init_history(self):
        for var in self.iter_state():
            var.init_history()

    def time_update(self, ts, functions, verbose=True):
        if verbose:
            output('updating variables...')

        for var in self:
            var.time_update(ts, functions)

        if verbose:
            output('...done')

    def advance(self, ts):
        for var in self.iter_state():
            var.advance(ts)

class Variable(Struct):

    @staticmethod
    def from_conf(key, conf, fields):
        aux = conf.kind.split()
        if len(aux) == 2:
            kind, family = aux

        elif len(aux) == 3:
            kind, family = aux[0], '_'.join(aux[1:])

        else:
            raise ValueError('variable kind is 2 or 3 words! (%s)' % conf.kind)

        history = conf.get('history', None)
        if history is not None:
            try:
                history = int(history)
                assert_(history >= 0)

            except (ValueError, TypeError):
                raise ValueError('history must be integer >= 0! (got "%s")'
                                 % history)

        order = conf.get('order', None)
        if order is not None:
            order = int(order)

        primary_var_name = conf.get('dual', None)
        if primary_var_name is None:
            if hasattr(conf, 'like'):
                primary_var_name = get_default(conf.like, '(set-to-None)')

            else:
                primary_var_name = None

        special = conf.get('special', None)

        if family == 'field':
            try:
                fld = fields[conf.field]
            except IndexError:
                msg = 'field "%s" does not exist!' % conf.field
                raise KeyError(msg)

            if "DG" in fld.family_name:
                obj = DGFieldVariable(conf.name, kind, fld, order, primary_var_name,
                                special=special, key=key, history=history)
            else:
                obj = FieldVariable(conf.name, kind, fld, order, primary_var_name,
                                    special=special, key=key, history=history)

        else:
            raise ValueError('unknown variable family! (%s)' % family)

        return obj

    def __init__(self, name, kind, order=None, primary_var_name=None,
                 special=None, flags=None, **kwargs):
        Struct.__init__(self, name=name, locked=False, **kwargs)

        self.flags = set()
        if flags is not None:
            for flag in flags:
                self.flags.add(flag)

        self.indx = slice(None)
        self.n_dof = None
        self.step = 0
        self.dt = 1.0
        self.initial_condition = None
        self.dual_var_name = None
        self.eq_map = None

        if self.is_virtual():
            self.data = None

        else:
            self.data = []
            self.data.append(None)

        self._set_kind(kind, order, primary_var_name, special=special)

    def _set_kind(self, kind, order, primary_var_name, special=None):
        if kind == 'unknown':
            self.flags.add(is_state)
            self._order = order

            self.dof_name = self.name

        elif kind == 'test':
            if primary_var_name == self.name:
                raise ValueError('primary variable for %s cannot be %s!'
                                 % (self.name, primary_var_name))

            self.flags.add(is_virtual)
            msg = 'test variable %s: related unknown missing' % self.name
            self.primary_var_name = get_default(primary_var_name, None, msg)
            self.dof_name = self.primary_var_name

        elif kind == 'parameter':
            self.flags.add(is_parameter)
            msg = 'parameter variable %s: related unknown missing' % self.name
            self.primary_var_name = get_default(primary_var_name, None, msg)
            if self.primary_var_name == '(set-to-None)':
                self.primary_var_name = None
                self.dof_name = self.name

            else:
                self.dof_name = self.primary_var_name

            if special is not None:
                self.special = special

        else:
            raise NotImplementedError('unknown variable kind: %s' % kind)

        self.kind = kind

    def _setup_dofs(self, n_nod, n_components, val_shape):
        """
        Setup number of DOFs and  DOF names.
        """
        self.n_nod = n_nod
        self.n_components = n_components
        self.val_shape = val_shape

        self.n_dof = self.n_nod * self.n_components

        self.dofs = [self.dof_name + ('.%d' % ii)
                     for ii in range(self.n_components)]

    def get_primary(self):
        """
        Get the corresponding primary variable.

        Returns
        -------
        var : Variable instance
            The primary variable, or `self` for state
            variables or if `primary_var_name` is None, or None if no other
            variables are defined.
        """
        if self.is_state():
            var = self

        elif self.primary_var_name is not None:
            if ((self._variables is not None)
                and (self.primary_var_name in self._variables.names)):
                var = self._variables[self.primary_var_name]

            else:
                var = None

        else:
            var = self

        return var

    def get_dual(self):
        """
        Get the dual variable.

        Returns
        -------
        var : Variable instance
            The primary variable for non-state variables, or the dual
            variable for state variables.
        """
        if self.is_state():
            if ((self._variables is not None)
                and (self.dual_var_name in self._variables.names)):
                var = self._variables[self.dual_var_name]

            else:
                var = None

        else:
            if ((self._variables is not None)
                and (self.primary_var_name in self._variables.names)):
                var = self._variables[self.primary_var_name]

            else:
                var = None

        return var

    def is_state(self):
        return is_state in self.flags

    def is_virtual(self):
        return is_virtual in self.flags

    def is_parameter(self):
        return is_parameter in self.flags

    def is_state_or_parameter(self):
        return (is_state in self.flags) or (is_parameter in self.flags)

    def is_kind(self, kind):
        return eval('self.is_%s()' % kind)

    def is_real(self):
        return self.dtype in real_types

    def is_complex(self):
        return self.dtype in complex_types

    def is_finite(self, step=0, derivative=None, dt=None):
        return nm.isfinite(self(step=step, derivative=derivative, dt=dt)).all()

    def get_primary_name(self):
        if self.is_state():
            name = self.name

        else:
            name = self.primary_var_name

        return name

    def init_history(self):
        """Initialize data of variables with history."""
        if self.history is None: return

        self.data = (self.history + 1) * [None]
        self.step = 0

    def time_update(self, ts, functions):
        """Implemented in subclasses."""
        pass

    def advance(self, ts):
        """
        Advance in time the DOF state history. A copy of the DOF vector
        is made to prevent history modification.
        """
        if self.history is None: return

        self.step = ts.step + 1

        if self.history > 0:
            # Copy the current step data to the history data, shift history,
            # initialize if needed. The current step data are left intact.
            # Note: self.data can point into Variables.vec.
            for ii in range(1, self.history + 1):
                if self.data[ii] is None:
                    self.data[ii] = self.data[ii-1].copy()

            for ii in range(self.history, 0, -1):
                self.data[ii][self.indx] = self.data[ii - 1][self.indx]

            # Advance evaluate cache.
            for step_cache in six.itervalues(self.evaluate_cache):
                steps = sorted(step_cache.keys())
                for step in steps:
                    if step is None:
                        # Special caches with possible custom advance()
                        # function.
                        for key, val in six.iteritems(step_cache[step]):
                            if hasattr(val, '__advance__'):
                                val.__advance__(ts, val)

                    elif -step < self.history:
                        step_cache[step-1] = step_cache[step]

                if len(steps) and (steps[0] is not None):
                    step_cache.pop(steps[-1])

    def init_data(self, step=0):
        """
        Initialize the dof vector data of time step `step` to zeros.
        """
        if self.is_state_or_parameter():
            self.set_constant(val=0.0, step=step)

    def set_constant(self, val=0.0, step=0):
        """
        Set the variable dof vector data of time step `step` to a scalar `val`.
        """
        data = nm.empty((self.n_dof,), dtype=self.dtype)
        data.fill(val)
        self.set_data(data, step=step)

    def set_data(self, data=None, indx=None, step=0,
                 preserve_caches=False):
        """
        Set data (vector of DOF values) of the variable.

        Parameters
        ----------
        data : array
            The vector of DOF values.
        indx : int, optional
            If given, `data[indx]` is used.
        step : int, optional
            The time history step, 0 (default) = current.
        preserve_caches : bool
            If True, do not invalidate evaluate caches of the variable.
        """
        if not self.is_state_or_parameter():
            raise ValueError(
                'Only state or parameter variables can have values!'
            )

        if self.locked:
            raise ValueError(
                f'Variable {self.name} is locked! It can be set only using'
                ' the state manipulation functions of Variables.'
            )

        data = data.ravel()

        if indx is None:
            indx = slice(0, len(data))
        else:
            indx = slice(int(indx.start), int(indx.stop))
        n_data_dof = indx.stop - indx.start

        if self.n_dof != n_data_dof:
            msg = 'incompatible data shape! (%d (variable) == %d (data))' \
                  % (self.n_dof, n_data_dof)
            raise ValueError(msg)

        elif (step > 0) or (-step >= len(self.data)):
            raise ValueError('step %d out of range! ([%d, 0])'
                             % (step, -(len(self.data) - 1)))

        else:
            self.data[step] = data
            self.indx = indx

        if not preserve_caches:
            self.invalidate_evaluate_cache(step=step)

    def __call__(self, step=0, derivative=None, dt=None):
        """
        Return vector of degrees of freedom of the variable.

        Parameters
        ----------
        step : int, default 0
            The time step (0 means current, -1 previous, ...).
        derivative : None or 'dt'
            If not None, return time derivative of the DOF vector,
            approximated by the backward finite difference.

        Returns
        -------
        vec : array
            The DOF vector. If `derivative` is None: a view of the data vector,
             otherwise: required derivative of the DOF vector
             at time step given by `step`.

        Notes
        -----
        If the previous time step is requested in step 0, the step 0
        DOF vector is returned instead.
        """
        if derivative is None:
            if (self.step == 0) and (step == -1):
                data = self.data[0]

            else:
                data = self.data[-step]

            if data is None:
                raise ValueError('data of variable are not set! (%s, step %d)' \
                                 % (self.name, step))

            return data[self.indx]

        else:
            if self.history is None:
                msg = 'set history type of variable %s to use derivatives!'\
                      % self.name
                raise ValueError(msg)
            dt = get_default(dt, self.dt)

            return (self(step=step) - self(step=step-1)) / dt

    def get_initial_condition(self):
        if self.initial_condition is None:
            return 0.0
        else:
            return self.initial_condition

class FieldVariable(Variable):
    """
    A finite element field variable.

    field .. field description of variable (borrowed)
    """

    def __init__(self, name, kind, field, order=None, primary_var_name=None,
                 special=None, flags=None, history=None, **kwargs):
        Variable.__init__(self, name, kind, order, primary_var_name,
                          special, flags, history=history, **kwargs)

        self._set_field(field)

        self.has_field = True
        self.has_bc = True
        self._variables = None

        self.clear_evaluate_cache()

    def _set_field(self, field):
        """
        Set field of the variable.

        Takes reference to a Field instance. Sets dtype according to
        field.dtype. Sets `dim` attribute to spatial dimension.
        """
        self.is_surface = field.is_surface

        self.field = field
        self._setup_dofs(field.n_nod, field.n_components, field.val_shape)

        self.flags.add(is_field)
        self.dtype = field.dtype

        self.dim = field.domain.shape.dim


    def _get_setter(self, kind, functions, **kwargs):
        """
        Get the setter function of the variable and its arguments depending in
        the setter kind.
        """
        if not (hasattr(self, 'special') and (kind in self.special)):
            return

        setter_name = self.special[kind]
        setter = functions[setter_name]

        region = self.field.region
        nod_list = self.field.get_dofs_in_region(region)
        nods = nm.unique(nod_list)

        coors = self.field.get_coor(nods)

        if kind == 'setter':
            sargs = (kwargs.get('ts'), coors)

        elif kind == 'ic':
            sargs = (coors, )

        skwargs = {'region' : region, 'variable' : self}

        return setter, sargs, skwargs

    def get_field(self):
        return self.field

    def get_mapping(self, region, integral, integration,
                    get_saved=False, return_key=False):
        """
        Get the reference element mapping of the underlying field.

        See Also
        --------
        sfepy.discrete.common.fields.Field.get_mapping
        """
        if region is None:
            region = self.field.region

        out = self.field.get_mapping(region, integral, integration,
                                     get_saved=get_saved,
                                     return_key=return_key)
        return out

    def get_dof_conn(self, region_name, dct, trace_region=None):
        """
        Get active dof connectivity of a variable.

        Notes
        -----
        The primary and dual variables must have the same Region.
        """
        if self.is_virtual():
            var = self.get_primary()
            # No primary variable can occur in single term evaluations.
            var_name = var.name if var is not None else self.name

        else:
            var_name = self.name

        if trace_region is None:
            mregion_name = None
        else:
            mregion = self.field.domain.regions[region_name]
            region = mregion.get_mirror_region(trace_region)
            region_name = region.name
            mregion_name = mregion.name

        key = (var_name, region_name, dct, mregion_name)
        dc = self.adof_conns[key]

        return dc

    def get_dof_info(self, active=False):
        details = Struct(name='field_var_dof_details',
                         n_nod=self.n_nod,
                         dpn=self.n_components)
        if active:
            n_dof = self.n_adof

        else:
            n_dof = self.n_dof

        return n_dof, details

    def time_update(self, ts, functions):
        """
        Store time step, set variable data for variables with the setter
        function.
        """
        if ts is not None:
            self.dt = ts.dt

        if hasattr(self, 'special') and ('setter' in self.special):
            setter, sargs, skwargs = self._get_setter('setter', functions,
                                                      ts=ts)

            self.set_data(setter(*sargs, **skwargs))
            output('data of %s set by %s()' % (self.name, setter.name))

    def set_from_qp(self, data_qp, integral, step=0):
        """
        Set DOFs of variable using values in quadrature points
        corresponding to the given integral.
        """
        data_vertex = self.field.average_qp_to_vertices(data_qp, integral)

        # Field nodes values.
        data = self.field.interp_v_vals_to_n_vals(data_vertex)
        data = data.ravel()
        self.indx = slice(0, len(data))

        self.data[step] = data

    def set_from_mesh_vertices(self, data):
        """
        Set the variable using values at the mesh vertices.
        """
        ndata = self.field.interp_v_vals_to_n_vals(data)
        self.set_data(ndata)

    def set_from_function(self, fun, step=0):
        """
        Set the variable data (the vector of DOF values) using a function of
        space coordinates.

        Parameters
        ----------
        fun : callable
            The function of coordinates returning DOF values of shape
            `(n_coor, n_components)`.
        step : int, optional
            The time history step, 0 (default) = current.
        """
        _, vv = self.field.set_dofs(fun, self.field.region, self.n_components)
        self.set_data(vv.ravel(), step=step)

    def equation_mapping(self, bcs, var_di, ts, functions, problem=None,
                         warn=False):
        """
        Create the mapping of active DOFs from/to all DOFs.

        Sets n_adof.

        Returns
        -------
        active_bcs : set
            The set of boundary conditions active in the current time.
        """
        self.eq_map = EquationMap('eq_map', self.dofs, var_di)

        if var_di.shared_dofs_with is not None:
            svar = self._variables[var_di.shared_dofs_with]
            eq_map = self.eq_map
            seq_map = svar.eq_map
            eq_map.eq = seq_map.eq
            eq_map.eq_ebc = seq_map.eq_ebc
            eq_map.eqi = seq_map.eqi
            eq_map.master = seq_map.master
            eq_map.n_eq = seq_map.n_eq
            eq_map.slave = seq_map.slave
            eq_map.val_ebc = seq_map.val_ebc

            self.n_adof = eq_map.n_eq

            return {}

        if bcs is not None:
            bcs.canonize_dof_names(self.dofs)
            bcs.sort()

        active_bcs = self.eq_map.map_equations(bcs, self.field, ts, functions,
                                               problem=problem, warn=warn)
        self.n_adof = self.eq_map.n_eq

        return active_bcs

    def setup_initial_conditions(self, ics, di, functions, warn=False):
        """
        Setup of initial conditions.
        """
        ics.canonize_dof_names(self.dofs)
        ics.sort()

        self.initial_condition = nm.zeros((di.n_dof[self.name],),
                                          dtype=self.dtype)
        for ic in ics:
            region = ic.region
            dofs, val = ic.dofs

            if warn:
                clean_msg = ('warning: ignoring nonexistent' \
                             ' IC node (%s) in ' % self.name)
            else:
                clean_msg = None

            nod_list = self.field.get_dofs_in_region(region)
            if len(nod_list) == 0:
                continue

            fun = get_condition_value(val, functions, 'IC', ic.name)
            if isinstance(fun, Function):
                aux = fun
                fun = lambda coors: aux(coors, ic=ic)

            nods, vv = self.field.set_dofs(fun, region, len(dofs), clean_msg)
            eq = expand_nodes_to_equations(nods, dofs, self.dofs)

            self.initial_condition[eq] = nm.ravel(vv)

    def get_data_shape(self, integral, integration='cell', region_name=None):
        """
        Get element data dimensions for given approximation.

        Parameters
        ----------
        integral : Integral instance
            The integral describing used numerical quadrature.
        integration : 'cell', 'facet', 'facet_extra', 'point' or 'custom'
            The term integration mode.
        region_name : str
            The name of the region of the integral.

        Returns
        -------
        data_shape : 5 ints
            The `(n_el, n_qp, dim, n_en, n_comp)` for volume shape kind,
            `(n_fa, n_qp, dim, n_fn, n_comp)` for surface shape kind and
            `(n_nod, 0, 0, 1, n_comp)` for point shape kind.

        Notes
        -----
        - `n_el`, `n_fa` = number of elements/facets
        - `n_qp` = number of quadrature points per element/facet
        - `dim` = spatial dimension
        - `n_en`, `n_fn` = number of element/facet nodes
        - `n_comp` = number of variable components in a point/node
        - `n_nod` = number of element nodes
        """
        aux = self.field.get_data_shape(integral, integration=integration,
                                        region_name=region_name)
        data_shape = aux + (self.n_components,)

        return data_shape

    def clear_evaluate_cache(self):
        """
        Clear current evaluate cache.
        """
        self.evaluate_cache = {}

    def invalidate_evaluate_cache(self, step=0):
        """
        Invalidate variable data in evaluate cache for time step given
        by `step`  (0 is current, -1 previous, ...).

        This should be done, for example, prior to every nonlinear
        solver iteration.
        """
        for step_cache in six.itervalues(self.evaluate_cache):
            for key in list(step_cache.keys()):
                if key == step: # Given time step to clear.
                    step_cache.pop(key)

    def evaluate(self, mode='val',
                 region=None, integral=None, integration=None,
                 step=0, time_derivative=None, trace_region=None,
                 dt=None, bf=None):
        """
        Evaluate various quantities related to the variable according to
        `mode` in quadrature points defined by `integral`.

        The evaluated data are cached in the variable instance in
        `evaluate_cache` attribute.

        Parameters
        ----------
        mode : one of 'val', 'grad', 'div', 'cauchy_strain'
            The evaluation mode.
        region : Region instance, optional
            The region where the evaluation occurs. If None, the
            underlying field region is used.
        integral : Integral instance, optional
            The integral defining quadrature points in which the
            evaluation occurs. If None, the first order volume integral
            is created. Must not be None for surface integrations.
        integration : 'cell', 'facet', 'facet_extra', or 'point'
            The term integration type. If None, it is derived from
            the region kind.
        step : int, default 0
            The time step (0 means current, -1 previous, ...).
        time_derivative : None or 'dt'
            If not None, return time derivative of the data,
            approximated by the backward finite difference.
        trace_region : None or str
            If not None, evaluate of trace of the variable on a boundary
            region.
        dt : float, optional
            The time step to be used if `derivative` is `'dt'`. If None,
            the `dt` attribute of the variable is used.
        bf : Base function, optional
            The base function to be used in 'val' mode.

        Returns
        -------
        out : array
            The 4-dimensional array of shape
            `(n_el, n_qp, n_row, n_col)` with the requested data,
            where `n_row`, `n_col` depend on `mode`.
        """
        if integration == 'custom':
            msg = 'cannot use FieldVariable.evaluate() with custom integration!'
            raise ValueError(msg)

        step_cache = self.evaluate_cache.setdefault(mode, {})
        cache = step_cache.setdefault(step, {})

        field = self.field
        if region is None:
            region = field.region

        if trace_region is not None:
            mregion = region.get_mirror_region(trace_region)
            trace_region = region.name
            region = mregion

        if (region is not field.region) and not region.is_empty:
            assert_(field.region.contains(region))

        if integral is None:
            integral = Integral('aux_1', 1)

        if integration is None:
            integration = region.kind

        geo, _, key = field.get_mapping(region, integral, integration,
                                        return_key=True)
        key += (time_derivative, trace_region)

        dct = ('cell' if integration == 'facet_extra' else integration,
               region.tdim)

        if key in cache:
            out = cache[key]

        else:
            vec = self(step=step, derivative=time_derivative, dt=dt)
            conn = field.get_econn(dct, region, trace_region)

            shape = self.get_data_shape(integral, integration, region.name)

            if self.dtype == nm.float64:
                out = eval_real(vec, conn, geo, mode, shape, bf)

            else:
                out = eval_complex(vec, conn, geo, mode, shape, bf)

            cache[key] = out

        return out

    def get_state_in_region(self, region, reshape=True, step=0):
        """
        Get DOFs of the variable in the given region.

        Parameters
        ----------
        region : Region
            The selected region.
        reshape : bool
            If True, reshape the DOF vector to a 2D array with the individual
            components as columns. Otherwise a 1D DOF array of the form [all
            DOFs in region node 0, all DOFs in region node 1, ...] is returned.
        step : int, default 0
            The time step (0 means current, -1 previous, ...).

        Returns
        -------
        out : array
            The selected DOFs.
        """
        nods = self.field.get_dofs_in_region(region, merge=True)

        eq = nm.empty((len(nods) * self.n_components,), dtype=nm.int32)
        for idof in range(self.n_components):
            eq[idof::self.n_components] = self.n_components * nods \
                                          + idof + self.indx.start

        out = self.data[step][eq]
        if reshape:
            out.shape = (len(nods), self.n_components)

        return out

    def apply_ebc(self, vec, offset=0, force_values=None):
        """
        Apply essential (Dirichlet) and periodic boundary conditions to
        vector `vec`, starting at `offset`.
        """
        eq_map = self.eq_map
        ii = offset + eq_map.eq_ebc

        # EBC,
        if force_values is None:
            vec[ii] = eq_map.val_ebc

        else:
            if isinstance(force_values, dict):
                vec[ii] = force_values[self.name]

            else:
                vec[ii] = force_values

        # EPBC.
        vec[offset+eq_map.master] = vec[offset+eq_map.slave]

    def apply_ic(self, vec, offset=0, force_values=None):
        """
        Apply initial conditions conditions to vector `vec`, starting at
        `offset`.
        """
        ii = slice(offset, offset + self.n_dof)

        if force_values is None:
            vec[ii] = self.get_initial_condition()

        else:
            if isinstance(force_values, dict):
                vec[ii] = force_values[self.name]

            else:
                vec[ii] = force_values

    def has_ebc(self, vec=None, force_values=None):
        eq_map = self.eq_map
        ii = eq_map.eq_ebc
        if force_values is None:
            if not nm.allclose(vec[ii], eq_map.val_ebc):
                return False
        else:
            if isinstance(force_values, dict):
                if not nm.allclose(vec[ii], force_values[self.name]):
                    return False
            else:
                if not nm.allclose(vec[ii], force_values):
                    return False
        # EPBC.
        if not nm.allclose(vec[eq_map.master], vec[eq_map.slave]):
            return False

        return True

    def get_reduced(self, vec, offset=0, follow_epbc=False):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.

        Notes
        -----
        The full vector starts in `vec` at `offset`. If 'follow_epbc' is True,
        values of EPBC master DOFs are not simply thrown away, but added to the
        corresponding slave DOFs, just like when assembling. For vectors with
        state (unknown) variables it should be set to False, for assembled
        vectors it should be set to True.
        """
        eq_map = self.eq_map
        ii = offset + eq_map.eqi

        r_vec = vec[ii]

        if follow_epbc:
            master = offset + eq_map.master
            slave = eq_map.eq[eq_map.slave]
            ii = slave >= 0
            la.assemble1d(r_vec, slave[ii], vec[master[ii]])

        return r_vec

    def get_full(self, r_vec, r_offset=0, force_value=None,
                 vec=None, offset=0):
        """
        Get the full DOF vector satisfying E(P)BCs from a reduced DOF
        vector.

        Notes
        -----
        The reduced vector starts in `r_vec` at `r_offset`.
        Passing a `force_value` overrides the EBC values. Optionally,
        `vec` argument can be provided to store the full vector (in
        place) starting at `offset`.
        """
        if vec is None:
            vec = nm.empty(self.n_dof, dtype=r_vec.dtype)

        else:
            vec = vec[offset:offset+self.n_dof]

        eq_map = self.eq_map
        r_vec = r_vec[r_offset:r_offset+eq_map.n_eq]

        # EBC.
        vec[eq_map.eq_ebc] = get_default(force_value, eq_map.val_ebc)

        # Reduced vector values.
        vec[eq_map.eqi] = r_vec

        # EPBC.
        vec[eq_map.master] = vec[eq_map.slave]

        unused_dofs = self.field.get('unused_dofs')
        if unused_dofs is not None:
            vec[:] = self.field.restore_substituted(vec)

        return vec

    def create_output(self, vec=None, key=None, extend=True, fill_value=None,
                      linearization=None):
        """
        Convert the DOF vector to a dictionary of output data usable by
        Mesh.write().

        Parameters
        ----------
        vec : array, optional
            An alternative DOF vector to be used instead of the variable
            DOF vector.
        key : str, optional
            The key to be used in the output dictionary instead of the
            variable name.
        extend : bool
            Extend the DOF values to cover the whole domain.
        fill_value : float or complex
           The value used to fill the missing DOF values if `extend` is True.
        linearization : Struct or None
            The linearization configuration for higher order approximations.
        """
        linearization = get_default(linearization, Struct(kind='strip'))

        if vec is None:
            vec = self()

        key = get_default(key, self.name)

        aux = nm.reshape(vec,
                         (self.n_dof // self.n_components, self.n_components))

        out = self.field.create_output(aux, self.name, dof_names=self.dofs,
                                       key=key, extend=extend,
                                       fill_value=fill_value,
                                       linearization=linearization)

        return out

    def get_element_diameters(self, cells, mode, square=False):
        """Get diameters of selected elements."""
        field = self.field
        domain = field.domain

        cells = nm.array(cells)

        diameters = nm.empty((cells.shape[0],), dtype=nm.float64)

        integral = Integral('i_tmp', 1)

        vg, _ = field.get_mapping(field.region, integral, 'cell')

        diameters = domain.get_element_diameters(cells, vg.volume, mode,
                                                 square=square)

        return diameters

    def save_as_mesh(self, filename):
        """
        Save the field mesh and the variable values into a file for
        visualization. Only the vertex values are stored.
        """
        mesh = self.field.create_mesh(extra_nodes=False)
        vec = self()

        n_dof, dpn = self.n_dof, self.n_components
        aux = nm.reshape(vec, (n_dof // dpn, dpn))

        ext = self.field.extend_dofs(aux, 0.0)

        out = {}
        if self.field.approx_order != 0:
            out[self.name] = Struct(name='output_data',
                                    mode='vertex', data=ext,
                                    var_name=self.name, dofs=self.dofs)
        else:
            ext.shape = (ext.shape[0], 1, ext.shape[1], 1)
            out[self.name] = Struct(name='output_data',
                                    mode='cell', data=ext,
                                    var_name=self.name, dofs=self.dofs)

        mesh.write(filename, io='auto', out=out)

    def has_same_mesh(self, other):
        """
        Returns
        -------
        flag : int
            The flag can be either 'different' (different meshes), 'deformed'
            (slightly deformed same mesh), or 'same' (same).
        """
        f1 = self.field
        f2 = other.field

        c1 = f1.get_coor()
        c2 = f2.get_coor()

        if c1.shape != c2.shape:
            flag = 'different'

        else:
            eps = 10.0 * nm.finfo(nm.float64).eps

            if nm.allclose(c1, c2, rtol=eps, atol=0.0):
                flag = 'same'

            elif nm.allclose(c1, c2, rtol=0.1, atol=0.0):
                flag = 'deformed'

            else:
                flag = 'different'

        return flag

    def get_interp_coors(self, strategy='interpolation', interp_term=None):
        """
        Get the physical coordinates to interpolate into, based on the strategy
        used.
        """
        if strategy == 'interpolation':
            coors = self.field.get_coor()

        elif strategy == 'projection':
            region = self.field.region
            integral = Integral(term=interp_term)
            coors = get_physical_qps(region, integral)

        else:
            raise ValueError('unknown interpolation strategy! (%s)' % strategy)

        return coors

    def evaluate_at(self, coors, mode='val', strategy='general',
                    close_limit=0.1, get_cells_fun=None,
                    cache=None, ret_cells=False,
                    ret_status=False, ret_ref_coors=False, verbose=False):
        """
        Evaluate the variable in the given physical coordinates. Convenience
        wrapper around :func:`Field.evaluate_at()
        <sfepy.discrete.common.fields.Field.evaluate_at()>`, see its
        docstring for more details.
        """
        source_vals = self().reshape((self.n_nod, self.n_components))
        out = self.field.evaluate_at(coors, source_vals,
                                     mode=mode,
                                     strategy=strategy,
                                     close_limit=close_limit,
                                     get_cells_fun=get_cells_fun,
                                     cache=cache,
                                     ret_cells=ret_cells,
                                     ret_status=ret_status,
                                     ret_ref_coors=ret_ref_coors,
                                     verbose=verbose)

        return out

    def set_from_other(self, other, strategy='projection', close_limit=0.1):
        """
        Set the variable using another variable. Undefined values (e.g. outside
        the other mesh) are set to numpy.nan, or extrapolated.

        Parameters
        ----------
        strategy : 'projection' or 'interpolation'
            The strategy to set the values: the L^2 orthogonal projection (not
            implemented!), or a direct interpolation to the nodes (nodal
            elements only!)

        Notes
        -----
        If the other variable uses the same field mesh, the coefficients are
        set directly.
        """
        flag_same_mesh = self.has_same_mesh(other)

        if flag_same_mesh == 'same':
            self.set_data(other())
            return

        if strategy == 'interpolation':
            coors = self.get_interp_coors(strategy)

        elif strategy == 'projection':
            ## interp_term = Term() # TODO
            ## coors = self.get_interp_coors(strategy, interp_term)
            pass

        else:
            raise ValueError('unknown interpolation strategy! (%s)' % strategy)

        vals = other.evaluate_at(coors, strategy='general',
                                 close_limit=close_limit)

        if strategy == 'interpolation':
            self.set_data(vals)

        elif strategy == 'projection':
            raise NotImplementedError('unsupported strategy! (%s)' % strategy)

        else:
            raise ValueError('unknown interpolation strategy! (%s)' % strategy)


class DGFieldVariable(FieldVariable):
    """
    Fieald variable specificaly intended for use with DGFields, bypasses
    application of EBC and EPBC as this is done in DGField.

    Is instance checked in create_adof_conns.
    """

    def __init__(self, name, kind, field, order=None, primary_var_name=None,
                 special=None, flags=None, history=None, **kwargs):
        FieldVariable.__init__(self, name, kind, field, order=order,
                               primary_var_name=primary_var_name,
                               special=special, flags=flags,
                               history=history, **kwargs)

        from sfepy.discrete.dg.fields import DGField
        if isinstance(field, DGField):
            pass
        else:
            raise ValueError("Attempted to use DGFieldVariable with non DGField!")

    def apply_ebc(self, vec, offset=0, force_values=None):
        pass

    def get_full(self, r_vec, r_offset=0, force_value=None,
                 vec=None, offset=0):
        """
        Get the full DOF vector satisfying E(P)BCs from a reduced DOF
        vector.

        Notes
        -----
        The reduced vector starts in `r_vec` at `r_offset`.
        Passing a `force_value` overrides the EBC values. Optionally,
        `vec` argument can be provided to store the full vector (in
        place) starting at `offset`.
        """
        if vec is None:
            vec = nm.empty(self.n_dof, dtype=r_vec.dtype)

        else:
            vec = vec[offset:offset+self.n_dof]

        eq_map = self.eq_map
        r_vec = r_vec[r_offset:r_offset+eq_map.n_eq]

        # overide to hotfix second application of EBCs
        # # EBC.
        # vec[eq_map.eq_ebc] = get_default(force_value, eq_map.val_ebc)

        # Reduced vector values, for DG this is full vector as eq_map.eq
        # contains all dofs, cf. create_adof_conns
        vec[eq_map.eqi] = r_vec

        # EPBC.
        # vec[eq_map.master] = vec[eq_map.slave]

        unused_dofs = self.field.get('unused_dofs')
        if unused_dofs is not None:
            vec[:] = self.field.restore_substituted(vec)

        return vec
