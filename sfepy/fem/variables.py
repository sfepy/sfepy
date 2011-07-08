import time

import numpy as nm

from sfepy.base.base import real_types, complex_types
from sfepy.base.base import assert_, get_default, get_default_attr
from sfepy.base.base import output, OneTypeList, Container, Struct
import sfepy.linalg as la
from sfepy.fem.meshio import convert_complex_output
from sfepy.fem.integrals import Integral
from sfepy.fem.dof_info \
     import DofInfo, EquationMap, LCBCOperators, \
            expand_nodes_to_equations, make_global_lcbc_operator, is_active_bc
from sfepy.fem.mappings import get_physical_qps
from sfepy.fem.evaluate_variable import eval_real, eval_real_extra, eval_complex

is_state = 0
is_virtual = 1
is_parameter = 2
is_field = 10

def create_adof_conn(eq, dc, indx):
    """Given a dof connectivity and equation mapping, create the active dof
    connectivity."""
    aux = eq[dc]
    adc = aux + nm.asarray( nm.where( aux >= 0, indx.start, 0 ),
                            dtype = nm.int32 )
    return adc

##
# 14.07.2006, c
class Variables( Container ):
    """
    Container holding instances of Variable.
    """

    @staticmethod
    def from_conf(conf, fields):
        """
        This method resets the variable counters for automatic order!
        """
        Variable.reset()

        obj = Variables()
        for key, val in conf.iteritems():
            var = Variable.from_conf(key, val, fields)

            obj[var.name] = var

        obj.setup_dtype()
        obj.setup_ordering()

        return obj

    def __init__(self, variables=None):
        Container.__init__(self, OneTypeList(Variable),
                           state = set(),
                           virtual = set(),
                           parameter = set(),
                           has_virtual_dcs = False,
                           has_lcbc = False,
                           has_eq_map = False,
                           ordered_state = [],
                           ordered_virtual = [])

        if variables is not None:
            for var in variables:
                self[var.name] = var

            self.setup_ordering()

        self.setup_dtype()

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

    def setup_dtype( self ):
        """Setup data types of state variables - all have to be of the same
        data type, one of nm.float64 or nm.complex128."""
        dtypes = {nm.complex128 : 0, nm.float64 : 0}
        for var in self.iter_state(ordered=False):
            dtypes[var.dtype] += 1

        if dtypes[nm.float64] and dtypes[nm.complex128]:
            raise ValueError( "All variables must have the same dtype!" )

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

    def setup_ordering(self):
        """
        Setup ordering of variables.
        """
        self.link_duals()

        orders = []
        for var in self:
            try:
                orders.append(var._order)
            except:
                pass
        orders.sort()

        self.ordered_state = [None] * len(self.state)
        for var in self.iter_state(ordered=False):
            ii = orders.index(var._order)
            self.ordered_state[ii] = var.name

        self.ordered_virtual = [None] * len(self.virtual)
        ii = 0
        for var in self.iter_state(ordered=False):
            if var.dual_var_name is not None:
                self.ordered_virtual[ii] = var.dual_var_name
                ii += 1

    ##
    # 26.07.2007, c
    def get_names( self, kind = None ):
        if kind is None:
            names = [var.name for var in self]
        else:
            names = [var.name for var in self if var.is_kind( kind )]
        return names

    def has_virtuals(self):
        return len(self.virtual) > 0

    def setup_dof_info(self, make_virtual=False):
        """
        Setup global DOF information.
        """
        self.di = DofInfo('state_dof_info')
        for var_name in self.ordered_state:
            self.di.append_variable(self[var_name])
        
        if make_virtual:
            self.vdi = DofInfo('virtual_dof_info')
            for var_name in self.ordered_virtual:
                self.vdi.append_variable(self[var_name])

        else:
            self.vdi = self.di

    def setup_lcbc_operators(self, lcbcs, ts=None, functions=None):
        """
        Prepare linear combination BC operator matrix.
        """
        if lcbcs is None:
            self.lcdi = self.adi
            return

        self.lcbcs = lcbcs
        lcbc_of_vars = self.lcbcs.group_by_variables()

        # Assume disjoint regions.
        lcbc_ops = {}
        offset = 0
        for var_name, bcs in lcbc_of_vars.iteritems():
            var = self[var_name]

            lcbc_op = var.create_lcbc_operators(bcs, offset,
                                                ts=ts, functions=functions)
            lcbc_ops[var_name] = lcbc_op

            if lcbc_op is not None:
                offset += lcbc_op.n_op

        self.op_lcbc, self.lcdi = make_global_lcbc_operator(lcbc_ops, self.adi)

        self.has_lcbc = self.op_lcbc is not None

    ##
    # 04.10.2007, c
    def get_lcbc_operator( self ):
        if self.has_lcbc:
            return self.op_lcbc
        else:
            raise ValueError( 'no LCBC defined!' )

    def equation_mapping(self, ebcs, epbcs, ts, functions, problem=None):
        """
        Create the mapping of active DOFs from/to all DOFs for all state
        variables.

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

            ## print var.eq_map
            ## pause()

        self.adi = DofInfo('active_state_dof_info')
        for var_name in self.ordered_state:
            self.adi.append_variable(self[var_name], active=True)

        if self.has_virtual_dcs:
            self.avdi = DofInfo('active_virtual_dof_info')
            for var_name in self.ordered_virtual:
                self.avdi.append_variable(self[var_name], active=True)

        else:
            self.avdi = self.adi

        self.has_eq_map = True

        return active_bcs

    def get_matrix_shape(self):
        if not self.has_eq_map:
            raise ValueError('call equation_mapping() first!')

        return (self.avdi.ptr[-1], self.adi.ptr[-1])

    def setup_initial_conditions(self, ics, functions):
        self.ics = ics
        self.ic_of_vars = self.ics.group_by_variables()

        for var_name in self.di.var_names:
            var = self[var_name]

            ics = self.ic_of_vars.get(var.name, None)
            if ics is None: continue

            var.setup_initial_conditions(ics, self.di, functions)

    ##
    # c: 09.01.2008, r: 09.01.2008
    def get_nodes_of_global_dofs( self, igdofs ):
        """not stripped..."""
        di = self.di
        
        nods = nm.empty( (0,), dtype = nm.int32 )
        for ii in self.state:
            var = self[ii]
            indx = di.indx[var.name]
            igdof = igdofs[(igdofs >= indx.start) & (igdofs < indx.stop)]
            ivdof = igdof - indx.start
            inod = ivdof / var.n_components
            nods = nm.concatenate( (nods, inod) )
##             print var.name, indx
##             print igdof
##             print ivdof
##             print inod
##             pause()
        return nods

    def setup_adof_conns( self ):
        """Translate dofs to active dofs.
        Active dof connectivity key = (variable.name, region.name, type, ig)"""
        self.adof_conns = {}
        for var in self:
            var.setup_adof_conns(self.adof_conns, self.adi)

    def create_state_vector( self ):
        vec = nm.zeros( (self.di.ptr[-1],), dtype = self.dtype )
        return vec

    def create_stripped_state_vector( self ):
        vec = nm.zeros( (self.adi.ptr[-1],), dtype = self.dtype )
        return vec

    def apply_ebc(self, vec, force_values=None):
        """
        Apply essential (Dirichlet) and periodic boundary conditions
        defined for the state variables to vector `vec`.
        """
        for var in self.iter_state():
            var.apply_ebc(vec, self.di.indx[var.name].start, force_values)

    def apply_ic(self, vec, force_values=None):
        """
        Apply initial conditions defined for the state variables to
        vector `vec`.
        """
        for var in self.iter_state():
            var.apply_ic(vec, self.di.indx[var.name].start, force_values)

    def strip_state_vector(self, vec, follow_epbc=True):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.

        If 'follow_epbc' is True, values of EPBC master dofs are not
        simply thrown away, but added to the corresponding slave dofs,
        just like when assembling.
        """
        svec = nm.empty((self.adi.ptr[-1],), dtype=self.dtype)
        for var in self.iter_state():
            aindx = self.adi.indx[var.name]
            svec[aindx] = var.get_reduced(vec, self.di.indx[var.name].start,
                                          follow_epbc)
        return svec

    def make_full_vec(self, svec, force_value=None):
        """
        Make a full DOF vector satisfying E(P)BCs from a reduced DOF
        vector.

        Passing a `force_value` overrides the EBC values.
        """
        self.check_vector_size(svec, stripped=True)

        if self.has_lcbc:
            svec = self.op_lcbc * svec

        vec = self.create_state_vector()
        for var in self.iter_state():
            indx = self.di.indx[var.name]
            aindx = self.adi.indx[var.name]
            var.get_full(svec, aindx.start, force_value, vec, indx.start)

        return vec

    ##
    # 14.03.2007, c
    def has_ebc( self, vec, force_values = None ):
        for var_name in self.di.var_names:
            eq_map = self[var_name].eq_map
            i0 = self.di.indx[var_name].start
            ii = i0 + eq_map.eq_ebc
            if force_values is None:
                if not nm.allclose( vec[ii], eq_map.val_ebc ):
                    return False
            else:
                if isinstance( force_values, dict ):
                    if not nm.allclose( vec[ii], force_values[var_name] ):
                        return False
                else:
                    if not nm.allclose( vec[ii], force_values ):
                        return False
            # EPBC.
            if not nm.allclose( vec[i0+eq_map.master], vec[i0+eq_map.slave] ):
                return False
        return True

    ##
    # 26.07.2007, c
    def get_indx( self, var_name, stripped = False, allow_dual = False ):
        var = self[var_name]

        if not var.is_state():
            if allow_dual and var.is_virtual():
                var_name = var.primary_var_name
            else:
                msg = '%s is not a state part' % var_name
                raise IndexError( msg )
        
        if stripped:
            return self.adi.indx[var_name]
        else:
            return self.di.indx[var_name]

    def check_vector_size(self, vec, stripped=False):
        """
        Check whether the shape of the DOF vector corresponds to the
        total number of DOFs of the state variables.

        Parameters
        ----------
        vec : array
            The vector of DOF values.
        stripped : bool
            If True, the size of the DOF vector should be reduced,
            i.e. without DOFs fixed by boundary conditions.
        """
        if not stripped:
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

    def get_state_part_view(self, state, var_name, stripped=False):
        self.check_vector_size(state, stripped=stripped)
        return state[self.get_indx( var_name, stripped )]

    def set_state_part(self, state, part, var_name, stripped=False):
        self.check_vector_size(state, stripped=stripped)
        state[self.get_indx( var_name, stripped )] = part

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
        if vec is not None:
            self.check_vector_size(vec)

        out = {}
        for var in self.iter_state():
            if vec is None:
                out[var.name] = var()

            else:
                out[var.name] = vec[self.di.indx[var.name]]

        return out

    def set_data(self, data, step=0, ignore_unknown=False):
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
        """
        if data is None: return

        if isinstance(data, dict):

            for key, val in data.iteritems():
                try:
                    var = self[key]

                except (ValueError, IndexError):
                    if ignore_unknown:
                        pass

                    else:
                        raise KeyError('unknown variable! (%s)' % key)

                else:
                    var.data_from_any(val, step=step)

        elif isinstance(data, nm.ndarray):
            self.data_from_state(data)

        else:
            raise ValueError('unknown data class! (%s)' % data.__class__)

    def data_from_state(self, state=None):
        self.check_vector_size(state)

        for ii in self.state:
            var = self[ii]
            var.data_from_state( state, self.di.indx[var.name] )

    def non_state_data_from_state(self, var_names, state, var_names_state):
        self.check_vector_size(state)

        if isinstance( var_names, str ):
            var_names = [var_names]
            var_names_state = [var_names_state]

        for ii, var_name in enumerate( var_names ):
            var_name_state = var_names_state[ii]
            if self[var_name_state].is_state():
                self[var_name].data_from_data( state,
                                               self.di.indx[var_name_state] )
            else:
                msg = '%s is not a state part' % var_name_state
                raise IndexError( msg )


    def state_to_output(self, vec, fill_value=None, var_info=None,
                        extend=True):
        """Convert a state vector to a dictionary of output data usable by
        Mesh.write()."""
        di = self.di

        if var_info is None:
            self.check_vector_size(vec)

            var_info = {}
            for name in di.var_names:
                var_info[name] = (False, name)

        out = {}
        for key, indx in di.indx.iteritems():
            var = self[key]

            if key not in var_info.keys(): continue
            is_part, name = var_info[key]

            if is_part:
                aux = vec

            else:
                aux = vec[indx]

            out.update(var.create_output(aux, extend=extend,
                                         fill_value=fill_value))

        return out

    ##
    # c: 27.11.2006, r: 22.05.2008
    def iter_state( self, ordered = True ):

        if ordered:
            for ii in self.ordered_state:
                yield self[ii]

        else:
            for ii in self.state:
                yield self[ii]

    def init_history(self):
        for var in self.iter_state():
            var.init_history()

    def time_update(self, ts, functions):
        output('updating variables...')
        for var in self:
            var.time_update(ts, functions)
        output('...done')

    def advance( self, ts ):
        for var in self.iter_state():
            var.advance( ts )


##
# 11.07.2006, c
class Variable( Struct ):
    _count = 0
    _orders = []
    _all_var_names = set()

    @staticmethod
    def reset():
        Variable._count = 0
        Variable._orders = []
        Variable._all_var_names = set()

    def from_conf(key, conf, fields):
        aux = conf.kind.split()
        if len(aux) == 2:
            kind, family = aux

        elif len(aux) == 3:
            kind, family = aux[0], '_'.join(aux[1:])

        else:
            raise ValueError('variable kind is 2 or 3 words! (%s)' % conf.kind)

        history = get_default_attr( conf, 'history', None )
        assert_( (history is None) or (history in ['previous', 'full']) )

        n_components = conf.get_default_attr('n_components', None)

        order = conf.get_default_attr('order', None)
        if order is not None:
            order = int(order)

        primary_var_name = conf.get_default_attr('dual', None)
        if primary_var_name is None:
            if hasattr(conf, 'like'):
                primary_var_name = get_default(conf.like, '(set-to-None)')

            else:
                primary_var_name = None

        special = conf.get_default_attr('special', None)

        if family == 'field':
            try:
                fld = fields[conf.field]
            except IndexError:
                msg = 'field "%s" does not exist!' % conf.field
                raise KeyError( msg )

            if n_components is None:
                # Workaround until new syntax for Variable is introduced.
                n_components = fld.shape[0]
                
            obj = FieldVariable(conf.name, kind, fld, n_components,
                                order, primary_var_name,
                                special=special, key=key, history=history)

        elif family == 'constant':
            obj = ConstantVariable(conf.name, kind, order, primary_var_name,
                                   conf.field, special=special,
                                   key=key, history=history)

        else:
            raise ValueError('unknown variable family! (%s)' % family)

        return obj
    from_conf = staticmethod( from_conf )

    def __init__(self, name, kind, n_components, order=None,
                 primary_var_name=None, special=None, flags=None, **kwargs):
        Struct.__init__(self, name=name, n_components=n_components,
                        **kwargs)

        self.flags = set()
        if flags is not None:
            for flag in flags:
                self.flags.add(flag)

        self.data = []
        self.data.append(None)
        self.indx = None
        self.n_dof = None
        self.step = 0
        self.dt = 1.0
        self.initial_condition = None
        self.dual_var_name = None

        if self.is_virtual():
            self.data = None

        self._set_kind(kind, order, primary_var_name, special=special)
        Variable._all_var_names.add(name)

    def _set_kind(self, kind, order, primary_var_name, special=None):
        if kind == 'unknown':
            self.flags.add(is_state)
            if order is not None:
                if order in Variable._orders:
                    raise ValueError('order %d already used!' % order)
                else:
                    self._order = order
                    Variable._orders.append(order)

            else:
                self._order = Variable._count
                Variable._orders.append(self._order)
            Variable._count += 1

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
            self.flags.add( is_parameter )
            msg = 'parameter variable %s: related unknown missing' % self.name
            self.primary_var_name = get_default(primary_var_name, None, msg)
            if self.primary_var_name == '(set-to-None)':
                self.primary_var_name = None
            self.dof_name = self.primary_var_name

            if special is not None:
                self.special = special

        else:
            raise NotImplementedError('unknown variable kind: %s' % kind)

        self.kind = kind

    def _setup_dofs(self, n_nod):
        """
        Setup number of DOFs and  DOF names.
        """
        self.n_nod = n_nod

        self.n_dof = self.n_nod * self.n_components

        if self.dof_name is None:
            dof_name = 'aux'
        else:
            dof_name = self.dof_name
        self.dofs = [dof_name + ('.%d' % ii) for ii in range(self.n_components)]

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
            if self._variables is not None:
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
            var = self._variables[self.dual_var_name]

        else:
            var = self._variables[self.primary_var_name]

        return var

    ##
    # 11.07.2006, c
    def is_state( self ):
        return is_state in self.flags

    ##
    # 11.07.2006, c
    def is_virtual( self ):
        return is_virtual in self.flags

    ##
    # 26.07.2007, c
    def is_parameter( self ):
        return is_parameter in self.flags

    def is_state_or_parameter( self ):
        return (is_state in self.flags) or (is_parameter in self.flags)

    ##
    # 26.07.2007, c
    def is_kind( self, kind ):
        return eval( 'self.is_%s()' % kind )

    ##
    # 26.07.2006, c
    def is_non_state_field( self ):
        return (is_field in self.flags)\
               and not (self.is_state() or self.is_virtual())

    def is_real( self ):
        return self.dtype in real_types

    def is_complex( self ):
        return self.dtype in complex_types

    def init_data(self, step=0):
        """
        Initialize the dof vector data of time step `step` to zeros.
        """
        if self.is_state_or_parameter():
            self.data[step] = nm.zeros((self.n_dof,), dtype=self.dtype)

    def set_constant(self, val):
        """
        Set the variable to a constant value.
        """
        data = nm.empty((self.n_dof,), dtype=self.dtype)
        data.fill(val)
        self.data_from_any(data)

    def get_primary_name(self):
        if self.is_state():
            name = self.name

        else:
            name = self.primary_var_name

        return name

    def init_history(self):
        """Initialize data of variables with history."""
        if self.history is None: return

        self.data.append( None )
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
        if self.history == 'previous':
            self.data[:] = [None, self.data[0].copy()]

            # Advance evaluate cache.
            for cache in self.evaluate_cache.itervalues():
                for key in cache.keys():
                    if key[4] == -1: # Previous time step.
                        key0 = list(key)
                        key0[4] = 0
                        key0 = tuple(key0)

                        if key0 in cache:
                            cache[key] = cache[key0]
                            cache.pop(key0)

                        else:
                            cache.pop(key)

        else:
            self.data.append(None)

    def data_from_state( self, state = None, indx = None, step = 0 ):
        """step: 0 = current,  """
        if (not self.is_state()) or (state is None): return

        self.data_from_any(state, indx, step)

    def data_from_data( self, data = None, indx = None, step = 0 ):
        if (not self.is_non_state_field()) or (data is None): return

        self.data_from_any(data, indx, step)

    def data_from_any(self, data=None, indx=None, step=0):
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

        else:
            self.data[step] = data
            self.indx = indx

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
            data = self.data[step]
            if data is None:
                if (self.step == 0) and (step == -1):
                    data = self.data[0]

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

    def get_initial_condition( self ):
        if self.initial_condition is None:
            return 0.0
        else:
            return self.initial_condition

    def get_full_state( self, step = 0 ):
        return self.data[step]

class CloseNodesIterator(Struct):

    def __init__(self, field, create_mesh=True, create_graph=True,
                 strategy=None):
        self.field = field
        self.coors = self.field.get_coor()

        if create_mesh or create_graph:
            self.mesh = self.field.create_mesh()

        if create_graph:
            self.graph = self.mesh.create_conn_graph()
            self.perm = self.get_permutation(strategy=strategy)
            self.strategy = strategy

        else:
            self.graph = None
            self.strategy = None

    def __call__(self, strategy=None):
        if strategy is None or (strategy != self.strategy):
            self.perm = self.get_permutation(strategy=strategy)
            self.strategy = strategy

        self.ii = 0
        return self

    def get_permutation(self, strategy=None):
        graph = self.graph

        n_nod = self.coors.shape[0]
        dtype = nm.int32

        ## tt = time.clock()

        if strategy is None:
            perm = nm.arange(n_nod, dtype=dtype)

        elif strategy == 'rcm':
            from sfepy.linalg import rcm
            perm = rcm(graph)
            ## print 'rcm', time.clock() - tt

        elif 'greedy' in strategy:
            ipop, iin = {'00' : (0, 0),
                         'e0' : (-1, 0),
                         '0e' : (0, -1),
                         'ee' : (-1, -1),
                         '01' : (0, 1),
                         }[strategy[-2:]]

            perm_i = nm.empty((n_nod,), dtype=dtype)
            perm_i.fill(-1)

            n_nod = perm_i.shape[0]
            num = graph.indptr[1:] - graph.indptr[:-1]

            ir = nm.argmin(num)
            perm_i[ir] = 0
            active = [ir]
            ii = 1
            while ii < n_nod:
                ir = active.pop(ipop)
                row = graph.indices[graph.indptr[ir]:graph.indptr[ir+1]]
##                 print ir, row
                ips = []
                for ip in row:
                    if perm_i[ip] < 0:
                        perm_i[ip] = ii
                        ii += 1
                        ips.append(ip)
                if iin >= 0:
                    active[iin:iin] = ips
                else:
                    active.extend(ips)

            perm = nm.empty_like(perm_i)
            perm[perm_i] = nm.arange(perm_i.shape[0], dtype=perm.dtype)

        ## print time.clock() - tt
             
        return perm

    def test_permutations(self, strategy='rcm'):
        from sfepy.linalg import permute_in_place, save_sparse_txt

        save_sparse_txt('graph', self.graph, fmt='%d %d %d\n')
        graph = self.graph.copy()

        perm = self.get_permutation('rcm')

        g_types = ['00', 'e0', '0e', 'ee', '01']
        g_names = ['greedy_%s' % ii for ii in g_types]
        g_perms = [self.get_permutation('greedy_%s' % ii) for ii in g_types]

        c1 = self.mesh.coors
        d1 = la.norm_l2_along_axis(c1[1:] - c1[:-1])
        d2 = la.norm_l2_along_axis(c1[perm][1:] - c1[perm][:-1])
        print d1.min(), d1.mean(), d1.max(), d1.std(), d1.var()
        print d2.min(), d2.mean(), d2.max(), d2.std(), d2.var()
        ds = []
        for g_perm in g_perms:
            d3 = la.norm_l2_along_axis(c1[g_perm][1:] - c1[g_perm][:-1])
            ds.append(d3)
            print d3.min(), d3.mean(), d3.max(), d3.std(), d3.var()

        permute_in_place(graph, perm)
        save_sparse_txt('graph_rcm', graph, fmt='%d %d %d\n')

        for ii, g_name in enumerate(g_names):
            graph = self.graph.copy()
            permute_in_place(graph, g_perms[ii])
            save_sparse_txt('graph_%s' % g_name, graph, fmt='%d %d %d\n')

        from matplotlib import pyplot as plt
        n_bins = 30
        plt.figure()
        plt.subplot(311)
        _, bins, ps = plt.hist(d1, n_bins, histtype='bar')
        plt.legend(ps[0:1], ['default'])
        plt.subplot(312)
        plt.hist(d2, bins, histtype='bar')
        plt.legend(ps[0:1], ['RCM'])
        plt.subplot(313)
        _, _, ps = plt.hist(nm.array(ds).T, bins, histtype='bar')
        plt.legend([ii[0] for ii in ps], g_names)
        plt.savefig('hist_distances_sub.pdf', transparent=True)

        plt.figure()
        _, _, ps = plt.hist(nm.array([d1, d2] + ds).T, n_bins, histtype='bar')
        plt.legend([ii[0] for ii in ps], ['default', 'RCM'] + g_names)
        plt.savefig('hist_distances.pdf', transparent=True)
        plt.show()

    def __iter__(self):
        return self

    def next(self):
        try:
            ii = self.perm[self.ii]
            val = self.coors[ii]
        except IndexError:
            raise StopIteration

        self.ii += 1

        return ii, val

## ##
## # 11.07.2006, c
## class FEVariable( Variable ):
##     """Finite element Variable
## field .. field description of variable (borrowed)
## """
class FieldVariable(Variable):
    """A finite element field variable.
    
    field .. field description of variable (borrowed)
    """

    def __init__(self, name, kind, field, n_components, order=None,
                 primary_var_name=None, special=None, flags=None, **kwargs):
        Variable.__init__(self, name, kind, n_components, order,
                          primary_var_name, special, flags, **kwargs)

        self._set_field(field)

        self.has_field = True
        self.has_bc = True
        self.has_lcbc = False
        self._variables = None

        self.clear_bases()
        self.clear_current_group()
        self.clear_evaluate_cache()

    def _set_field(self, field):
        """
        Set field of the variable.

        Takes reference to a Field instance. Sets dtype according to
        field.dtype. Sets `dim` attribute to spatial dimension.
        """
        from sfepy.fem.fields import SurfaceField

        if isinstance(field, SurfaceField):
            self.is_surface = True

        else:
            self.is_surface = False

        self.field = field
        self._setup_dofs(field.n_nod)

        self.flags.add(is_field)
        self.dtype = field.dtype

        self.dim = field.coors.shape[1]

    def get_field(self):
        return self.field

    def describe_geometry(self, geometry_type, region, integral, ig,
                          term_region=None):
        field = self.field

        if isinstance(region, str):
            region = field.region

        if term_region is None:
            term_region = region

        geo = field.describe_geometry(geometry_type, ig, region, term_region,
                                      integral)

        return geo

    def setup_adof_conns(self, adof_conns, adi):
        """
        Translate dof connectivity of the variable to active dofs.

        Active dof connectivity key:
            (variable.name, region.name, type, ig)
        """
        self.adof_conns = {}

        for key, dc in self.field.dof_conns.iteritems():
            var = self.get_primary()
            akey = (var.name,) + key[2:]
            if akey in adof_conns:
                self.adof_conns[akey] = adof_conns[akey]

            else:
                if var.name in adi.indx:
                    indx = adi.indx[var.name]
                    eq = var.eq_map.eq

                else: # Special or pure parameter variables.
                    indx = slice(0, var.n_dof)
                    eq = nm.arange(var.n_dof, dtype=nm.int32)

                self.adof_conns[akey] = create_adof_conn(eq, dc, indx)

        adof_conns.update(self.adof_conns)

    def get_global_node_tab(self, dc_type, ig, is_trace=False):

        if self.n_components == 1:

            if not is_trace:
                region_name = dc_type.region_name
                aig = ig

            else:
                aux = self.field.domain.regions[dc_type.region_name]
                region, _, ig_map = aux.get_mirror_region()
                region_name = region.name
                aig = ig_map[ig]

            key = (self.field.name, self.n_components, region_name,
                   dc_type.type, aig)
            dc = self.field.dof_conns[key]
            inod = self.field.get_vertices()
            nodtab = inod[dc];
        else:
            raise NotImplementedError

        return nodtab

    def get_dof_conn(self, dc_type, ig, active=False, is_trace=False):
        """Get active dof connectivity of a variable.
        
        Note that primary and dual variables must have same Region!"""
        if not active:
            dc = self.field.get_dof_conn(dc_type, ig)

        else:
            var = self.get_primary()

            if self.is_virtual():
                var_name = var.name

            else:
                var_name = self.name

        if not is_trace:
            region_name = dc_type.region_name
            aig = ig

        else:
            aux = self.field.domain.regions[dc_type.region_name]
            region, _, ig_map = aux.get_mirror_region()
            region_name = region.name
            aig = ig_map[ig]

        key = (var_name, region_name, dc_type.type, aig)
        dc = self.adof_conns[key]

        return dc

    def get_dof_info(self, active=False):
        details = Struct(name = 'field_var_dof_details',
                         n_nod = self.n_nod,
                         dpn = self.n_components)
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
            setter_name = self.special['setter']
            setter = functions[setter_name]

            region = self.field.region
            nod_list = self.field.get_dofs_in_region(region, clean=True)
            nods = nm.unique(nm.hstack(nod_list))

            coor = self.field.get_coor(nods)
            self.data_from_any(setter(ts, coor, region=region))
            output('data of %s set by %s()' % (self.name, setter_name))

    def data_from_qp(self, data_qp, integral, step=0):
        """
        Set DOFs of variable using values in quadrature points
        corresponding to the given integral.
        """
        data_vertex = self.field.average_qp_to_vertices(data_qp, integral)

        ##
        # Field nodes values - TODO!.
        #        data = self.field.interp_v_vals_to_n_vals(data_vertex)
        data = data_vertex.squeeze()
        self.indx = slice(0, len(data))

        self.data[step] = data

    def create_lcbc_operators(self, bcs, offset, ts=None, functions=None):
        if len(bcs) == 0: return None

        bcs.canonize_dof_names(self.dofs)
        bcs.sort()

        ops = LCBCOperators('lcbc:%s' % self.name, self.eq_map, offset)
        for bc in bcs:
            # Skip conditions that are not active in the current time.
            if not is_active_bc(bc, ts=ts, functions=functions):
                continue

            output('lcbc:', self.name, bc.name)

            ops.add_from_bc(bc, self.field)

        ops.finalize()

        self.has_lcbc = True

        return ops

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
        if bcs is not None:
            bcs.canonize_dof_names(self.dofs)
            bcs.sort()

        active_bcs = self.eq_map.map_equations(bcs, self.field, ts, functions,
                                               problem=problem, warn=warn)
        self.n_adof = self.eq_map.n_eq

        return active_bcs

    def setup_initial_conditions(self, ics, di, functions, warn=False):
        """Setup of initial conditions."""
        ics.canonize_dof_names(self.dofs)
        ics.sort()

        for ic in ics:
            region = ic.region
            dofs, val = ic.dofs

            if warn:
                clean_msg = ('warning: ignoring nonexistent' \
                             ' IC node (%s) in ' % self.name)
            else:
                clean_msg = None

            nod_list = self.field.get_dofs_in_region(region, clean=True,
                                                     warn=clean_msg)
            if len( nod_list ) == 0:
                continue

            vv = nm.empty( (0,), dtype = self.dtype )
            nods = nm.unique( nm.hstack( nod_list ) )
            coor = self.field.get_coor( nods )
            if type( val ) == str:
                fun = functions[val]
                vv = fun(coor, ic=ic)
            else:
                vv = nm.repeat( [val], nods.shape[0] * len( dofs ) )

            eq = expand_nodes_to_equations(nods, dofs, self.dofs)

            ic_vec = nm.zeros( (di.n_dof[self.name],), dtype = self.dtype )
            ic_vec[eq] = vv
            
            self.initial_condition = ic_vec

    def get_approximation(self, ig):
        return self.field.aps[ig]

    def assign_geometries(self, geometries):
        """
        Initialize the shared dict of geometries.
        """
        self.geometries = geometries

    def get_data_shape(self, ig, integral,
                       integration='volume', region_name=None):
        """
        Get element data dimensions for given approximation.

        Parameters
        ----------
        ig : int
            The element group index.
        integral : Integral instance
            The integral describing used numerical quadrature.
        integration : 'volume', 'surface', 'surface_extra' or 'point'
            The term integration type.
        region_name : str
            The name of surface region, required when `shape_kind` is
            'surface'.

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
        ap = self.field.aps[ig]
        if integration in ('surface', 'surface_extra'):
            data_shape = ap.get_s_data_shape(integral, region_name)

            if integration == 'surface_extra':
                n_en = ap.get_v_data_shape(integral)[-1]
                data_shape = data_shape[:-1] + (n_en,)

        elif integration == 'volume':
            data_shape = ap.get_v_data_shape(integral)

            # Override ap.region with the required region.
            region = self.field.domain.regions[region_name]
            data_shape = (region.get_n_cells(ig),) + data_shape[1:]

        elif integration == 'point':
            region = self.field.domain.regions[region_name]
            dofs = self.field.get_dofs_in_region(region, merge=True)
            data_shape = (dofs.shape[0], 0, 0, 1)

        else:
            raise NotImplementedError('unsupported integration! (%s)'
                                      % integration)

        data_shape += (self.n_components,)

        return data_shape

    def clear_bases(self):
        """
        Clear base functions, base function gradients and element data
        dimensions.
        """
        self.bfs = {}
        self.bfgs = {}
        self.data_shapes = {}

    def setup_bases(self, geo_key, ig, geo, integral, shape_kind='volume'):
        """
        Setup and cache base functions and base function gradients for
        given geometry. Also cache element data dimensions.
        """
        if geo_key not in self.bfs:
            ap = self.field.aps[ig]

            region_name = geo_key[1]

            self.data_shapes[geo_key] = self.get_data_shape(ig, integral,
                                                            shape_kind,
                                                            region_name)

            if shape_kind == 'surface':
                sd = ap.surface_data[region_name]
                key = sd.face_type

            elif shape_kind == 'volume':
                key = 'v'

            ebf = ap.get_base(key, 0, integral)
            bf = la.insert_strided_axis(ebf, 0, ap.econn.shape[0])
            self.bfs[geo_key] = bf

            if integral.kind == 'v':
                bfg = geo.variable(0)

            else:
                try:
                    bfg = geo.variable(3)

                except:
                    bfg = None

            self.bfgs[geo_key] = bfg

    def clear_current_group(self):
        """
        Clear current group data.
        """
        self._ap = None
        self._data_shape = None
        self._bf = self._bfg = None

    def set_current_group(self, geo_key, ig):
        """
        Set current group data, initialize current DOF counter to `None`.

        The current group data are the approximation, element data
        dimensions, base functions and base function gradients.
        """
        self._ap = self.field.aps[ig]
        self._data_shape = self.data_shapes[geo_key]
        self._bf = self.bfs[geo_key]
        self._bfg = self.bfgs[geo_key]

        self._idof = None
        self._inod = None
        self._ic = None

    def val(self, ic=None):
        """
        Return base function values in quadrature points.

        Parameters
        ----------
        ic : int, optional
            The index of variable component.
        """
        if self._inod is None:
            # Evaluation mode.
            out = self.val_qp(ic=ic)

        else:
            out = self._bf[..., self._inod : self._inod + 1]

        return out

    def val_qp(self, ic=None):
        """
        Return variable evaluated in quadrature points.

        Parameters
        ----------
        ic : int, optional
            The index of variable component.
        """
        vec = self()[ic::self.n_components]
        evec = vec[self._ap.econn]

        aux = la.insert_strided_axis(evec, 1, self._bf.shape[1])[..., None]

        out = la.dot_sequences(aux, self._bf, 'ATBT')

        return out

    def grad(self, ic=None, ider=None):
        """
        Return base function gradient (space elements) values in
        quadrature points.

        Parameters
        ----------
        ic : int, optional
            The index of variable component.
        ider : int, optional
            The spatial derivative index. If not given, the whole
            gradient is returned.
        """
        if ider is None:
            iders = slice(None)

        else:
            iders = slice(ider, ider + 1)

        if self._inod is None:
            out = self.grad_qp(ic=ic, ider=ider)

        else:
            out = self._bfg[..., iders, self._inod : self._inod + 1]

        return out

    def grad_qp(self, ic=None, ider=None):
        """
        Return variable gradient evaluated in quadrature points.

        Parameters
        ----------
        ic : int, optional
            The index of variable component.
        ider : int, optional
            The spatial derivative index. If not given, the whole
            gradient is returned.
        """
        if ider is None:
            iders = slice(None)

        else:
            iders = slice(ider, ider + 1)

        vec = self()[ic::self.n_components]
        evec = vec[self._ap.econn]

        aux = la.insert_strided_axis(evec, 1, self._bfg.shape[1])[..., None]
        out = la.dot_sequences(self._bfg[:, :, iders, :], aux)

        return out

    def iter_dofs(self):
        """
        Iterate over element DOFs (DOF by DOF).
        """
        n_en, n_c = self._data_shape[3:]

        for ii in xrange(n_en):
            self._inod = ii
            for ic in xrange(n_c):
                self._ic = ic
                self._idof = n_en * ic + ii
                yield self._idof

    def get_element_zeros(self):
        """
        Return array of zeros with correct shape and type for term
        evaluation.
        """
        n_el, n_qp = self._data_shape[:2]

        return nm.zeros((n_el, n_qp, 1,  1), dtype=self.dtype)

    def get_component_indices(self):
        """
        Return indices of variable components according to current term
        evaluation mode.

        Returns
        -------
        indx : list of tuples
            The list of `(ii, slice(ii, ii + 1))` of the variable
            components. The first item is the index itself, the second
            item is a convenience slice to index components of material
            parameters.
        """
        if self._ic is None:
            indx = [(ii, slice(ii, ii + 1)) for ii in range(self.n_components)]

        else:
            indx = [(ii, slice(ii, ii + 1)) for ii in [self._ic]]

        return indx

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
        for cache in self.evaluate_cache.itervalues():
            for key in cache.keys():
                if key[4] == step: # Given time step to clear.
                    cache.pop(key)

    def evaluate(self, ig, mode='val',
                 region=None, integral=None, integration=None,
                 step=0, time_derivative=None, is_trace=False,
                 dt=None):
        """
        Evaluate various quantities related to the variable according to
        `mode` in quadrature points defined by `integral`.

        The evaluated data are cached in the variable instance in
        `evaluate_cache` attribute.

        Parameters
        ----------
        ig : int
            The element group index.
        mode : one of 'val', 'grad', 'div', 'cauchy_strain'
            The evaluation mode.
        region : Region instance, optional
            The region where the evaluation occurs. If None, the
            underlying field region is used.
        integral : Integral instance, optional
            The integral defining quadrature points in which the
            evaluation occurs. If None, the first order volume integral
            is created. Must not be None for surface integrations.
        integration : one of 'volume', 'surface', 'surface_extra'
            The term integration type. If None, it is derived from
            `integral`.
        step : int, default 0
            The time step (0 means current, -1 previous, ...).
        derivative : None or 'dt'
            If not None, return time derivative of the data,
            approximated by the backward finite difference.
        is_trace : bool, default False
            Indicate evaluation of trace of the variable on a boundary
            region.
        dt : float, optional
            The time step to be used if `derivative` is `'dt'`. If None,
            the `dt` attribute of the variable is used.

        Returns
        -------
        out : array
            The 4-dimensional array of shape
            `(n_el, n_qp, n_row, n_col)` with the requested data,
            where `n_row`, `n_col` depend on `mode`.
        """
        cache = self.evaluate_cache.setdefault(mode, {})

        field = self.field
        if region is None:
            region = field.region

        if is_trace:
            region, ig_map, ig_map_i = region.get_mirror_region()
            ig = ig_map_i[ig]

        if region is not field.region:
            assert_(field.region.contains(region))

        if integral is None:
            if integration in ('surface', 'surface_extra'):
                msg = 'integral must be given for surface integration!'
                raise ValueError(msg)

            integral = Integral('aux_1', 'v', 1)

        if integration is None:
            integration = {'v' : 'volume', 's' : 'surface'}[integral.kind]

        geo, _, key = field.get_mapping(ig, region, integral, integration,
                                        return_key=True)
        key += (step, time_derivative)

        if key in cache:
            out = cache[key]

        else:
            vec = self(step=step, derivative=time_derivative, dt=dt)
            ap = field.aps[ig]
            conn = ap.get_connectivity(region, integration)

            shape = self.get_data_shape(ig, integral, integration, region.name)

            if self.dtype == nm.float64:
                if integration != 'surface_extra':
                    out = eval_real(vec, conn, geo, mode, shape)

                else:
                    out = eval_real_extra(vec, conn, geo, mode, shape)

            else:
                out = eval_complex(vec, conn, geo, mode, shape)

            cache[key] = out

        return out

    def get_state_in_region( self, region, igs = None, reshape = True,
                             step = 0 ):
        nods = self.field.get_dofs_in_region(region, merge=True, igs=igs)
##         print nods, len( nods )
##         pause()
        eq = nm.empty( (len( nods ) * self.n_components,), dtype = nm.int32 )
        for idof in range( self.n_components ):
            eq[idof::self.n_components] = self.n_components * nods \
                                          + idof + self.indx.start

        out = self.data[step][eq]
        if reshape:
            out.shape = (len( nods ), self.n_components)

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
            if isinstance(force_values, dict ):
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

    def get_reduced(self, vec, offset=0, follow_epbc=True):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.

        Notes
        -----
        The full vector starts in `vec` at `offset`. If 'follow_epbc' is
        True, values of EPBC master DOFs are not simply thrown away, but
        added to the corresponding slave DOFs, just like when
        assembling.
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

        return vec

    def extend_dofs(self, data, fill_value=None):
        """
        Extend DOFs to the whole domain using the `fill_value`, or the
        smallest value in `dofs` if `fill_value` is None.
        """
        return self.field.extend_dofs(data, fill_value=fill_value)

    def remove_extra_dofs(self, dofs):
        """
        Remove DOFs defined in higher order nodes (order > 1).
        """
        return self.field.remove_extra_dofs(dofs)


    def create_output(self, vec=None, key=None, extend=True, fill_value=None):
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
        """
        if vec is None:
            vec = self()

        key = get_default(key, self.name)

        aux = nm.reshape(vec,
                         (self.n_dof / self.n_components, self.n_components))

        if extend:
            ext = self.extend_dofs(aux, fill_value)

        else:
            ext = self.remove_extra_dofs(aux)

        out = {}

        if ext is not None:
            approx_order = self.field.get_output_approx_order()

            if approx_order != 0:
                # Has vertex data.
                out[key] = Struct(name='output_data', mode='vertex', data=ext,
                                  var_name=self.name, dofs=self.dofs)

            else:
                ext.shape = (ext.shape[0], 1, ext.shape[1], 1)
                out[key] = Struct(name='output_data', mode='cell', data=ext,
                                  var_name=self.name, dofs=self.dofs)

        out = convert_complex_output(out)

        return out

    def get_element_diameters(self, cells, mode, square=False):
        """Get diameters of selected elements."""
        field = self.field
        domain = field.domain

        cells = nm.array(cells)

        diameters = nm.empty((cells.shape[0],), dtype=nm.float64)

        igs = nm.unique(cells[:,0])
        for ig in igs:
            ap = field.aps[ig]
            vg = ap.describe_geometry(field, 'volume', field.region)

            ii = nm.where(cells[:,0] == ig)[0]
            aux = domain.get_element_diameters(ig, cells[ii,1].copy(), vg,
                                               mode, square=square)
            diameters[ii] = aux

        return diameters

    def save_as_mesh(self, filename):
        """
        Save the field mesh and the variable values into a file for
        visualization. Only the vertex values are stored.
        """
        mesh = self.field.create_mesh(extra_nodes=False)
        vec = self()

        n_nod, n_dof, dpn = mesh.n_nod, self.n_dof, self.n_components
        aux = nm.reshape(vec, (n_dof / dpn, dpn))

        ext = self.extend_dofs(aux, 0.0)

        out = {}
        if self.field.approx_order != 0:
            out[self.name] = Struct(name = 'output_data',
                                    mode = 'vertex', data = ext,
                                    var_name = self.name, dofs = self.dofs)
        else:
            ext.shape = (ext.shape[0], 1, ext.shape[1], 1)
            out[self.name] = Struct(name = 'output_data',
                                    mode = 'cell', data = ext,
                                    var_name = self.name, dofs = self.dofs)

        mesh.write(filename, io='auto', out=out)

    def set_from_mesh_vertices(self, data):
        """Set the variable using values at the mesh vertices."""
        ndata = self.field.interp_v_vals_to_n_vals(data)
        self.data_from_any(ndata)

##         print data.shape
##         print ndata.shape

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

    def evaluate_at(self, coors, strategy='kdtree',
                    close_limit=0.1, cache=None, ret_cells=False,
                    ret_status=False):
        """
        Evaluate self in the given physical coordinates. Convenience
        wrapper around :func:`Field.evaluate_at()`, see its docstring
        for more details.
        """
        source_vals = self().reshape((self.n_nod, self.n_components))
        out = self.field.evaluate_at(coors, source_vals, strategy=strategy,
                                     close_limit=close_limit, cache=cache,
                                     ret_cells=ret_cells, ret_status=ret_status)

        return out

    def set_from_other(self, other, strategy='projection',
                       search_strategy='kdtree', ordering_strategy='rcm',
                       close_limit=0.1):
        """
        Set the variable using another variable. Undefined values (e.g. outside
        the other mesh) are set to numpy.nan, or extrapolated.

        Parameters
        ----------
        strategy : 'projection' or 'interpolation'
            The strategy to set the values: the L^2 orthogonal projection, or
            a direct interpolation to the nodes (nodal elements only!)
        
        Notes
        -----
        If the other variable uses the same field mesh, the coefficients are
        set directly.
        
        If the other variable uses the same field mesh, only deformed slightly,
        it is advisable to provide directly the node ids as a hint where to
        start searching for a containing element; the order of nodes does not
        matter then.

        Otherwise (large deformation, unrelated meshes, ...) there are
        basically two ways:
        a) query each node (its coordinates) using a KDTree of the other nodes
        - this completely disregards the connectivity information;
        b) iterate the mesh nodes so that the subsequent ones are close to each
        other - then also the elements of the other mesh should be close to each
        other: the previous one can be used as a start for the directional
        neighbour element crawling to the target point.

        Not sure which way is faster, depends on implementation efficiency and
        the particular meshes.
        """
        flag_same_mesh = self.has_same_mesh(other)

        if flag_same_mesh == 'same':
            self.data_from_any(other())
            return

        if strategy == 'interpolation':
            coors = self.get_interp_coors(strategy)

        elif strategy == 'projection':
            ## interp_term = Term() # TODO
            ## coors = self.get_interp_coors(strategy, interp_term)
            pass

        else:
            raise ValueError('unknown interpolation strategy! (%s)' % strategy)

        if search_strategy == 'kdtree':
            tt = time.clock()
            iter_nodes = CloseNodesIterator(self.field, create_graph=False)
            output('iterator: %f s' % (time.clock()-tt))
            
        elif search_strategy == 'crawl':
            tt = time.clock()
            iter_nodes = CloseNodesIterator(self.field, strategy='rcm')
            output('iterator: %f s' % (time.clock()-tt))

            iter_nodes.test_permutations()

        else:
            raise ValueError('unknown search strategy! (%s)' % search_strategy)

        perm = iter_nodes.get_permutation(iter_nodes.strategy)

        vals = other.evaluate_at(coors[perm], strategy=search_strategy,
                                 close_limit=close_limit)

        if strategy == 'interpolation':
            self.data_from_any(vals)

        elif strategy == 'projection':
            self.data_from_projection(vals)

        else:
            raise ValueError('unknown interpolation strategy! (%s)' % strategy)

class MultiplierVariable(FieldVariable):
    """
    A multiplier variable.

    This class can represent Lagrange multipliers defined in nodes of a
    field. Boundary conditions cannot be applied.
    """
    def __init__(self, name, kind, field, n_components, order=None,
                 flags=None, **kwargs):
        Variable.__init__(self, name, kind, n_components, order,
                          flags=flags, **kwargs)

        self._set_field(field)

        self.has_field = False
        self.has_bc = False
        self._variables = None

    def _set_field(self, field):
        """
        Set field of the variable.

        Takes reference to a Field instance. Sets dtype according to
        field.dtype.
        """
        self.field = field
        self._setup_dofs(field.n_nod)

        self.flags.add(is_field)
        self.dtype = field.dtype

    def get_dof_info(self, active=False):
        details = Struct(name = 'multiplier_var_dof_details',
                         n_nod = self.n_nod,
                         dpn = self.n_components)

        return self.n_dof, details

    def equation_mapping(self, bcs, var_di, ts, functions, problem=None,
                         warn=False):
        """
        Trivial mapping (no boundary conditions). Sets n_adof.

        Returns
        -------
        active_bcs : set
            The empty set.
        """
        if bcs is not None:
            raise ValueError('MultiplierVariable cannot have BC!')

        self.eq_map = EquationMap('eq_map', self.dofs, var_di)
        self.eq_map.map_equations(bcs, self.field, ts, functions, warn=warn)
        self.n_adof = self.eq_map.n_eq

        return set()

    def setup_adof_conns(self, adof_conns, adi):
        """
        The multiplier variables have no connectivity, so do nothing. It
        is up to user to allocate the global matrix entries properly.
        """
        self.adof_conns = {}

class ConstantVariable(Variable):
    """A constant variable.
    """
    def __init__(self, name, kind, order, primary_var_name,
                 dtype, special=None, flags=None, **kwargs):
        Variable.__init__(self, name, kind, order, primary_var_name,
                          flags, **kwargs)

        dtypes = {'real' : nm.float64, 'complex' : nm.complex128}
        self.dtype = dtypes[dtype]

        self.n_dof = 1

        self.has_field = False
        self.has_bc = False

    def get_dof_info(self):
        details = Struct(name = 'constant_var_dof_details')
        return self.n_dof, details

    def setup_extra_data(self, geometry, info, is_trace):
        pass

    def setup_dof_conns(self, dof_conns, dc_type, region):
        dct = dc_type.type
        if region is not None:
            region_name = region.name
        else:
            region_name = None

        key = (self.name, region_name, dct, None)
        dof_conns[key] = nm.zeros((1,), dtype=nm.int32)
