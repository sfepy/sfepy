from collections import deque

from sfepy.base.base import *
import sfepy.base.la as la
from sfepy.fem.mesh import make_inverse_connectivity, find_nearest_nodes, \
     TreeItem
from sfepy.fem.integrals import Integral
from extmods.fem import raw_graph, evaluate_at
from sfepy.fem.utils import extend_cell_data
from sfepy.fem.conditions import Conditions
from sfepy.fem.dof_info \
     import DofInfo, EquationMap, LCBCOperators, \
            expand_nodes_to_equations, make_global_lcbc_operator

is_state = 0
is_virtual = 1
is_parameter = 2
is_field = 10

def create_dof_conn(conn, dpn):
    """Given element a node connectivity, create the dof connectivity."""
    if dpn == 1:
        dc = conn.copy()
    else:
        n_el, n_ep = conn.shape
        n_ed = n_ep * dpn
        dc = nm.empty( (n_el, n_ed), dtype = conn.dtype )
        for ic in range( n_ed ):
            inod = ic / dpn
            idof = ic % dpn
##                    iloc = ic
            iloc = n_ep * idof + inod # Hack: For DBD order.
            dc[:,iloc] = dpn * conn[:,inod] + idof

    return dc

def create_a_dof_conn(eq, dc, indx):
    """Given a dof connectivity and equation mapping, create the active dof
    connectivity."""
    aux = eq[dc]
    adc = aux + nm.asarray( nm.where( aux >= 0, indx.start, 0 ),
                            dtype = nm.int32 )
    return adc

##
# c: 26.07.2006, r: 15.04.2008
def zero_conf_ebc( conf ):
    new = {}
    for key, bcs in conf.iteritems():
        newbc = copy( bcs )
        newbc.dofs = {}
        for dd, val in bcs.dofs.iteritems():
            newbc.dofs[dd] = 0.0
        new[key] = newbc
    return new

def _fix_scalar_dc(dc1, dc2):
    aux = nm.empty((dc2.shape[0], 1), dtype=nm.int32)
    aux.fill(dc1)
    return aux

##
# 14.07.2006, c
class Variables( Container ):
    """
    Container holding instances of Variable.
    """

    @staticmethod
    def from_conf(conf, fields):
        Variable.reset()

        obj = Variables()
        for key, val in conf.iteritems():
            var = Variable.from_conf(key, val, fields)

            obj[var.name] = var

        # Move into equations.
        obj.setup_dtype()
        obj.setup_ordering()
        # Remove.
        obj.domain = fields[0].domain

        return obj

    def __init__(self):
        Container.__init__(self, OneTypeList(Variable),
                           state = set(),
                           virtual = set(),
                           parameter = set(),
                           has_virtual_dcs = False,
                           has_lcbc = False,
                           has_eq_map = False,
                           dual_map = None,
                           ordered_state = None,
                           ordered_virtual = None)

    def __setitem__(self, ii, var):
        Container.__setitem__(self, ii, var)

        if var.is_state():
            self.state.add(var.name)

        elif var.is_virtual():
            self.virtual.add(var.name)

        elif var.is_parameter():
            self.parameter.add(var.name)

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

        else:
            self.dtype = nm.complex128

    def link_duals(self):
        """
        Link state variables with corresponding virtual variables.

        Usually, when solving a PDE in the weak form, each state
        variable has to have a corresponding virtual variable.

        As it is possible to use state variables in place of parameter
        variables in term evaluation, the above requirement is relaxed -
        it is possible to define more state variables than virtual ones.
        """
        
        self.dual_map = {}
        for ii in self.virtual:
            vvar = self[ii]
            try:
                self[vvar.primary_var_name].dual_var_name = vvar.name
            except IndexError:
                msg = 'variable %s is not active!' % vvar.primary_var_name
                raise ValueError(msg)

            self.dual_map[vvar.name] = vvar.primary_var_name

	if ((len(self.dual_map) != len(set(self.dual_map.values())))
	    or (len(self.virtual) > len(self.state))):
                msg = 'state and virtual variables do not match!'
                raise ValueError(msg)

    def setup_ordering(self):
	"""
        Setup ordering of variables.
        """
        self.ordered_state = [None] * len(self.state)
        for var in self.iter_state(ordered=False):
            self.ordered_state[var._order] = var.name

        if self.dual_map is None:
            self.link_duals()

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

    def get_primary_names(self, var_names):
        """Return a dictionary of names of primary variables corresponding
        to the variables named in var_names."""
        primary_names = {}
        for name in var_names:
            var = self[name]
            if var.is_state():
                primary_names[name] = var.name
            else:
                primary_names[name] = var.primary_var_name

        return primary_names

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

    def setup_lcbc_operators(self, lcbc, regions):
        """
        Prepare linear combination BC operator matrix.
        """
        if lcbc is None: return

        self.lcbcs = Conditions.from_conf(lcbc)
        lcbc_of_vars = self.lcbcs.group_by_variables()

        # Assume disjoint regions.
        lcbc_ops = {}
        offset = 0
        for var_name, bcs in lcbc_of_vars.iteritems():
            var = self[var_name]

            lcbc_op = var.create_lcbc_operators( bcs, regions, offset )
            lcbc_ops[var_name] = lcbc_op

            if lcbc_op is not None:
                offset += lcbc_op.n_op

        self.op_lcbc = make_global_lcbc_operator(lcbc_ops, self.adi)
        self.has_lcbc = self.op_lcbc is not None

    ##
    # 04.10.2007, c
    def get_lcbc_operator( self ):
        if self.has_lcbc:
            return self.op_lcbc
        else:
            raise ValueError( 'no LCBC defined!' )

    def equation_mapping(self, ebc, epbc, regions, ts, functions,
                         vregions=None):

        if vregions is None:
            vregions = regions

        self.ebcs = Conditions.from_conf(ebc)
        self.epbcs = Conditions.from_conf(epbc)

        ##
        # Assing EBC, PBC to variables and regions.
        self.bc_of_vars = self.ebcs.group_by_variables()
        self.bc_of_vars = self.epbcs.group_by_variables(self.bc_of_vars)

        ##
        # List EBC nodes/dofs for each variable.
        for var_name in self.di.var_names:
            var = self[var_name]
            bcs = self.bc_of_vars.get(var.name, None)

            var_di = self.di.get_info(var_name)
            var.equation_mapping(bcs, regions, var_di, ts, functions)

            if self.has_virtual_dcs:
                vvar = self[var.dual_var_name]
                vvar_di = self.vdi.get_info(var_name)
                vvar.equation_mapping(bcs, vregions, vvar_di, ts, functions)

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

    def setup_initial_conditions(self, conf_ics, regions, functions):
        self.ics = Conditions.from_conf(conf_ics)
        self.ic_of_vars = self.ics.group_by_variables()

        for var_name in self.di.var_names:
            var = self[var_name]

            ics = self.ic_of_vars.get(var.name, None)
            if ics is None: continue

            var.setup_initial_conditions(ics, regions, self.di, functions)

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

    def get_mirror_region(self, region,
                          return_ig_map=False, return_ig_map_i=False):
        for info in iter_dict_of_lists(self.conn_info):
            if isinstance(region, str):
                out = info.mirror_map[region]
            else:
                out = info.mirror_map[region.name]

            if return_ig_map and return_ig_map_i:
                return out

            if return_ig_map:
                return out[:2]

            elif return_ig_map_i:
                return (out[0], out[2])

            else:
                return out[0]

    def setup_extra_data(self):
        """Dof connectivity key = (field.name, region.name, type, ig)"""
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
##             print key, ii
##             print info
            for var_name in info.all_vars:
                self[var_name].setup_extra_data(info.ps_tg, info, info.is_trace)
            
    def setup_dof_conns(self, make_virtual=False, single_term=False):
        """Dof connectivity key = (field.name, region.name, type, ig)"""
        output('setting up dof connectivities...')
        tt = time.clock()

        self.has_virtual_dcs = make_virtual == True

        dof_conns = {}
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
##             print key, ii
##             print info

            if info.primary is not None:
                var = self[info.primary]
                var.setup_extra_data(info.ps_tg, info, info.is_trace)
                var.setup_dof_conns(dof_conns, info.dc_type, info.get_region())

            if info.has_virtual and (ii == 0):
                # This is needed regardless make_virtual.
                var = self[info.virtual]
                var.setup_extra_data(info.v_tg, info, False)

                if make_virtual or single_term or (info.primary is None):
                    var.setup_dof_conns(dof_conns, info.dc_type,
                                        info.get_region(can_trace=False))

##         print dof_conns
##         pause()

        self.dof_conns = dof_conns
        output( '...done in %.2f s' % (time.clock() - tt) )

    def setup_a_dof_conns( self ):
        """Translate dofs to active dofs.
        Active dof connectivity key = (variable.name, region.name, type, ig)"""
        adof_conns = {}
        for key, dc in self.dof_conns.iteritems():
            for var in self.iter_state():
                if var.has_dof_conn_key(key):
                    indx = self.adi.indx[var.name]
                    eq = var.eq_map.eq
                    akey = (var.name,) + key[1:]
                    adof_conns[akey] = create_a_dof_conn(eq, dc, indx)

##         print adof_conns
##         pause()

        self.adof_conns = adof_conns

    def get_a_dof_conn(self, var_name, is_dual, dc_type, ig, is_trace=False):
        """Get active dof connectivity of a variable.
        
        Note that primary and dual variables must have same Region!"""
        if is_dual and not self.has_virtual_dcs:
            var_name = self.dual_map[var_name]

        if self[var_name].has_field:
            if not is_trace:
                region_name = dc_type.region_name
                aig = ig

            else:
                region, ig_map = self.get_mirror_region(dc_type.region_name,
                                                        return_ig_map_i=True)
                region_name = region.name
                aig = ig_map[ig]

        else: # scalar variable, assume no traces.
            if dc_type.type == 'scalar':
                region_name = aig = None

            else:
                region_name = dc_type.region_name
                aig = None

        key = (var_name, region_name, dc_type.type, aig)
        dc = self.adof_conns[key]
        return dc
        
    def create_matrix_graph( self, var_names = None, vvar_names = None ):
        """
        Create tangent matrix graph. Order of dof connectivities is not
        important here...
        """
        if not self.has_virtuals():
            output( 'no matrix (no test variables)!' )
            return None

        shape = (self.avdi.ptr[-1], self.adi.ptr[-1])
        output( 'matrix shape:', shape )
        if nm.prod( shape ) == 0:
            output( 'no matrix (zero size)!' )
            return None

        adcs = self.adof_conns

        # Only volume dof connectivities are used, with the exception of trace
        # surface dof connectivities.
        shared = set()
        rdcs = []
        cdcs = []
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
            dct = info.dc_type.type
            if not (dct in ('volume', 'scalar') or info.is_trace):
                continue

            rn, cn = info.virtual, info.state
            if (rn is None) or (cn is None):
                continue
            
            rreg_name = info.get_region_name(can_trace=False)
            creg_name = info.get_region_name()

            for rig, cig in info.iter_igs():
                rkey = (self.dual_map[rn], rreg_name, dct, rig)
                ckey = (cn, creg_name, dct, cig)
                
                dc_key = (rkey, ckey)
##                 print dc_key
                if not dc_key in shared:
                    try:
                        rdcs.append(adcs[rkey])
                        cdcs.append(adcs[ckey])
                    except:
                        debug()
                    shared.add(dc_key)

##         print shared
        for ii in range(len(rdcs)):
            if (rdcs[ii].ndim == 1) and (cdcs[ii].ndim == 2):
                rdcs[ii] = _fix_scalar_dc(rdcs[ii], cdcs[ii])

            elif (cdcs[ii].ndim == 1) and (rdcs[ii].ndim == 2):
                cdcs[ii] = _fix_scalar_dc(cdcs[ii], rdcs[ii])

            elif (cdcs[ii].ndim == 1) and (rdcs[ii].ndim == 1):
                rdcs[ii] = nm.array(rdcs[ii], ndmin=2)
                cdcs[ii] = nm.array(cdcs[ii], ndmin=2)

##             print rdcs[ii], cdcs[ii]
##         pause()

        if not shared:
            # No virtual, state variable -> no matrix.
            output( 'no matrix (empty dof connectivities)!' )
            return None
        
        output( 'assembling matrix graph...' )
        tt = time.clock()

#	shape = nm.array( shape, dtype = nm.long )
        ret, prow, icol = raw_graph( int( shape[0] ), int( shape[1] ),
                                    len( rdcs ), rdcs, cdcs )
        output( '...done in %.2f s' % (time.clock() - tt) )
        nnz = prow[-1]
        output( 'matrix structural nonzeros: %d (%.2e%% fill)' \
                % (nnz, float( nnz ) / nm.prod( shape ) ) )
##         print ret, prow, icol, nnz
	
        data = nm.zeros( (nnz,), dtype = self.dtype )
        matrix = sp.csr_matrix( (data, icol, prow), shape )
##         matrix.save( 'matrix', format = '%d %d %e\n' )
##         pause()

        return matrix

    def create_state_vector( self ):
        vec = nm.zeros( (self.di.ptr[-1],), dtype = self.dtype )
        return vec

    def create_stripped_state_vector( self ):
        vec = nm.zeros( (self.adi.ptr[-1],), dtype = self.dtype )
        return vec

    ##
    # 22.11.2005, c
    # 25.07.2006
    # 19.09.2006
    # 18.10.2006
    def apply_ebc( self, vec, force_values = None ):
        """Apply essential (Dirichlet) boundary conditions."""
        for var_name in self.bc_of_vars.iterkeys():
            eq_map = self[var_name].eq_map
            i0 = self.di.indx[var_name].start
            ii = i0 + eq_map.eq_ebc
##             print ii, eq_map.val_ebc
##             pause()
            if force_values is None:
                vec[ii] = eq_map.val_ebc
            else:
                if isinstance( force_values, dict ):
                    vec[ii] = force_values[var_name]
                else:
                    vec[ii] = force_values
            # EPBC.
            vec[i0+eq_map.master] = vec[i0+eq_map.slave]

    def apply_ic( self, vec, force_values = None ):
        """Apply initial conditions."""
        for var in self.iter_state():
            ii = self.di.indx[var.name]

            if force_values is None:
                vec[ii] = var.get_initial_condition()
            else:
                if isinstance( force_values, dict ):
                    vec[ii] = force_values[var_name]
                else:
                    vec[ii] = force_values

    ##
    # 27.11.2005, c
    # 09.12.2005
    # 25.07.2006
    # 18.10.2006
    def update_vec( self, vec, delta ):
        for var_name in self.di.var_names:
            eq_map = self[var_name].eq_map
            i0 = self.di.indx[var_name].start
            ii = i0 + eq_map.eqi
##            print ii.shape, delta[adi.indx[var_name]].shape
            vec[ii] -= delta[self.adi.indx[var_name]]
            # EPBC.
            vec[i0+eq_map.master] = vec[i0+eq_map.slave]

    def strip_state_vector( self, vec, follow_epbc = True ):
        """
        Strip a full vector by removing EBC dofs. If 'follow_epbc' is True,
        values of EPBC master dofs are not simply thrown away, but added to the
        corresponding slave dofs, just like when assembling.
        """
        svec = nm.empty( (self.adi.ptr[-1],), dtype = self.dtype )
        for var_name in self.di.var_names:
            eq_map = self[var_name].eq_map
            i0 = self.di.indx[var_name].start
            ii = i0 + eq_map.eqi
##            print ii.shape, delta[adi.indx[var_name]].shape
            aindx = self.adi.indx[var_name]
            svec[aindx] = vec[ii]

            if follow_epbc:
                """ In [10]: a
                    Out[10]: array([0, 1, 2, 3, 4])

                    In [11]: a[[0,0,0,1,1]] += 1

                    In [12]: a
                    Out[12]: array([1, 2, 2, 3, 4])
                    """
                # svec[aindx.start + eq_map.eq[eq_map.slave]] += vec[eq_map.master]
                for ii, im in enumerate( eq_map.master ):
                    i1 = aindx.start + eq_map.eq[eq_map.slave[ii]]
                    if i1 < 0: continue
                    svec[i1] += vec[im]
#                    print ii, i1, im, eq_map.slave[ii], svec[i1], vec[im]

        return svec

    def make_full_vec( self, svec, var_name = None, force_value = None ):
        """
        Make a full vector satisfying E(P)BC
        from a stripped vector. For a selected variable if var_name is set.
        """
        def _make_full_vec( vec, svec, eq_map ):
            # EBC.
            ii = eq_map.eq_ebc
            if force_value is None:
                vec[ii] = eq_map.val_ebc
            else:
                vec[ii] = force_value

            # Stripped vector values.
            ii = eq_map.eqi
            vec[ii] = svec

            # EPBC.
            vec[eq_map.master] = vec[eq_map.slave]

        if self.has_lcbc:
            svec = self.op_lcbc * svec

        if var_name is None:
            vec = self.create_state_vector()

            for var_name in self.di.var_names:
                eq_map = self[var_name].eq_map
                _make_full_vec( vec[self.di.indx[var_name]],
                                svec[self.adi.indx[var_name]], eq_map )
        else:
            vec = nm.empty( (self.di.n_dof[var_name],), dtype = self.dtype )
            eq_map = self[var_name].eq_map
            _make_full_vec( vec, svec, eq_map )

        return vec

    ##
    # 14.03.2007, c
    def has_ebc( self, vec, force_values = None ):
        for var_name in self.bc_of_vars.iterkeys():
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

    ##
    # 26.07.2006, c
    # 12.04.2007
    # 26.07.2007
    def get_state_part_view( self, state, var_name, stripped = False ):
        return state[self.get_indx( var_name, stripped )]

    ##
    # 26.07.2006, c
    # 12.04.2007
    # 26.07.2007
    def set_state_part( self, state, part, var_name, stripped = False ):
        state[self.get_indx( var_name, stripped )] = part


    def set_data(self, data, step = 0):
        """Set data (vectors of values) of variables.

        Arguments:
        data .. state vector or dictionary of {variable_name : data vector}
        step .. time history step, 0 = current.
        """
        if data is None: return
        
        if isinstance(data, dict):

            for key, val in data.iteritems():
                try:
                    var = self[key]
                except (ValueError, IndexError):
                    raise KeyError('unknown variable! (%s)' % key)

                var.data_from_any(val, step=step)

        elif isinstance(data, nm.ndarray):
            self.data_from_state(data)

        else:
            raise ValueError('unknown data class! (%s)' % data.__class__)

    ##
    # 24.07.2006, c
    # 25.07.2006
    # 04.08.2006
    def data_from_state( self, state = None ):
        for ii in self.state:
            var = self[ii]
            var.data_from_state( state, self.di.indx[var.name] )

    ##
    # 26.07.2006, c
    # 02.08.2006
    # 04.08.2006
    def non_state_data_from_state( self, var_names, state, var_names_state ):
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


    def state_to_output( self, vec, fill_value = None, var_info = None,
                       extend = True ):
        """Convert a state vector to a dictionary of output data usable by
        Mesh.write()."""

        n_nod = self.domain.shape.n_nod
        di = self.di

        if var_info is None:
            var_info = {}
            for name in di.var_names:
                var_info[name] = (False, name)

        out = {}
        for key, indx in di.indx.iteritems():
            var = self[key]

            if key not in var_info.keys(): continue
            is_part, name = var_info[key]

            details = di.details[key]
            dpn = details.dpn

            if is_part:
                aux = nm.reshape( vec, (di.n_dof[key] / dpn, dpn) )
            else:
                aux = nm.reshape( vec[indx], (di.n_dof[key] / dpn, dpn) )

            if var.field.approx_order != '0':
                # Has vertex data.
                if extend:
                    ext = var.extend_data( aux, n_nod, fill_value )
                else:
                    ext = var.remove_extra_data( aux )

                out[name] = Struct( name = 'output_data',
                                    mode = 'vertex', data = ext,
                                    var_name = key, dofs = var.dofs )
            else:
                if extend:
                    ext = extend_cell_data(aux, self.domain, var.field.region,
                                           val=fill_value)
                else:
                    ext = aux

                ext.shape = (ext.shape[0], 1, ext.shape[1], 1)
                out[name] = Struct( name = 'output_data',
                                    mode = 'cell', data = ext,
                                    var_name = key, dofs = var.dofs )

        out = self.convert_complex_output( out )
        
        return out

    def convert_complex_output( self, out_in ):
        out = {}
        for key, val in out_in.iteritems():

            if val.data.dtype in  complex_types:
                rval = copy( val )
                rval.data = val.data.real
                out['real(%s)' % key] = rval

                ival = copy( val )
                ival.data = val.data.imag
                out['imag(%s)' % key] = ival

            else:
                out[key] = val

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

    def init_state( self, state ):
        for var in self.iter_state():
            var.init_state( state, self.di.indx[var.name] )

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
    _all_vars = {}

    @staticmethod
    def reset():
        Variable._count = 0
        Variable._orders = []
        Variable._all_vars = {}

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

        self.data = deque()
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
        Variable._all_vars[name] = self

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
            obj.flags.add( is_other )
            msg = 'unknown variable kind: %s' % kind
            raise NotImplementedError( msg )
            
        self.kind = kind

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
    ##
    # 26.07.2007, c
    def is_other( self ):
        return is_other in self.flags

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

    def init_state( self, state, indx ):
        """Initialize data of variables with history."""
        if self.history is None: return

        self.data.append( None )
        self.step = 0
        self.data_from_state( state, indx, step = 0 )

    def time_update(self, ts, functions):
        """Implemented in subclasses."""
        pass

    def advance( self, ts ):
        if self.history is None: return

        self.step = ts.step + 1
        if self.history == 'previous':
            self.data.rotate()
        else:
            self.data.append( None )

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

    def __call__( self, step = 0, derivative = None, dt = None ):
        """Returns:
             if `derivative` is None: a view of the data vector,
             otherwise: required derivative of the data vector
             at time step given by `step`.

           Supports only the backward difference w.r.t. time."""
        if derivative is None:
            return self.data[step][self.indx]
        else:
            if self.history is None:
                msg = 'set history type of variable %s to use derivatives!'\
                      % self.name
                raise ValueError( msg )
            dt = get_default( dt, self.dt )
##            print self.name, step, dt
            return (self( step = step ) - self( step = step-1 )) / dt
            
    def get_initial_condition( self ):
        if self.initial_condition is None:
            return 0.0
        else:
            return self.initial_condition

    def get_full_state( self, step = 0 ):
        return self.data[step]

    def get_indx( self ):
        return self.indx

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
            print 'rcm', time.clock() - tt
            
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

        self.set_field(field)

        self.has_field = True
        self.has_bc = True

    def set_field(self, field):
        """
        Set field of the variable.

        Takes reference to a Field instance. Sets dtype according to
        field.dtype.
        """
        self.field = field
        self.n_nod = field.n_nod
        self.n_dof = self.n_nod * self.n_components

        if self.dof_name is None:
            dof_name = 'aux'
        else:
            dof_name = self.dof_name
        self.dofs = [dof_name + ('.%d' % ii) for ii in range(self.n_components)]

        self.flags.add(is_field)
        self.dtype = field.dtype

        self.current_ap = None

    def get_field(self):
        return self.field

    def get_dof_info(self, active=False):
        details = Struct(name = 'field_var_dof_details',
                         n_nod = self.n_nod,
                         dpn = self.n_components)
        if active:
            n_dof = self.n_adof

        else:
            n_dof = self.n_dof
            
        return n_dof, details

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type
        
        if geometry != None:
            geometry_flag = geometry.find('Surface') >= 0
        else:
            geometry_flag = False     
            
        if (dct == 'surface') or (geometry_flag):
            reg = info.get_region()
            reg.select_cells_of_surface(reset=False)

            self.field.aps.setup_surface_data(reg)

        elif dct == 'edge':
            raise NotImplementedError('dof connectivity type %s' % dct)

        elif dct == 'point':
            self.field.aps.setup_point_data(self.field, info.region)

        elif dct not in ('volume', 'scalar'):
            raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def setup_dof_conns(self, dof_conns, dc_type, region):
        """Setup dof connectivities of various kinds as needed by terms."""
        dpn = self.n_components
        field = self.field
        dct = dc_type.type

        ##
        # Expand nodes into dofs.
        can_point = True
        for region_name, ig, ap in field.aps.iter_aps(igs=region.igs):
            region_name = region.name # True region name.
            key = (field.name, region_name, dct, ig)
            if key in dof_conns: continue

            if dct == 'volume':
                dc = create_dof_conn(ap.econn, dpn)
                dof_conns[key] = dc

            elif dct == 'surface':
                sd = ap.surface_data[region_name]
                dc = create_dof_conn(sd.econn, dpn)
                dof_conns[key] = dc

            elif dct == 'edge':
                raise NotImplementedError('dof connectivity type %s' % dct)
                
            elif dct == 'point':
                if can_point:
                    # Point data only in the first group to avoid multiple
                    # assembling of nodes on group boundaries.
                    conn = ap.point_data[region_name]
                    dc = create_dof_conn(conn, dpn)
                    dof_conns[key] = dc
                    can_point = False

            else:
                raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def has_dof_conn_key(self, key):
        return self.field.name == key[0]

    def time_update(self, ts, functions):
        """
        Store time step, set variable data for variables with the setter
        function.
        """
        self.dt = ts.dt

        if hasattr(self, 'special') and ('setter' in self.special):
            setter_name = self.special['setter']
            setter = functions[setter_name]

            region = self.field.region
            nod_list = region.get_field_nodes(self.field, clean=True)
            nods = nm.unique1d(nm.hstack(nod_list))

            coor = self.field.get_coor(nods)
            self.data_from_any(setter(ts, coor, region=region))
            output('data of %s set by %s()' % (self.name, setter_name))

    def data_from_qp(self, data_qp, integral_name, step=0):
        """u_n = \sum_e (u_{e,avg} * volume_e) / \sum_e volume_e
               = \sum_e \int_{volume_e} u / \sum volume_e"""
        domain = self.field.domain
        if domain.shape.n_el != data_qp.shape[0]:
            msg = 'incomatible shape! (%d == %d)' % (domain.shape.n_el,
                                                     data_qp.shape[0])
            raise ValueError(msg)

        n_vertex = domain.shape.n_nod
        dim = data_qp.shape[2]

        nod_vol = nm.zeros((n_vertex,), dtype=nm.float64)
        data_vertex = nm.zeros((n_vertex, dim), dtype=nm.float64)
        for region_name, ig, ap in self.field.aps.iter_aps():
            ap_key = (integral_name, region_name, ig)
            aux, vg = self.get_approximation(ap_key, 'Volume')

            volume = nm.squeeze(vg.variable(2))
            iels = domain.regions[region_name].cells[ig]
            
            data_e = nm.zeros((volume.shape[0], 1, dim, 1), dtype=nm.float64)
            vg.integrate(data_e, data_qp[iels])

            ir = nm.arange(dim, dtype=nm.int32)

            conn = domain.groups[ig].conn
            for ii, cc in enumerate(conn):
                # Assumes unique nodes in cc!
                ind2, ind1 = nm.meshgrid(ir, cc)
                data_vertex[ind1,ind2] += data_e[iels[ii],0,:,0]
                nod_vol[cc] += volume[ii]
        data_vertex /= nod_vol[:,nm.newaxis]

        ##
        # Field nodes values - TODO!.
        #        data = self.field.interp_v_vals_to_n_vals(data_vertex)
        data = data_vertex.squeeze()
        self.indx = slice(0, len(data))

        self.data[step] = data

    def create_lcbc_operators(self, bcs, regions, offset):
        if len(bcs) == 0: return None

        bcs.canonize_dof_names(self.dofs)
        bcs.sort()

        ops = LCBCOperators('lcbc:%s' % self.name, self.eq_map, offset)
        for bc in bcs:
            output('lcbc:', self.name, bc.name)

            ops.add_from_bc(bc, self.field, regions)

        ops.finalize()
        
        return ops

    def equation_mapping(self, bcs, regions, var_di, ts, functions, warn=False):
        """Set n_adof."""
        
        self.eq_map = EquationMap('eq_map', self.dofs, var_di)
        if bcs is not None:
            bcs.canonize_dof_names(self.dofs)
            bcs.sort()

        self.eq_map.map_equations(bcs, self.field, regions, ts,
                                  functions, warn=warn)
        self.n_adof = self.eq_map.n_eq

    def setup_initial_conditions(self, ics, regions, di, functions, warn=False):
        """Setup of initial conditions."""
        ics.canonize_dof_names(self.dofs)
        ics.sort()

        for ic in ics:
            dofs, val = ic.dofs

            try:
                region = regions[ic.region]
            except IndexError:
                msg = "no region '%s' used in IC %s!" % (ic.region, ic)
                raise IndexError( msg )

            if warn:
                clean_msg = ('warning: ignoring nonexistent' \
                             ' IC node (%s) in ' % self.name)
            else:
                clean_msg = None

            nod_list = region.get_field_nodes(self.field, clean=True,
                                              warn=clean_msg)
            if len( nod_list ) == 0:
                continue

            vv = nm.empty( (0,), dtype = self.dtype )
            nods = nm.unique1d( nm.hstack( nod_list ) )
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

    def get_approximation( self, key, kind = 'Volume', is_trace = False ):
        return self.field.aps.get_approximation(key, kind, is_trace)

    ##
    # c: 28.11.2006, r: 15.01.2008
    def get_data_shapes( self, key, kind = 'Volume' ):
        iname, ig = key[0], key[-1]
        ap = self.field.aps.aps_per_group[ig]
        if kind == 'Volume':
            shape = ap.get_v_data_shape( iname )
        else:
            region_name = key[1]
            shape = ap.get_s_data_shape( iname, region_name )

        return shape

    def get_state_in_region( self, region, igs = None, reshape = True,
                             step = 0 ):
        nods = region.get_field_nodes( self.field, merge = True, igs = igs )
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

    def extend_data( self, data, n_nod, val = None ):
        """Extend data (with value val) to cover whole domain."""
        cnt_vn = self.field.cnt_vn
        indx = cnt_vn[cnt_vn >= 0]

        if val is None:
            if data.shape[1] > 1: # Vector.
                val = nm.amin( nm.abs( data ) )
            else: # Scalar.
                val = nm.amin( data )

        extdata = nm.empty( (n_nod, data.shape[1]), dtype = self.dtype )
        extdata.fill( val )
        extdata[indx] = data[:indx.size]

        return extdata

    ##
    # c: 12.05.2008, r: 12.05.2008
    def remove_extra_data( self, data ):
        """Removes data in extra nodes."""
        cnt_vn = self.field.cnt_vn
        indx = self.field.remap[cnt_vn[cnt_vn >= 0]]
        newdata = data[indx]
        return newdata

    def get_element_diameters(self, cells, mode, square=False):
        """Get diameters of selected elements."""
        field = self.field
        domain = field.domain

        cells = nm.array(cells)

        diameters = nm.empty((cells.shape[0],), dtype=nm.float64)

        igs = nm.unique1d(cells[:,0])
        for ig in igs:
            ap = field.aps.aps_per_group[ig]
            vg = ap.describe_geometry(field, 'Volume', field.region)

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

        out = {}
        if self.field.approx_order != '0':
            ext = self.extend_data(aux, n_nod, 0.0)

            out[self.name] = Struct(name = 'output_data',
                                    mode = 'vertex', data = ext,
                                    var_name = self.name, dofs = self.dofs)
        else:
            ext = extend_cell_data(aux, self.field.domain, self.field.region,
                                   val=0.0)

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
            coors = self.field.aps.get_physical_qps(region, integral)

        else:
            raise ValueError('unknown interpolation strategy! (%s)' % strategy)

        return coors

    def evaluate_at(self, coors, strategy='kdtree', flag_same_mesh='different',
                    close_limit=0.1, cache=None, ret_cells=False,
                    ret_status=False):
        """
        Evaluate self in the given physical coordinates.
        """
        # Assume different meshes -> general interpolation.
        mesh = self.field.create_mesh()
        scoors = mesh.coors

        output('interpolating from %d nodes to %d nodes...' % (scoors.shape[0],
                                                               coors.shape[0]))

        if cache is None:
            offsets, iconn = make_inverse_connectivity(mesh.conns, mesh.n_nod,
                                                       ret_offsets=True)
        else:
            offsets, iconn = cache.offsets, cache.iconn

        if strategy == 'kdtree':
            if cache is None:
                from scipy.spatial import cKDTree as KDTree
                ## from scipy.spatial import KDTree

                tt = time.clock()
                ctree = KDTree(scoors)
                output('ctree: %f s' % (time.clock()-tt))

            else:
                ctree = cache.ctree

            tt = time.clock()

            vals = nm.empty((coors.shape[0], self.n_components),
                            dtype=self.dtype)
            cells = nm.empty((coors.shape[0], 2), dtype=nm.int32)
            status = nm.empty((coors.shape[0],), dtype=nm.int32)
            source_vals = self()

            ics = ctree.query(coors)[1]
            ics = nm.asarray(ics, dtype=nm.int32)

            vertex_coorss, nodess, orders, mtx_is = [], [], [], []
            conns, conns0 = [], []
            for ap in self.field.aps:
                ps = ap.interp.poly_spaces['v']
                if ps.order == 0:
                    # Use geometry element space and connectivity to locate an
                    # element a point is in.
                    ps = ap.interp.gel.interp.poly_spaces['v']
                    assert_(ps.order == 1)

                    orders.append(0) # Important!
                    iels = ap.region.cells[ap.ig]
                    conn = ap.region.domain.groups[ap.ig].conn
                    conns.append(conn)

                else:
                    orders.append(ps.order)
                    conns.append(ap.econn)

                vertex_coorss.append(ps.geometry.coors)
                nodess.append(ps.nodes)
                mtx_is.append(ps.get_mtx_i())

                # Always the true connectivity for extracting source values.
                conns0.append(ap.econn)

            orders = nm.array(orders, dtype=nm.int32)

            evaluate_at(vals, cells, status, coors, source_vals,
                        ics, offsets, iconn,
                        scoors, conns0, conns,
                        vertex_coorss, nodess, orders, mtx_is,
                        1, close_limit, 1e-15, 100, 1e-8)

            output('interpolator: %f s' % (time.clock()-tt))

        elif strategy == 'crawl':
            raise NotImplementedError

        else:
            raise ValueError('unknown search strategy! (%s)' % strategy)

        output('...done')

        if ret_status:
            return vals, cells, status

        elif ret_cells:
            return vals, cells

        else:
            return vals

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
            interp_term = Term() # TODO
            coors = self.get_interp_coors(strategy, interp_term)

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
                                 flag_same_mesh=flag_same_mesh,
                                 close_limit=close_limit)
        
        if strategy == 'interpolation':
            self.data_from_any(vals)

        elif strategy == 'projection':
            self.data_from_projection(vals)

        else:
            raise ValueError('unknown interpolation strategy! (%s)' % strategy)

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

    def has_dof_conn_key(self, key):
        return self.name == key[0]
