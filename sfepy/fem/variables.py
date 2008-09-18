from collections import deque

from sfepy.base.base import *
import sfepy.base.la as la
from extmods.fem import raw_graph

is_state = 0
is_virtual = 1
is_parameter = 2
is_other = 3
is_field = 10

##
# 11.07.2007, c
# 19.07.2007
def create_dof_conn( conn, dpn, imax ):
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

    try:
        imax = max( imax, nm.amax( dc.ravel() ) )
    except: # empty dc (a non-first point dc - e.g. springs)
        pass

    return dc, imax

##
# c: 11.07.2007, r: 04.02.2008
def create_a_dof_conns( eq, iterator, indx ):
    adcs = {}
    for key, dc in iterator():
        if isinstance( dc, dict ):
            adcss = create_a_dof_conns( eq, dc.iteritems, indx )
            for subkey, subdc in adcss.iteritems():
                adcs[(key, subkey)] = subdc
        elif dc is not None:
            aux = eq[dc]
            adcs[key] = aux + nm.asarray( nm.where( aux >= 0, indx.start, 0 ),
					  dtype = nm.int32 )
        else:
            adcs[key] = None
    return adcs

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

##
# 27.11.2006, c
def invert_data_shapes( var_shapes ):
    shapes = {}
    for name, groups in var_shapes.iteritems():
        for ig, shape in groups.iteritems():
            shapes.setdefault( ig, {} )[name] = shape
##                 if not shapes.has_key( ig ):
##                     shapes[ig] = {}
##                 shapes[ig][name] = shape
    return shapes

##
# 15.03.2007, c
# 04.06.2007
def resolve_chains( master_slave, chains ):
    """Resolve EPBC chains - e.g. in corner nodes."""

    for chain in chains:
        slave = chain[-1]
        master_slave[chain[:-1]] = slave + 1
        master_slave[slave] = - chain[0] - 1 # Any of masters...

##
# 04.06.2007, c
def group_chains( chain_list ):

    chains = []
    while len( chain_list ):
        chain = set( chain_list.pop( 0 ) )
##         print ':', chain
        ii = 0
        while ii < len( chain_list ):
            c1 = sorted( chain_list[ii] )
#            print '--', ii, c1, chain
            is0 = c1[0] in chain
            is1 = c1[1] in chain

            if is0 and is1:
                chain_list.pop( ii )
            elif is0 or is1:
                chain.update( c1 )
                chain_list.pop( ii )
                ii = 0
            else:
                ii += 1
#            print ii, chain, chain_list
##         print '->', chain
##         print chain_list
##         pause()

        chains.append( list( chain ) )

#    print 'EPBC chain groups:', chains
    aux = {}
    for chain in chains:
        aux.setdefault( len( chain ), [0] )[0] += 1
#    print 'EPBC chain counts:', aux

    return chains

##
# 19.07.2006
class DofInfo( Struct ):
    pass

##
# 14.07.2006, c
class Variables( Container ):

    ##
    # c: 14.07.2006, r: 22.05.2008
    def from_conf( conf, fields ):
        objs = OneTypeList( Variable )
        state, virtual, parameter, other = [], [], [], []
        for ii, (key, val) in enumerate( conf.iteritems() ):
            var = Variable.from_conf( key, val, fields )
            if var.is_state():
                state.append( ii )
            elif var.is_virtual():
                virtual.append( ii )
            elif var.is_parameter():
                parameter.append( ii )
            elif var.is_other():
                other.append( ii )
            objs.append( var )

        obj = Variables( objs,
                         state = state,
                         virtual = virtual,
                         parameter = parameter,
                         other = other,
                         domain = fields[0].domain,
                         fields = fields,
                         has_virtual_d_cs = False,
                         has_lcbc = False )

        indx = nm.array( [var._order
                          for var in obj.iter_state( ordered = False )] )
        obj.ordered_state = [obj.state[ii] for ii in indx.argsort()]

        obj.link_duals()

        obj.ordered_virtual = []
        for var in obj.iter_state( ordered = True ):
            for ii in obj.virtual:
                if obj[ii].primary_var_name == var.name:
                    obj.ordered_virtual.append( ii )
                    break

        obj.setup_dtype()

        return obj
    from_conf = staticmethod( from_conf )

    ##
    # c: 05.09.2007, r: 10.04.2008
    def link_duals( self ):
        for ii in self.virtual:
            vvar = self[ii]
            try:
                self[vvar.primary_var_name].dual_var_name = vvar.name
            except ValueError:
                output( 'variable %s is not active!' % vvar.primary_var_name )
                raise

    def setup_dtype( self ):
        """Setup data types of state variables - all have to be of the same
        data type, one of nm.float64 or nm.complex128."""
        dtypes = {nm.complex128 : 0, nm.float64 : 0}
        for var in self.iter_state():
            dtypes[var.dtype] += 1

        if dtypes[nm.float64] and dtypes[nm.complex128]:
            raise ValueError( "All variables must have the same dtype!" )

        elif dtypes[nm.float64]:
            self.dtype = nm.float64

        else:
            self.dtype = nm.complex128

    ##
    # 26.07.2007, c
    def get_names( self, kind = None ):
        if kind is None:
            names = [var.name for var in self]
        else:
            names = [var.name for var in self if var.is_kind( kind )]
        return names

    ##
    # c: 07.10.2005, r: 22.05.2008
    def setup_dof_info( self, make_virtual = False ):
        """Sets also i_dof_max."""
        def _setup_dof_info( iterable ):
            ptr = [0]
            n_nod = []
            dpn = []
            vnames = []
            for ii in iterable:
                var = self[ii]
                var.i_dof_max = var.field.n_nod * var.dpn

                n_nod.append( var.field.n_nod )
                dpn.append( var.dpn )
                ptr.append( ptr[-1] + dpn[-1] * n_nod[-1] )
                vnames.append( var.name )

            di = DofInfo(
                name = 'dof_info',
                ptr = nm.array( ptr, dtype = nm.int32 ),
                n_nod = nm.array( n_nod, dtype = nm.int32 ),
                dpn = nm.array( dpn, dtype = nm.int32 ),
                vnames = vnames, indx = {}, n_dofs = {}
            )

            for ii, name in enumerate( di.vnames ):
                di.indx[name] = slice( int( di.ptr[ii] ), int( di.ptr[ii+1] ) )
                di.n_dofs[name] = di.ptr[ii+1] - di.ptr[ii]
            return di

        self.di = _setup_dof_info( self.ordered_state )
        if make_virtual:
            self.vdi = _setup_dof_info( self.ordered_virtual )
        else:
            self.vdi = self.di

    ##
    # c: 16.10.2006, r: 15.04.2008
    def _list_bc_of_vars( self, bc_defs, is_ebc = True ):

        bc_of_vars = dict_from_keys_init( (key for key in self.di.vnames), list )
        if bc_defs is None: return bc_of_vars

        for key, bc in bc_defs.iteritems():
##             print key
##             print bc
            for dofs, val in bc.dofs.iteritems():
                vname = dofs.split( '.' )[0]
                if bc_of_vars.has_key( vname ):
                    var = self[vname]
                    vbc = copy( bc )
                    if is_ebc:
                        vbc.dofs = (var._canonize( dofs ), val)
                    else:
                        vbc.dofs = (var._canonize( dofs ), var._canonize( val ))
                    bc_of_vars[vname].append( (key, vbc) )
                else:
                    output( 'BC: ignoring nonexistent dof(s) %s' % dofs )
        return bc_of_vars

    ##
    # c: 03.10.2007, r: 18.02.2008
    def setup_lcbc_operators( self, lcbc, regions ):
        if lcbc is None: return

        self.has_lcbc =True

        lcbc_of_vars = self._list_bc_of_vars( lcbc )

        # Assume disjoint regions.
        lcbc_ops = {}
        for var_name, bcs in lcbc_of_vars.iteritems():
            var = self[var_name]
            lcbc_ops[var_name] = var.create_lcbc_operators( bcs, regions )

        ops_lc = []
        eq_lcbc = nm.empty( (0,), dtype = nm.int32 )
        n_groups = 0
        for var_name, lcbc_op in lcbc_ops.iteritems():
            if lcbc_op is None: continue
#            print var_name, lcbc_op

            indx = self.adi.indx[var_name]
            aux = nm.where( lcbc_op.eq_lcbc >= 0, indx.start, 0 )
            eq_lcbc = nm.hstack( (eq_lcbc, lcbc_op.eq_lcbc + aux) )
            ops_lc.extend( lcbc_op.ops_lc )

            n_rigid_dof = lcbc_op.n_rigid_dof
            dim = lcbc_op.dim
            n_groups += lcbc_op.n_groups

        if n_groups == 0:
            self.has_lcbc = False
            return
            
        n_dof = self.adi.ptr[-1]

        ii = nm.nonzero( eq_lcbc )[0]
        n_constrained = ii.shape[0]
        n_dof_not_rigid = n_dof - n_constrained
        n_dof_reduced = n_dof_not_rigid + n_groups * n_rigid_dof
        print n_dof, n_dof_reduced, n_constrained, n_dof_not_rigid

        mtx_lc = sp.lil_matrix( (n_dof, n_dof_reduced), dtype = nm.float64 )
        ir = nm.where( eq_lcbc == 0 )[0]
        ic = nm.arange( n_dof_reduced, dtype = nm.int32 )
        mtx_lc[ir,ic] = 1.0
        for ii, op_lc in enumerate( ops_lc ):
            indx = nm.where( eq_lcbc == (ii + 1) )[0]
            icols = slice( n_dof_not_rigid + n_rigid_dof * ii,
                           n_dof_not_rigid + n_rigid_dof * (ii + 1) )
            mtx_lc[indx,icols] = op_lc

        mtx_lc = mtx_lc.tocsr()
##         import pylab
##         from sfepy.base.plotutils import spy
##         spy( mtx_lc )
##         pylab.show()
##         print mtx_lc
        nnz = n_dof - n_constrained + n_constrained * dim
        print nnz, mtx_lc.getnnz()
        assert nnz >= mtx_lc.getnnz()

        self.op_lcbc = mtx_lc

    ##
    # 04.10.2007, c
    def get_lcbc_operator( self ):
        if self.has_lcbc:
            return self.op_lcbc
        else:
            print 'no LCBC defined!'
            raise ValueError

    ##
    # c: 01.11.2005, r: 12.05.2008
    def equation_mapping( self, ebc, epbc, regions, ts, funmod,
                         vregions = None ):

        if vregions is None:
            vregions = regions

        ##
        # Assing EBC, PBC to variables and regions.
        self.bc_of_vars = self._list_bc_of_vars( ebc )
        dict_extend( self.bc_of_vars, self._list_bc_of_vars( epbc, is_ebc = False ) )

        ##
        # List EBC nodes/dofs for each variable.
        for var_name, bcs in self.bc_of_vars.iteritems():
            var = self[var_name]
            var.equation_mapping( bcs, regions, self.di, ts, funmod )
            if self.has_virtual_d_cs:
                vvar = self[var.dual_var_name]
                vvar.equation_mapping( bcs, vregions, self.vdi, ts, funmod )

##             print var.eq_map
##             pause()

        ##
        # Adjust by offsets - create active dof info.
        def _create_a_dof_info( di ):
            adi = DofInfo(
                name = 'active_dof_info',
                ptr = nm.array( di.ptr, dtype = nm.int32 ),
                n_nod = nm.array( di.n_nod, dtype = nm.int32 ),
                dpn = nm.array( di.dpn, dtype = nm.int32 ),
                vnames = di.vnames, indx = {}, n_dofs = {}
            )
            for ii, key in enumerate( adi.vnames ):
                adi.n_dofs[key] = self[key].eq_map.n_eq
                adi.ptr[ii+1] = adi.ptr[ii] + adi.n_dofs[key]
                adi.indx[key] = slice( int( adi.ptr[ii] ), int( adi.ptr[ii+1] ) )
            return adi

        self.adi = _create_a_dof_info( self.di )
        if self.has_virtual_d_cs:
            self.avdi = _create_a_dof_info( self.vdi )
        else:
            self.avdi = self.adi

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
            inod = ivdof / var.dpn
            nods = nm.concatenate( (nods, inod) )
##             print var.name, indx
##             print igdof
##             print ivdof
##             print inod
##             pause()
        return nods

    ##
    # c: 23.11.2005, r: 26.02.2008
    def setup_dof_conns( self, make_virtual = False ):
        output( 'setting up dof connectivities...' )
        tt = time.clock()

        for ii in self.state:
            var = self[ii]
            var.setup_dof_conns()

        if make_virtual:
            for ii in self.virtual:
                var = self[ii]
                var.setup_dof_conns()
            self.has_virtual_d_cs = True
        else:
            self.has_virtual_d_cs = False

        output( '...done in %.2f s' % (time.clock() - tt) )

    ##
    # 08.08.2006, c
    # 11.10.2006
    # 20.02.2007
    # 22.02.2007
    # 11.07.2007
    # 05.09.2007
    def setup_a_dof_conns( self ):
        """Translate dofs to active dofs."""
        def _setup_a_dof_conns( iterable, adi ):
            adof_conns = dict_from_keys_init( [self[ii].name for ii in iterable],
                                          Struct )
            for ii in iterable:
                var = self[ii]
                indx = adi.indx[var.name]
                eq = var.eq_map.eq
                adof_conns[var.name].name = 'adof_conns'
                it =  var.iter_dof_conns( 'volume' )
                adof_conns[var.name].volume_d_cs = create_a_dof_conns( eq, it, indx )
                it =  var.iter_dof_conns( 'surface' )
                adof_conns[var.name].surface_d_cs = create_a_dof_conns( eq, it, indx )
                it =  var.iter_dof_conns( 'edge' )
                adof_conns[var.name].edge_d_cs = create_a_dof_conns( eq, it, indx )
                it =  var.iter_dof_conns( 'point' )
                adof_conns[var.name].point_d_cs = create_a_dof_conns( eq, it, indx )
            return adof_conns

        self.adof_conns = _setup_a_dof_conns( self.state, self.adi )
        if self.has_virtual_d_cs:
            self.avdof_conns = _setup_a_dof_conns( self.virtual, self.avdi )
        else:
            self.avdof_conns = self.adof_conns

##         print self.adof_conns.values()[0]
##         pause()

    ##
    # c: 10.12.2007, r: 15.01.2008
    def get_a_dof_conn( self, var_name, is_dual, dc_type, ig ):
        """Note that primary and dual variables must have same Region!"""
        kind, region_name = dc_type

        var = self[var_name]
        if is_dual:
            if not self.has_virtual_d_cs:
                var_name = var.primary_var_name
            adcs = self.avdof_conns[var_name]
        else:
            adcs = self.adof_conns[var_name]

        if kind == 'volume':
            try:
                dc = adcs.volume_d_cs[ig]
            except:
                debug()
        else:
            if kind == 'surface':
                dcs = adcs.surface_d_cs
            elif kind == 'edge':
                dcs = adcs.edge_d_cs
            elif kind == 'point':
                dcs = adcs.point_d_cs
            else:
                print 'uknown dof connectivity kind:', kind
                raise ValueError
            dc = dcs[(ig, region_name)]
        return dc
        
    ##
    # 05.09.2007, c
    def fix_coarse_grid_a_dof_conns( self, iemaps, var_name ):
        """Volume only!"""
        dcs = self.adof_conns[var_name].volume_d_cs
        for ig, dc in dcs.iteritems():
            iemap = iemaps[ig]
            dcs[ig] = dc[iemap]

    def create_matrix_graph( self, var_names = None, vvar_names = None ):
        """
        Create tangent matrix graph. Order of dof connectivities is not
        important here...
        """
        def _prepare_dc_lists( adof_conns, var_names = None ):
            if var_names is None:
                var_names = adof_conns.iterkeys()

            gdcs = {}
            for var_name in var_names:
                adcs = adof_conns[var_name]
                for ig, dc in adcs.volume_d_cs.iteritems():
##                     print dc
                    gdcs.setdefault( ig, [] ).append( dc )
            return gdcs

        shape = (self.avdi.ptr[-1], self.adi.ptr[-1])
        output( 'matrix shape:', shape )
        if nm.prod( shape ) == 0:
            output( 'no matrix!' )
            return None

        cgdcs = _prepare_dc_lists( self.adof_conns, var_names )
##         print cgdcs
##         pause()
        if self.has_virtual_d_cs:
            rgdcs = _prepare_dc_lists( self.avdof_conns, vvar_names )
        else:
            rgdcs = cgdcs

        ##
        # Make all permutations per element group.
        rdcs = []
        cdcs = []
        for ig in rgdcs.iterkeys():
            rgdc, cgdc = rgdcs[ig], cgdcs[ig]
            for perm in la.cycle( [len( rgdc ), len( cgdc )] ):
                rdcs.append( rgdc[perm[0]] )
                cdcs.append( cgdc[perm[1]] )
#                print ' ', perm, '->', rdcs[-1].shape, cdcs[-1].shape

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

    ##
    # 24.07.2006, c
    # 25.07.2006
    # 04.08.2006
    def data_from_state( self, state = None ):
        for ii in self.state:
            var = self[ii]
            var.data_from_state( state, self.di.indx[var.name] )

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

    ##
    # 27.11.2005, c
    # 09.12.2005
    # 25.07.2006
    # 18.10.2006
    def update_vec( self, vec, delta ):
        for var_name in self.bc_of_vars.iterkeys():
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
        for var_name in self.bc_of_vars.iterkeys():
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
        def _make_full_vec( vec, eq_map ):
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

            for var_name in self.bc_of_vars.iterkeys():
                eq_map = self[var_name].eq_map
                _make_full_vec( vec[self.di.indx[var_name]], eq_map )
        else:
            vec = nm.empty( (self.di.n_dofs[var_name],), dtype = self.dtype )
            eq_map = self[var_name].eq_map
            _make_full_vec( vec, eq_map )

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
                output( '%s is not a state part' % var_name )
                raise IndexError
        
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
                output( '%s is not a state part' % var_name_state )
                raise IndexError


    def state_to_output( self, vec, fill_value = None, var_info = None,
                       extend = True ):
        """Works for vertex data only."""

        n_nod, di = self.domain.shape.n_nod, self.di

        if var_info is None:
            var_info = {}
            for name in di.vnames:
                var_info[name] = (False, name)

        out = {}
        for key, indx in di.indx.iteritems():
            if key not in var_info.keys(): continue
            is_part, name = var_info[key]
            
            dpn = di.dpn[di.vnames.index( key )]

            if is_part:
                aux = nm.reshape( vec, (di.n_dofs[key] / dpn, dpn) )
            else:
                aux = nm.reshape( vec[indx], (di.n_dofs[key] / dpn, dpn) )

            if extend:
                ext = self[key].extend_data( aux, n_nod, fill_value )
            else:
                ext = self[key].remove_extra_data( aux )
#            print ext.shape
#            pause()
            out[name] = Struct( name = 'output_data',
                                mode = 'vertex', data = ext,
                                var_name = key, dofs = self[key].dofs )

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
    # 24.08.2006, c
    # 20.09.2006
    def get_fields_of_vars( self, var_names ):
        field_names = {}
        for var_name in var_names:
            if not self.has_key( var_name ):
                raise RuntimeError, 'undefined variable %s' % var_name
            field_names[var_name] = self[var_name].field.name
        return field_names

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

    def advance( self, ts ):
        for var in self.iter_state():
            var.advance( ts )


##
# 11.07.2006, c
class Variable( Struct ):

    def from_conf( key, conf, fields ):
        flags = set()
        kind, family = conf.kind.split()

        history = get_default_attr( conf, 'history', None )
        assert (history is None) or (history in ['previous', 'full'])

        obj = Variable( flags, name = conf.name, key = key,
                        kind = kind, family = family, history = history )

        if kind == 'unknown':
            obj.flags.add( is_state )
            if hasattr( conf, 'order' ):
                obj._order = int( conf.order )
            else:
                output( 'unnown variable %s: order missing' % conf.name )
                raise ValueError
            obj.dof_name = obj.name
        elif kind == 'test':
            obj.flags.add( is_virtual )
            if hasattr( conf, 'dual' ):
                obj.primary_var_name = conf.dual
            else:
                output( 'test variable %s: related unknown missing' % conf.name )
                raise ValueError
            obj.dof_name = obj.primary_var_name
        elif kind == 'parameter':
            obj.flags.add( is_parameter )
            if hasattr( conf, 'like' ):
                obj.primary_var_name = conf.like
            else:
                output( 'parameter variable %s: related unknown missing'\
                        % conf.name )
                raise ValueError
            obj.dof_name = obj.primary_var_name
        else:
            obj.flags.add( is_other )
            print 'unknown variable family: %s' % family
            raise NotImplementedError

        if family == 'field':
            try:
                fld = fields[conf.field]
            except:
                output( 'field "%s" does not exist!' % conf.field )
                raise

            obj.set_field( fld )

        return obj
    from_conf = staticmethod( from_conf )

    def __init__( self, flags, data = None, indx = 0, **kwargs ):
        Struct.__init__( self, **kwargs )

        self.flags = set()
        for flag in flags:
            self.flags.add( flag )

        self.data = deque()
        self.data.append( data )
        self.indx = None
        self.n_dof = None
        self.current_ap = None
        self.step = None

        if self.is_virtual():
            self.data = None

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

        self.indx = slice( int( indx.start ), int( indx.stop ) )
        self.n_dof = indx.stop - indx.start
        self.data[step] = state

    def data_from_data( self, data = None, indx = None, step = 0 ):
        if (not self.is_non_state_field()) or (data is None): return

        self.data[step] = data
        if indx is None:
            self.indx = slice( 0, len( data ) )
        else:
            self.indx = slice( int( indx.start ), int( indx.stop ) )
        self.n_dof = self.indx.stop - self.indx.start

    def set_field( self, field ):
        """Takes reference to a Field instance. Sets dtype according to
        field.dtype."""
        self.field = field
        self.dpn = nm.product( field.dim )

        if self.dof_name is None:
            dof_name = 'aux'
        else:
            dof_name = self.dof_name
        self.dofs = [dof_name + ('.%d' % ii) for ii in range( self.dpn )]

        self.flags.add( is_field )
        self.dtype = field.dtype

    ##
    # c: 08.08.2006, r: 15.01.2008
    def setup_dof_conns( self ):
        dpn = self.dpn
        field = self.field

        dof_conns = Struct( name = 'dof_conns', volume_d_cs = {},
                           surface_d_cs = {}, edge_d_cs = {}, point_d_cs = {} )
        imax = -1
        ##
        # Expand nodes into dofs.
        for region_name, ig, ap in field.aps.iter_aps():

            dc, imax = create_dof_conn( ap.econn, dpn, imax )
            dof_conns.volume_d_cs[ig] = dc

            if ap.surface_data:
                dcs = {}
                for key, sd in ap.surface_data.iteritems():
                    dc, imax2 = create_dof_conn( sd.econn, dpn, 0 )
                    assert imax2 <= imax
                    dcs[key] = dc
                dof_conns.surface_d_cs[ig] = dcs

            else:
                dof_conns.surface_d_cs[ig] = None

            if ap.edge_data:
                raise NotImplementedError
            else:
                dof_conns.edge_d_cs[ig] = None

            if ap.point_data:
                dcs = {}
                for key, conn in ap.point_data.iteritems():
                    dc, imax2 = create_dof_conn( conn, dpn, 0 )
                    # imax2 can be greater than imax, as all spring nodes
                    # are assigned to the first group!!!
##                     print conn
##                     print dc
##                     print key, dc.shape
##                     pause()
                    dcs[key] = dc
                dof_conns.point_d_cs[ig] = dcs

            else:
                dof_conns.point_d_cs[ig] = None

##         print dof_conns
##         pause()
        self.dof_conns = dof_conns
        assert self.i_dof_max >= imax

    ##
    # 20.02.2007, c
    # 11.07.2007
    def iter_dof_conns( self, kind ):
        if kind == 'volume':
            return self.dof_conns.volume_d_cs.iteritems
        elif kind == 'surface':
            return self.dof_conns.surface_d_cs.iteritems
        elif kind == 'edge':
            return self.dof_conns.edge_d_cs.iteritems
        elif kind == 'point':
            return self.dof_conns.point_d_cs.iteritems
        else:
            print 'uknown dof connectivity kind:', kind
            raise ValueError

    ##
    # c: 25.02.2008, r: 25.02.2008
    def _canonize( self, dofs ):
        vname, dd = dofs.split( '.' )
        if dd == 'all':
            cdofs = self.dofs
        elif dd[0] == '[':
            cdofs = [vname + '.' + ii.strip()
                     for ii in dd[1:-1].split( ',' )]
        else:
            cdofs = [dofs]
        return cdofs

    ##
    # c: 18.10.2006, r: 15.04.2008
    def expand_nodes_to_equations( self, nods, dofs ):
        """dofs must be already canonized - it is done in
        Variables._list_bc_of_vars()"""
        eq = nm.array( [], dtype = nm.int32 )
        for dof in dofs:
            idof = self.dofs.index( dof )
            eq = nm.concatenate( (eq, self.dpn * nods + idof) )
        return eq

    ##
    # c: 03.10.2007, r: 25.02.2008
    def create_lcbc_operators( self, bcs, regions ):
        if len( bcs ) == 0: return None

        eq = self.eq_map.eq
        n_dof = self.eq_map.n_eq
        
        n_groups = len( bcs )
        eq_lcbc = nm.zeros( (n_dof,), dtype = nm.int32 )
        
        ops_lc = []
        for ii, (key, bc) in enumerate( bcs ):
            print self.name, bc.name

            region = regions[bc.region]
            print region.name

            nmaster = region.get_field_nodes( self.field, merge = True )
            print nmaster.shape

            dofs, kind = bc.dofs
            meq = eq[self.expand_nodes_to_equations( nmaster, dofs )]
            assert nm.all( meq >= 0 )
            
            eq_lcbc[meq] = ii + 1
##             print meq, meq.shape
##             print nm.where( eq_lcbc )[0]
            
            mcoor = self.field.get_coor( nmaster )[:,:-1]
            n_nod, dim = mcoor.shape

#            print mcoor, mcoor.shape

            mtx_e = nm.tile( nm.eye( dim, dtype = nm.float64 ), (n_nod, 1) )
            if dim == 2:
                mtx_r = nm.empty( (dim * n_nod, 1), dtype = nm.float64 )
                mtx_r[0::dim,0] = -mcoor[:,1]
                mtx_r[1::dim,0] = mcoor[:,0]
                n_rigid_dof = 3
            elif dim == 3:
                mtx_r = nm.zeros( (dim * n_nod, dim), dtype = nm.float64 )
                mtx_r[0::dim,1] = mcoor[:,2]
                mtx_r[0::dim,2] = -mcoor[:,1]
                mtx_r[1::dim,0] = -mcoor[:,2]
                mtx_r[1::dim,2] = mcoor[:,0]
                mtx_r[2::dim,0] = mcoor[:,1]
                mtx_r[2::dim,1] = -mcoor[:,0]
                n_rigid_dof = 6
            else:
                print 'dimension in [2,3]: %d' % dim
                raise ValueError

            op_lc = nm.hstack( (mtx_r, mtx_e) )
##             print op_lc, op_lc.shape

            ops_lc.append( op_lc )

        return Struct( eq_lcbc = eq_lcbc,
                       ops_lc = ops_lc,
                       n_groups = n_groups,
                       n_rigid_dof = n_rigid_dof,
                       dim = dim )

    def equation_mapping( self, bcs, regions, di, ts, funmod, warn = False ):
        """EPBC: master and slave dofs must belong to the same field (variables
        can differ, though)."""
        # Sort by ebc definition name.
        bcs.sort( cmp = lambda i1, i2: cmp( i1[0], i2[0] ) )
##         print bcs
##         pause()
        
        self.eq_map = eq_map = Struct()

        eq_map.eq = nm.arange( di.n_dofs[self.name], dtype = nm.int32 )
        eq_map.val_ebc = nm.empty( (0,), dtype = self.dtype )
        if len( bcs ) == 0:
            ##
            # No ebc for this field.
            eq_map.eqi = nm.arange( di.n_dofs[self.name], dtype = nm.int32 )
            eq_map.eq_ebc = nm.empty( (0,), dtype = nm.int32 )
            eq_map.n_eq = eq_map.eqi.shape[0]
            eq_map.n_ebc = eq_map.eq_ebc.shape[0]
            eq_map.master = nm.empty( (0,), dtype = nm.int32 )
            eq_map.slave = nm.empty( (0,), dtype = nm.int32 )
            return

        field = self.field

        eq_ebc = nm.zeros( (di.n_dofs[self.name],), dtype = nm.int32 )
        val_ebc = nm.zeros( (di.n_dofs[self.name],), dtype = self.dtype )
        master_slave = nm.zeros( (di.n_dofs[self.name],), dtype = nm.int32 )
        chains = []
        for key, bc in bcs:
            if key[:3] == 'ebc':
                ntype = 'EBC'
                rname = bc.region
            else:
                ntype = 'EPBC'
                rname = bc.region[0]

            try:
                region = regions[rname]
            except IndexError:
                print "no region '%s' used in BC %s!" % (rname, bc)
                raise

##             print ir, key, bc
##             debug()
            # Get master region nodes.
            master_nod_list = region.get_field_nodes( field )
            for master in master_nod_list[:]:
                if master is None:
                    master_nod_list.remove( master )
                    if warn:
                        output( 'warning: ignoring nonexistent %s'\
                                + ' node (%s) in %s'\
                                % (ntype, self.name, region.name) )

            if len( master_nod_list ) == 0:
                continue

            if ntype == 'EBC': # EBC.
                dofs, val = bc.dofs
                ##
                # Evaluate EBC values.
                vv = nm.empty( (0,), dtype = self.dtype )
                nods = nm.unique1d( nm.hstack( master_nod_list ) )
                coor = field.get_coor( nods )
                if type( val ) == str:
                    fun = getattr( funmod, val )
                    vv = fun( bc, ts, coor )
                else:
                    vv = nm.repeat( [val], nods.shape[0] * len( dofs ) )
##                 print nods
##                 print vv
                eq = self.expand_nodes_to_equations( nods, dofs )
                # Duplicates removed here...
                eq_ebc[eq] = 1
                if vv is not None: val_ebc[eq] = vv

            else: # EPBC.
                region = regions[bc.region[1]]
                slave_nod_list = region.get_field_nodes( field )
                for slave in slave_nod_list[:]:
                    if slave is None:
                        slave_nod_list.remove( slave )
                        if warn:
                            output( 'warning: ignoring nonexistent EPBC'\
                                    + ' slave node (%s) in %s'\
                                    % (self.name, region.name) )
                if len( slave_nod_list ) == 0:
                    continue

##                 print master_nod_list
##                 print slave_nod_list

                nmaster = nm.unique1d( nm.hstack( master_nod_list ) )
                nslave = nm.unique1d( nm.hstack( slave_nod_list ) )
##                 print nmaster + 1
##                 print nslave + 1
                if nmaster.shape != nslave.shape:
                    raise 'EPBC list lengths do not match!\n(%s,\n %s)' %\
                          (nmaster, nslave)

                mcoor = field.get_coor( nmaster )
                scoor = field.get_coor( nslave )
                fun = getattr( funmod, bc.match )
                i1, i2 = fun( mcoor, scoor )
##                 print nm.c_[mcoor[i1], scoor[i2]]
##                 print nm.c_[nmaster[i1], nslave[i2]] + 1

                meq = self.expand_nodes_to_equations( nmaster[i1], bc.dofs[0] )
                seq = self.expand_nodes_to_equations( nslave[i2], bc.dofs[1] )

                m_assigned = nm.where( master_slave[meq] != 0 )[0]
                s_assigned = nm.where( master_slave[seq] != 0 )[0]
                if m_assigned.size or s_assigned.size: # Chain EPBC.
##                     print m_assigned, meq[m_assigned]
##                     print s_assigned, seq[s_assigned]

                    aux = master_slave[meq[m_assigned]]
                    sgn = nm.sign( aux )
                    om_chain = zip( meq[m_assigned], (aux - sgn) * sgn )
##                     print om_chain
                    chains.extend( om_chain )

                    aux = master_slave[seq[s_assigned]]
                    sgn = nm.sign( aux )
                    os_chain = zip( seq[s_assigned], (aux - sgn) * sgn )
##                     print os_chain
                    chains.extend( os_chain )

                    m_chain = zip( meq[m_assigned], seq[m_assigned] )
##                     print m_chain
                    chains.extend( m_chain )

                    msd = nm.setdiff1d( s_assigned, m_assigned )
                    s_chain = zip( meq[msd], seq[msd] )
##                     print s_chain
                    chains.extend( s_chain )

                    msa = nm.union1d( m_assigned, s_assigned )
                    ii = nm.setdiff1d( nm.arange( meq.size ), msa )
                    master_slave[meq[ii]] = seq[ii] + 1
                    master_slave[seq[ii]] = - meq[ii] - 1
                else:
                    master_slave[meq] = seq + 1
                    master_slave[seq] = - meq - 1
##                 print 'ms', master_slave
##                 print chains

##         print master_slave
        chains = group_chains( chains )
        resolve_chains( master_slave, chains )

        ii = nm.argwhere( eq_ebc == 1 )
        eq_map.eq_ebc = nm.atleast_1d( ii.squeeze() )
        eq_map.val_ebc = nm.atleast_1d( val_ebc[ii].squeeze() )
        eq_map.master = nm.argwhere( master_slave > 0 ).squeeze()
        eq_map.slave = master_slave[eq_map.master] - 1
##         print eq_map.master
##         print eq_map.slave
##         pause()

        assert (eq_map.eq_ebc.shape == eq_map.val_ebc.shape)
##         print eq_map.eq_ebc.shape
##         pause()
        eq_map.eq[eq_map.eq_ebc] = -2
        eq_map.eq[eq_map.master] = -1
        eq_map.eqi = nm.compress( eq_map.eq >= 0, eq_map.eq )
        eq_map.eq[eq_map.eqi] = nm.arange( eq_map.eqi.shape[0], dtype = nm.int32 )
        eq_map.eq[eq_map.master] = eq_map.eq[eq_map.slave]
        eq_map.n_eq = eq_map.eqi.shape[0]
        eq_map.n_ebc = eq_map.eq_ebc.shape[0]
        eq_map.n_epbc = eq_map.master.shape[0]
##         print eq_map
##         pause()

    ##
    # c: 24.07.2006, r: 15.01.2008
    def get_approximation( self, key, kind = 'Volume' ):
        iname, region_name, ig = key
##         print tregion_name, aregion_name, ig
#        print self.field.aps.aps_per_group

        aps = self.field.aps
        geometries = aps.geometries
        ap = aps.aps_per_group[ig]
        g_key = (iname, kind, region_name, ap.name)
        try:
            return ap, geometries[g_key]
        except KeyError:
            print g_key
            print geometries
            raise

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

    def __call__( self, step = 0 ):
        """Returns: a view of the state vector."""
        return self.data[step][self.indx]

    def get_full_state( self, step = 0 ):
        return self.data[step]

    def get_indx( self ):
        return self.indx

    def get_state_in_region( self, region, igs = None, reshape = True,
                             step = 0 ):
        nods = region.get_field_nodes( self.field, merge = True, igs = igs )
##         print nods, len( nods )
##         pause()
        eq = nm.empty( (len( nods ) * self.dpn,), dtype = nm.int32 )
        for idof in range( self.dpn ):
            eq[idof::self.dpn] = self.dpn * nods + idof + self.indx.start

        out = self.data[step][eq]
        if reshape:
            out.shape = (len( nods ), self.dpn)

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
        
## ##
## # 11.07.2006, c
## class FEVariable( Variable ):
##     """Finite element Variable
## field .. field description of variable (borrowed)
## """
