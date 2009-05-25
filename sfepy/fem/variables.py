from collections import deque

from sfepy.base.base import *
import sfepy.base.la as la
from sfepy.fem.mesh import make_inverse_connectivity, find_nearest_nodes, \
     TreeItem
from sfepy.fem.integrals import Integral
from extmods.fem import raw_graph, inverse_element_mapping
from extmods.geometry import SurfaceGeometry

is_state = 0
is_virtual = 1
is_parameter = 2
is_other = 3
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

def create_lcbc_rigid( coors ):
    """Create transformation matrix for rigid LCBCs."""
    n_nod, dim = coors.shape

    mtx_e = nm.tile( nm.eye( dim, dtype = nm.float64 ), (n_nod, 1) )

    if dim == 2:
        mtx_r = nm.empty( (dim * n_nod, 1), dtype = nm.float64 )
        mtx_r[0::dim,0] = -coors[:,1]
        mtx_r[1::dim,0] = coors[:,0]
        n_rigid_dof = 3

    elif dim == 3:
        mtx_r = nm.zeros( (dim * n_nod, dim), dtype = nm.float64 )
        mtx_r[0::dim,1] = coors[:,2]
        mtx_r[0::dim,2] = -coors[:,1]
        mtx_r[1::dim,0] = -coors[:,2]
        mtx_r[1::dim,2] = coors[:,0]
        mtx_r[2::dim,0] = coors[:,1]
        mtx_r[2::dim,1] = -coors[:,0]
        n_rigid_dof = 6

    else:
        msg = 'dimension in [2,3]: %d' % dim
        raise ValueError( msg )

    return n_rigid_dof, nm.hstack( (mtx_r, mtx_e) )

def create_lcbc_no_penetration( normals ):
##     print normals
    ii = nm.abs( normals ).argmax( 1 )
##     print ii
    n_nod, dim = normals.shape

    irs = set( range( dim ) )

    data = []
    rows = []
    cols = []
    for idim in xrange( dim ):
        ic = nm.where( ii == idim )[0]
        if len( ic ) == 0: continue
##         print ic
##         print idim

        ir = list( irs.difference( [idim] ) )
        nn = nm.empty( (len( ic ), dim - 1), dtype = nm.float64 )
        for ik, il in enumerate( ir ):
            nn[:,ik] = - normals[ic,il] / normals[ic,idim]
        
        irn = dim * ic + idim
        ics = [(dim - 1) * ic + ik for ik in xrange( dim - 1 )]
        for ik in xrange( dim - 1 ):
            rows.append( irn )
            cols.append( ics[ik] )
            data.append( nn[:,ik] )

        ones = nm.ones( (nn.shape[0],), dtype = nm.float64 )
        for ik, il in enumerate( ir ):
            rows.append( dim * ic + il )
            cols.append( ics[ik] )
            data.append( ones )
            
##     print rows
##     print cols
##     print data
        
    rows = nm.concatenate( rows )
    cols = nm.concatenate( cols )
    data = nm.concatenate( data )

    n_np_dof = n_nod * (dim - 1)
    op = sp.coo_matrix( (data, (rows, cols)), shape = (n_nod * dim, n_np_dof) )
    op = op.tocsr()

##     import pylab
##     from sfepy.base.plotutils import spy
##     spy( op )
##     print op
##     pylab.show()

    return n_np_dof, op

def compute_nodal_normals( nodes, region, field ):
    """Nodal normals are computed by simple averaging of element normals of
    elements every node is contained in. """
    dim = field.dim[0]

    fa = region.domain.get_neighbour_lists( True )[2]
    region.setup_face_indices( fa )
    region.select_cells_of_surface()

    normals = nm.zeros( (nodes.shape[0], dim),
                        dtype = nm.float64 )
    mask = nm.zeros( (nodes.max()+1,), dtype = nm.int32 )
    imap = nm.empty_like( mask )
    imap.fill( nodes.shape[0] ) # out-of-range index for normals.
    imap[nodes] = nm.arange( nodes.shape[0], dtype = nm.int32 )
    
    for ig, fis in region.fis.iteritems():
        ap = field.aps[ig]
        n_fa = fis.shape[0]
        n_fp = ap.efaces.shape[1]
        face_type = 's%d' % n_fp

        faces = ap.efaces[fis[:,1]]
        ee = ap.econn[fis[:,0]]
        econn = nm.empty( faces.shape, dtype = nm.int32 )
        for ir, face in enumerate( faces ):
            econn[ir] = ee[ir,face]
        mask[econn] += 1

        integral = Integral( name = 'i', kind = 's',
                             quad_name = 'custom',
                             mode = 'custom' )
        integral.vals = ap.interp.nodes[face_type].bar_coors
        # Unit normals -> weights = ones.
        integral.weights = nm.ones( (n_fp,), dtype = nm.float64 )
        integral.setup()
        integral.create_qp()

        bf_sg, weights = ap.get_base( face_type, 1,
                                      integral = integral,
                                      base_only = False )

        sg = SurfaceGeometry( n_fa, n_fp, dim, n_fp )
        sg.describe( field.aps.coors, econn, bf_sg, weights )

        e_normals = sg.variable( 0 ).squeeze()

        # normals[imap[econn]] += e_normals
        im = imap[econn]
        for ii, en in enumerate( e_normals ):
            normals[im[ii]] += en

    # All nodes must have a normal.
    if not nm.all( mask[nodes] > 0 ):
        raise ValueError( 'region %s has not complete faces!' % region.name )

    normals /= la.norm_l2_along_axis( normals )[:,nm.newaxis]

    return normals

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
                         has_virtual_dcs = False,
                         has_lcbc = False,
                         has_eq_map = False )

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
        self.dual_map = {}
        for ii in self.virtual:
            vvar = self[ii]
            try:
                self[vvar.primary_var_name].dual_var_name = vvar.name
            except ValueError:
                msg = 'variable %s is not active!' % vvar.primary_var_name
                raise ValueError( msg )

            self.dual_map[vvar.name] = vvar.primary_var_name

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
        offset = 0
        for var_name, bcs in lcbc_of_vars.iteritems():
            var = self[var_name]
            lcbc_op = var.create_lcbc_operators( bcs, regions, offset )
            lcbc_ops[var_name] = lcbc_op
            if lcbc_op is not None:
                offset += lcbc_op.n_op

        n_dof = self.adi.ptr[-1]
        eq_lcbc = nm.zeros( (n_dof,), dtype = nm.int32 )
        n_dof_new = 0
        for var_name, lcbc_op in lcbc_ops.iteritems():
#            print var_name, lcbc_op
            if lcbc_op is None: continue

            indx = self.adi.indx[var_name]
            eq_lcbc[indx] = lcbc_op.eq_lcbc

            n_dof_new += nm.sum( lcbc_op.n_transformed_dof )

        if n_dof_new == 0:
            self.has_lcbc = False
            return
            
        ii = nm.nonzero( eq_lcbc )[0]
        n_constrained = ii.shape[0]
        n_dof_free = n_dof - n_constrained
        n_dof_reduced = n_dof_free + n_dof_new
        output( 'dofs: total %d, free %d, constrained %d, new %d'\
                % (n_dof, n_dof_free, n_constrained, n_dof_new) )
        output( ' -> reduced %d' % (n_dof_reduced) )
        mtx_lc = sp.lil_matrix( (n_dof, n_dof_reduced), dtype = nm.float64 )
        ir = nm.where( eq_lcbc == 0 )[0]
        ic = nm.arange( n_dof_reduced, dtype = nm.int32 )
        mtx_lc[ir,ic] = 1.0

        rows = []
        cols = []
        data = []
        for var_name, lcbc_op in lcbc_ops.iteritems():
            if lcbc_op is None: continue
            for ii, op_lc in enumerate( lcbc_op.ops_lc ):
                indx = nm.where( eq_lcbc == lcbc_op.markers[ii] )[0]
                icols = nm.arange(n_dof_free + lcbc_op.ics[ii],
                                  n_dof_free + lcbc_op.ics[ii+1])

                if isinstance(op_lc, sp.spmatrix):
                    lr, lc, lv = sp.find(op_lc)
                    rows.append(indx[lr])
                    cols.append(icols[lc])
                    data.append(lv)

                else:
                    irs, ics = nm.meshgrid(indx, icols)
                    rows.append(irs.ravel())
                    cols.append(ics.ravel())
                    data.append(op_lc.T.ravel())

        rows = nm.concatenate(rows)
        cols = nm.concatenate(cols)
        data = nm.concatenate(data)

        mtx_lc2 = sp.coo_matrix((data, (rows, cols)), shape=mtx_lc.shape)
        mtx_lc = (mtx_lc + mtx_lc2).tocsr()

##         import pylab
##         from sfepy.base.plotutils import spy
##         spy( mtx_lc )
##         print mtx_lc
##         pylab.show()

        self.op_lcbc = mtx_lc

    ##
    # 04.10.2007, c
    def get_lcbc_operator( self ):
        if self.has_lcbc:
            return self.op_lcbc
        else:
            raise ValueError( 'no LCBC defined!' )

    ##
    # c: 01.11.2005, r: 12.05.2008
    def equation_mapping( self, ebc, epbc, regions, ts, funmod,
                         vregions = None ):

        if vregions is None:
            vregions = regions

        ##
        # Assing EBC, PBC to variables and regions.
        self.bc_of_vars = self._list_bc_of_vars( ebc )
        dict_extend( self.bc_of_vars,
                     self._list_bc_of_vars( epbc, is_ebc = False ) )

        ##
        # List EBC nodes/dofs for each variable.
        for var_name, bcs in self.bc_of_vars.iteritems():
            var = self[var_name]
            var.equation_mapping( bcs, regions, self.di, ts, funmod )
            if self.has_virtual_dcs:
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
        if self.has_virtual_dcs:
            self.avdi = _create_a_dof_info( self.vdi )
        else:
            self.avdi = self.adi

        self.has_eq_map = True

    def setup_initial_conditions( self, conf_ics, regions, funmod ):
        self.ic_of_vars = self._list_bc_of_vars( conf_ics )

        for var_name, ics in self.ic_of_vars.iteritems():
            if len( ics ) == 0:
                continue
            var = self[var_name]
            var.setup_initial_conditions( ics, regions, self.di, funmod )

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

    def _setup_extra_data(self, var, geometry, info, is_trace, shared):
        dct = info.dc_type[0]

        if (dct == 'surface') or (geometry.find('Surface') >= 0):
            reg = info.get_region()
            if reg.name not in shared:
                shared.add(reg.name)

                fa = self.domain.get_neighbour_lists(force_faces=True)[2]
                reg.setup_face_indices(fa)
                reg.select_cells_of_surface()

            var.field.aps.setup_surface_data(reg)

        elif dct == 'edge':
            raise NotImplementedError('dof connectivity type %s' % dct)

        elif dct == 'point':
            var.field.aps.setup_point_data(var.field, info.region)

        elif dct != 'volume':
            raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def setup_extra_data(self):
        """Dof connectivity key = (field.name, region.name, type, ig)"""
        surface_regions = set()
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
##             print key, ii
##             print info
            for var_name in info.all_vars:
                self._setup_extra_data(self[var_name], info.ps_tg,
                                       info, info.is_trace, surface_regions)
            
    def setup_dof_conns(self, make_virtual=False, single_term=False):
        """Dof connectivity key = (field.name, region.name, type, ig)"""
        output('setting up dof connectivities...')
        tt = time.clock()

        self.has_virtual_dcs = make_virtual == True

        dof_conns = {}
        surface_regions = set()
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
##             print key, ii
##             print info

            if info.primary is not None:
                self._setup_extra_data(self[info.primary], info.ps_tg,
                                       info, info.is_trace, surface_regions)
                var = self[info.primary]
                var.setup_dof_conns(dof_conns, info.dc_type, info.get_region())

            if info.has_virtual and (ii == 0):
                # This is needed regardless make_virtual.
                self._setup_extra_data(self[info.virtual], info.v_tg, info,
                                       False, surface_regions)

                if make_virtual or single_term or (info.primary is None):
                    var = self[info.virtual]
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
                if var.field.name == key[0]:
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

        if not is_trace:
            key = (var_name, dc_type[1], dc_type[0], ig)

        else:
            region, ig_map = self.get_mirror_region(dc_type[1],
                                                    return_ig_map_i=True)
            key = (var_name, region.name, dc_type[0], ig_map[ig])

        dc = self.adof_conns[key]
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
        if not self.has_virtuals():
            output( 'no matrix!' )
            return None

        shape = (self.avdi.ptr[-1], self.adi.ptr[-1])
        output( 'matrix shape:', shape )
        if nm.prod( shape ) == 0:
            output( 'no matrix!' )
            return None

        adcs = self.adof_conns

        # Only volume dof connectivities are used, with the exception of trace
        # surface dof connectivities.
        shared = set()
        rdcs = []
        cdcs = []
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
            dct = info.dc_type[0]
            if not (dct == 'volume' or info.is_trace):
                continue

            rn, cn = info.virtual, info.state
            if (rn is None) or (cn is None):
                continue
            
            rreg = info.get_region(can_trace=False)
            creg = info.get_region()

            for rig, cig in info.iter_igs():
                rkey = (self.dual_map[rn], rreg.name, dct, rig)
                ckey = (cn, creg.name, dct, cig)
                
                dc_key = (rkey, ckey)
##                 print dc_key
                if not dc_key in shared:
                    try:
                        rdcs.append(adcs[rkey])
                        cdcs.append(adcs[ckey])
                    except:
                        debug()
                    shared.add(dc_key)
##         print rdcs
##         print cdcs
##         print shared

        if not shared:
            # No virtual, state variable -> no matrix.
            output( 'no matrix!' )
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

            for var_name in self.bc_of_vars.iterkeys():
                eq_map = self[var_name].eq_map
                _make_full_vec( vec[self.di.indx[var_name]],
                                svec[self.adi.indx[var_name]], eq_map )
        else:
            vec = nm.empty( (self.di.n_dofs[var_name],), dtype = self.dtype )
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
                raise RuntimeError( 'undefined variable %s' % var_name )
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

    def time_update( self, ts ):
        for var in self:
            var.time_update( ts )

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
        assert_( (history is None) or (history in ['previous', 'full']) )

        obj = Variable( flags, name = conf.name, key = key,
                        kind = kind, family = family, history = history )

        if kind == 'unknown':
            obj.flags.add( is_state )
            if hasattr( conf, 'order' ):
                obj._order = int( conf.order )
            else:
                msg = 'unnown variable %s: order missing' % conf.name
                raise ValueError( msg )
            obj.dof_name = obj.name
        elif kind == 'test':
            obj.flags.add( is_virtual )
            if hasattr( conf, 'dual' ):
                obj.primary_var_name = conf.dual
            else:
                msg = 'test variable %s: related unknown missing' % conf.name
                raise ValueError( msg )
            obj.dof_name = obj.primary_var_name
        elif kind == 'parameter':
            obj.flags.add( is_parameter )
            if hasattr( conf, 'like' ):
                obj.primary_var_name = conf.like
            else:
                msg = 'parameter variable %s: related unknown missing'\
                        % conf.name
                raise ValueError( msg )
            obj.dof_name = obj.primary_var_name
        else:
            obj.flags.add( is_other )
            msg = 'unknown variable family: %s' % family
            raise NotImplementedError( msg )

        if family == 'field':
            try:
                fld = fields[conf.field]
            except:
                msg = 'field "%s" does not exist!' % conf.field
                raise KeyError( msg )

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
        self.step = 0
        self.dt = 1.0
        self.initial_condition = None

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

    def time_update( self, ts ):
        self.dt = ts.dt

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

    def data_from_any( self, data = None, indx = None, step = 0 ):
        self.data[step] = data
        if indx is None:
            self.indx = slice( 0, len( data ) )
        else:
            self.indx = slice( int( indx.start ), int( indx.stop ) )
        self.n_dof = self.indx.stop - self.indx.start

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

    def set_field( self, field ):
        """Takes reference to a Field instance. Sets dtype according to
        field.dtype."""
        self.field = field
        self.dpn = nm.product( field.dim )
        self.n_nod = field.n_nod

        if self.dof_name is None:
            dof_name = 'aux'
        else:
            dof_name = self.dof_name
        self.dofs = [dof_name + ('.%d' % ii) for ii in range( self.dpn )]

        self.flags.add( is_field )
        self.dtype = field.dtype

    def setup_dof_conns(self, dof_conns, dc_type, region):
        """Setup dof connectivities of various kinds as needed by terms."""
        dpn = self.dpn
        field = self.field
        dct = dc_type[0]

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
    def expand_nodes_to_equations( self, nods, dofs=None ):
        """dofs must be already canonized - it is done in
        Variables._list_bc_of_vars()"""
        if dofs is None:
            dofs = self.dofs
            
        eq = nm.array( [], dtype = nm.int32 )
        for dof in dofs:
            idof = self.dofs.index( dof )
            eq = nm.concatenate( (eq, self.dpn * nods + idof) )
        return eq

    ##
    # c: 03.10.2007, r: 25.02.2008
    def create_lcbc_operators( self, bcs, regions, offset ):
        if len( bcs ) == 0: return None

        eq = self.eq_map.eq
        n_dof = self.eq_map.n_eq
        
        eq_lcbc = nm.zeros( (n_dof,), dtype = nm.int32 )
        
        ops_lc = []
        n_transformed_dof = []
        markers = []
        for ii, (key, bc) in enumerate( bcs ):
            print self.name, bc.name
            region = regions[bc.region]
            print region.name

            nmaster = region.get_field_nodes( self.field, merge = True )
            print nmaster.shape

            dofs, kind = bc.dofs
            meq = eq[self.expand_nodes_to_equations( nmaster, dofs )]
            assert_( nm.all( meq >= 0 ) )

            markers.append( offset + ii + 1 )
            eq_lcbc[meq] = markers[-1]
##             print meq, meq.shape
##             print nm.where( eq_lcbc )[0]
            
            if kind == 'rigid':
                mcoor = self.field.get_coor( nmaster )
                n_dof, op_lc = create_lcbc_rigid( mcoor )

                # Strip unconstrained dofs.
                n_nod, dim = mcoor.shape
                aux = dim * nm.arange(n_nod)
                indx = [aux + self.dofs.index(dof) for dof in dofs]
                indx = nm.array(indx).T.ravel()

                op_lc = op_lc[indx]

            elif kind == 'no_penetration':
                dim = self.field.dim[0]
                assert_( len( dofs ) == dim )

                normals = compute_nodal_normals( nmaster, region, self.field )

                if hasattr( bc, 'filename' ):
                    mesh = self.field.domain.mesh
                    nn = nm.zeros_like( mesh.coors )
                    nmax = region.all_vertices.shape[0]
                    nn[nmaster[:nmax]] = normals[:nmax]
                    out = {'normals' : Struct( name = 'output_data',
                                               mode = 'vertex', data = nn )}
                    mesh.write( bc.filename, out = out, io = 'auto' )

                n_dof, op_lc = create_lcbc_no_penetration( normals )

            else:
                raise ValueError( 'unknown LCBC kind! (%s)' % kind )

            n_transformed_dof.append( n_dof )
            ops_lc.append( op_lc )

        ics = nm.cumsum( nm.r_[0, n_transformed_dof]  )
        n_op = len( ops_lc )
        
        return Struct( eq_lcbc = eq_lcbc,
                       ops_lc = ops_lc,
                       n_op = n_op,
                       ics = ics,
                       markers = markers,
                       n_transformed_dof = n_transformed_dof )

    def clean_node_list( self, nod_list, ntype, region_name, warn = False ):
        for nods in nod_list[:]:
            if nods is None:
                nod_list.remove( nods )
                if warn:
                    output( 'warning: ignoring nonexistent %s'\
                            + ' node (%s) in %s'\
                            % (ntype, self.name, region.name) )
        return nod_list

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
                msg = "no region '%s' used in BC %s!" % (rname, bc)
                raise IndexError( msg )

##             print ir, key, bc
##             debug()
            # Get master region nodes.

            fn = region.get_field_nodes( field )
            master_nod_list = self.clean_node_list( fn, ntype, region.name,
                                                    warn = warn )
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
                fn = region.get_field_nodes( field )
                slave_nod_list = self.clean_node_list( fn, ntype, region.name,
                                                       warn = warn )
##                 print master_nod_list
##                 print slave_nod_list

                nmaster = nm.unique1d( nm.hstack( master_nod_list ) )
                nslave = nm.unique1d( nm.hstack( slave_nod_list ) )
##                 print nmaster + 1
##                 print nslave + 1
                if nmaster.shape != nslave.shape:
                    msg = 'EPBC list lengths do not match!\n(%s,\n %s)' %\
                          (nmaster, nslave)
                    raise ValueError( msg )

                mcoor = field.get_coor( nmaster )
                scoor = field.get_coor( nslave )
                fun = getattr( funmod, bc.match )
                i1, i2 = fun( mcoor, scoor )
##                print nm.c_[mcoor[i1], scoor[i2]]
##                print nm.c_[nmaster[i1], nslave[i2]] + 1

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

        assert_( (eq_map.eq_ebc.shape == eq_map.val_ebc.shape) )
##         print eq_map.eq_ebc.shape
##         pause()
        eq_map.eq[eq_map.eq_ebc] = -2
        eq_map.eq[eq_map.master] = -1
        eq_map.eqi = nm.compress( eq_map.eq >= 0, eq_map.eq )
        eq_map.eq[eq_map.eqi] = nm.arange( eq_map.eqi.shape[0],
                                           dtype = nm.int32 )
        eq_map.eq[eq_map.master] = eq_map.eq[eq_map.slave]
        eq_map.n_eq = eq_map.eqi.shape[0]
        eq_map.n_ebc = eq_map.eq_ebc.shape[0]
        eq_map.n_epbc = eq_map.master.shape[0]
##         print eq_map
##         pause()

    def setup_initial_conditions( self, ics, regions, di, funmod, warn = False ):
        """Setup of initial conditions."""
        for key, ic in ics:
            dofs, val = ic.dofs

            try:
                region = regions[ic.region]
            except IndexError:
                msg = "no region '%s' used in IC %s!" % (ic.region, ic)
                raise IndexError( msg )

            fn = region.get_field_nodes( self.field )
            nod_list = self.clean_node_list( fn, 'IC', region.name,
                                             warn = warn )
            if len( nod_list ) == 0:
                continue

            vv = nm.empty( (0,), dtype = self.dtype )
            nods = nm.unique1d( nm.hstack( nod_list ) )
            coor = self.field.get_coor( nods )
            if type( val ) == str:
                fun = getattr( funmod, val )
                vv = fun( ic, coor )
            else:
                vv = nm.repeat( [val], nods.shape[0] * len( dofs ) )

            eq = self.expand_nodes_to_equations( nods, dofs )

            ic_vec = nm.zeros( (di.n_dofs[self.name],), dtype = self.dtype )
            ic_vec[eq] = vv
            
            self.initial_condition = ic_vec

    def get_approximation( self, key, kind = 'Volume', is_trace = False ):
        aps = self.field.aps
        geometries = aps.geometries

        iname, region_name, ig = key

        if is_trace:
            g_key = (iname, kind, region_name, ig)
            try:
                ap, geometry = aps.geometries[g_key]
            except KeyError:
                msg = 'no trace geometry %s in %s' % (key, aps.geometries)
                raise KeyError( msg )

        else:
            ap = aps.aps_per_group[ig]
            g_key = (iname, kind, region_name, ap.name)

            try:
                geometry = geometries[g_key]
            except KeyError:
                msg = 'no geometry %s in %s' % (g_key, geometries)
                raise KeyError( msg )

        return ap, geometry

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

    def interp_to_points( self, points, mesh, ctree = None, iconn=None ):
        """
        Interpolate self into given points. Works for scalar variables only!
        """
        if iconn is None:
            iconn = make_inverse_connectivity(mesh.conns, mesh.n_nod,
                                              combine_groups=True)
        field = self.field
        vdim = field.dim[0]

        vals = nm.empty((vdim, len(points)), dtype=self.dtype)
        coor = mesh.coors
        conns = mesh.conns

        tts = [0.0, 0.0, 0.0, 0.0]
        tt0 = time.clock()
        for ii, point in enumerate(points):
##             print ii, point
            if ctree is None:
                tt = time.clock()
                ic = find_nearest_nodes(coor, point)
                tts[0] += time.clock() - tt
            else:
                tt = time.clock()
                if isinstance(ctree, TreeItem):
                    ic = ctree.find_nearest_node(coor, point)
                else:
                    ic = ctree.query(point)[1]
                tts[0] += time.clock() - tt

            els = iconn[ic]
            bf = None
            for ig, iel in els:
##                 print ic, ig, iel
                tt1 = time.clock()
                nodes = conns[ig][iel]
                ecoor = coor[nodes]
##                 print ecoor

                interp = field.aps[ig].interp
                ref_coors = interp.nodes['v'].bar_coors
                base_fun = interp.base_funs['v'].fun
                tts[2] += time.clock() - tt1

                tt = time.clock()
                n_v, dim = ecoor.shape
                if n_v == (dim + 1):
                    bc = la.barycentric_coors(point, ecoor)
                    xi = nm.dot(bc.T, ref_coors)
                else: # Tensor-product and other.
                    xi = nm.empty((ecoor.shape[1],), dtype=nm.float64)
                    inverse_element_mapping(xi, point, ecoor, ref_coors,
                                            100, 1e-8)
                tts[1] += time.clock() - tt

                try:
                    # Verify that we are inside the element.
                    bf = base_fun.value(nm.atleast_2d(xi), base_fun.nodes,
                                        suppress_errors=False)
                except AssertionError:
                    continue
                break
##             print xi, bf

            if bf is None:
                # Point outside the mesh.
                vals[:,ii] = nm.nan
            else:
                # For scalar fields only!!!
                vals[:,ii] = nm.dot(bf,self()[nodes])
        tts[-1] = time.clock() - tt0
        print tts
#        print tts[0], tts[3]
        
        return vals

## ##
## # 11.07.2006, c
## class FEVariable( Variable ):
##     """Finite element Variable
## field .. field description of variable (borrowed)
## """
