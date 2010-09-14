from sfepy.base.base import *
import sfepy.base.la as la
import fea
from mesh import Mesh, make_point_cells
import sfepy.terms as terms
import extmods.geometry as gm

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

def _fix_scalar_dc(dc1, dc2):
    aux = nm.empty((dc2.shape[0], 1), dtype=nm.int32)
    aux.fill(dc1)
    return aux

def fields_from_conf(conf, regions):
    fields = {}
    for key, val in conf.iteritems():
        field = Field.from_conf(val, regions)
        fields[field.name] = field

    return fields

def setup_extra_data(conn_info):
    """
    Setup extra data required for non-volume integration.
    """
    for key, ii, info in iter_dict_of_lists(conn_info, return_keys=True):
        ## print key, ii
        ## print info
        for var in info.all_vars:
            field = var.get_field()
            field.setup_extra_data(info.ps_tg, info, info.is_trace)

def setup_dof_conns(conn_info, dof_conns=None,
                    make_virtual=False, verbose=True):
    """
    Dof connectivity key:
        (field.name, var.n_components, region.name, type, ig)
    """
    if verbose:
        output('setting up dof connectivities...')
        tt = time.clock()

    dof_conns = get_default(dof_conns, {})

    for key, ii, info in iter_dict_of_lists(conn_info, return_keys=True):
        ## print key, ii
        ## print info

        if info.primary is not None:
            var = info.primary
            field = var.get_field()
            field.setup_extra_data(info.ps_tg, info, info.is_trace)
            field.setup_dof_conns(dof_conns, var.n_components,
                                  info.dc_type, info.get_region())

        if info.has_virtual and not info.is_trace:
            # This is needed regardless make_virtual.
            var = info.virtual
            field = var.get_field()
            field.setup_extra_data(info.v_tg, info, False)
            field.setup_dof_conns(dof_conns, var.n_components,
                                  info.dc_type,
                                  info.get_region(can_trace=False))

    ## print dof_conns
    ## pause()
    if verbose:
        output('...done in %.2f s' % (time.clock() - tt))

    return dof_conns

##
# 14.07.2006, c
class Field( Struct ):

    @staticmethod
    def from_conf(conf, regions):
        """To refactor... very hackish now."""
        space = conf.get_default_attr('space', 'H1')
        poly_space_base = conf.get_default_attr('poly_space_base', 'lagrange')

        if isinstance(conf.region, tuple):
            region_name, kind = conf.region
            region = regions[region_name]
            if kind == 'surface':
                obj = SurfaceField(conf.name, conf.dtype, conf.shape, region,
                                   space=space,
                                   poly_space_base=poly_space_base,
                                   approx_order=conf.approx_order)

            else:
                raise ValueError('unknown field kind! (%s)', kind)

        else:
            obj = Field(conf.name, conf.dtype, conf.shape, regions[conf.region],
                        space=space,
                        poly_space_base=poly_space_base,
                        approx_order=conf.approx_order)

        return obj

    def __init__(self, name, dtype, shape, region,
                 space='H1', poly_space_base='lagrange', approx_order=1):
        """Create a Field.

        Parameters
        ----------
        name : str
            Object name.
        dtype : numpy.dtype
            Field data type: float64 or complex128.
        shape : int/tuple/str
            Field shape: 1 or (1,) or 'scalar', space dimension (2, or
            (2,) or 3 or (3,)) or 'vector'. The field shape determines
            the shape of the FE base functions and can be different from
            a FieldVariable instance shape. (TODO)

        region : Region
            The region where the field is defined.
        space : str
            The function space name.
        poly_space_base : str
            The name of polynomial space base.
        approx_order : int/str
            FE approximation order, e.g. 0, 1, 2, '1B' (1 with bubble).
        """
        if isinstance(shape, str):
            try:
                shape = {'scalar' : (1,),
                         'vector' : (region.domain.shape.dim,)}[shape]
            except KeyError:
                raise ValueError('unsupported field shape! (%s)', shape)

        elif isinstance(shape, int):
            shape = (shape,)

        Struct.__init__(self,
                        name = name,
                        dtype = dtype,
                        shape = shape,
                        region = region,
                        space = space,
                        poly_space_base = poly_space_base)
        self.domain = self.region.domain

        self.clear_dof_conns()

        self.set_approx_order(approx_order)

        # To refactor below...
        self.setup_bases()
        self.create_interpolant()
        self.setup_approximations()
##         print self.aps
##         pause()
        self.setup_global_base()
        self.setup_coors()

    def set_approx_order(self, approx_order):
        """Set a uniform approximation order."""
        
        ao_msg = 'unsupported approximation order! (%s)'
        force_bubble = False

        try:
            ao = int(approx_order)
        except ValueError:
            if approx_order[-1] == 'B':
                ao = int(approx_order[:-1])
                force_bubble = True
            else:
                raise ValueError(ao_msg % approx_order)

        if ao < 0:
            raise ValueError(ao_msg % approx_order)

        self.approx_order = '%s' % approx_order
        self._ao = ao
        self.force_bubble = force_bubble

    def setup_bases(self):
        """Setup FE bases according to self.approx_order and region cell
        types. Assumes one cell type for the whole region!"""
        gel = self.domain.groups[self.region.igs[0]].gel
        dim, n_ep = gel.dim, gel.n_vertex

        if n_ep == (dim + 1): # simplex
            kind = 'P'
        else: # tensor product
            kind = 'Q'

        self.gel = gel
        self.base_name = '%d_%d_%s%s' % (dim, n_ep, kind, self.approx_order)

    def create_interpolant(self):
        self.interp = fea.Interpolant(self.base_name, self.gel)

    def setup_approximations(self):
        self.aps = fea.Approximations(self.interp, self.region)

    ##
    #
    def igs( self ):
        return self.aps.igs

    ##
    # 19.07.2006, c
    def setup_global_base( self ):

        self.aps.describe_nodes()
        self.aps.setup_nodes()

        aux = self.aps.setup_global_base()
        self.n_nod, self.remap, self.cnt_vn, self.cnt_en = aux

##         print self.n_nod, self.cnt_vn, self.cnt_en
#        pause()

    ##
    # 19.07.2006, c
    def setup_coors( self ):
        """Coordinates of field nodes."""
        self.aps.setup_coors( self.domain.mesh, self.cnt_vn )


    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type
        
        if geometry != None:
            geometry_flag = 'surface' in geometry
        else:
            geometry_flag = False     
            
        if (dct == 'surface') or (geometry_flag):
            reg = info.get_region()
            reg.select_cells_of_surface(reset=False)

            self.aps.setup_surface_data(reg)

        elif dct == 'edge':
            raise NotImplementedError('dof connectivity type %s' % dct)

        elif dct == 'point':
            self.aps.setup_point_data(self, info.region)

        elif dct not in ('volume', 'scalar'):
            raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def clear_dof_conns(self):
        self.dof_conns = {}

    def setup_dof_conns(self, dof_conns, dpn, dc_type, region):
        """Setup dof connectivities of various kinds as needed by terms."""
        dct = dc_type.type

        ##
        # Expand nodes into dofs.
        can_point = True
        for ig, ap in self.aps.iter_aps(igs=region.igs):
            region_name = region.name # True region name.
            key = (self.name, dpn, region_name, dct, ig)
            if key in dof_conns:
                self.dof_conns[key] = dof_conns[key]

                if dct == 'point':
                    can_point = False
                continue

            if dct == 'volume':
                dc = create_dof_conn(ap.econn, dpn)
                self.dof_conns[key] = dc

            elif dct == 'surface':
                sd = ap.surface_data[region_name]
                dc = create_dof_conn(sd.econn, dpn)
                self.dof_conns[key] = dc

            elif dct == 'edge':
                raise NotImplementedError('dof connectivity type %s' % dct)

            elif dct == 'point':
                if can_point:
                    # Point data only in the first group to avoid multiple
                    # assembling of nodes on group boundaries.
                    conn = ap.point_data[region_name]
                    dc = create_dof_conn(conn, dpn)
                    self.dof_conns[key] = dc
                    can_point = False

            else:
                raise ValueError('unknown dof connectivity type! (%s)' % dct)

        dof_conns.update(self.dof_conns)

    ##
    # c: 02.01.2008, r: 02.01.2008
    def get_extra_nodes_as_simplices( self, iextra = None ):
        dim = self.domain.mesh.dim
        if iextra is None:
            noft = self.aps.node_offset_table
            iextra = nm.arange( noft[1,0], noft[-1,-1], dtype = nm.int32 )
        extra = make_point_cells( iextra, dim )
        return {2 : '2_3', 3 : '3_4'}[dim], -nm.ones_like( iextra ), extra

    def create_mesh(self, extra_nodes=True):
        """
        Create a mesh from the field region, optionally including the field
        extra nodes.
        """
        mesh = self.domain.mesh

        if self.approx_order != '0':
            conns, mat_ids, descs = [], [], []
            for ig, ap in self.aps.iter_aps():
                region = ap.region
                group = region.domain.groups[ig]
                if extra_nodes:
                    conn = ap.econn
                else:
                    offset = group.shape.n_ep
                    conn = ap.econn[:,:offset]
                conns.append(conn)
                mat_ids.append(mesh.mat_ids[ig])
                descs.append(mesh.descs[ig])

            mesh = Mesh.from_data(self.name,
                                  self.aps.coors, None, conns,
                                  mat_ids, descs)

        return mesh

    ##
    # c: 19.07.2006, r: 27.02.2008
    def write_mesh( self, name_template, field_name = None ):
        """Extra nodes are written as zero-size simplices (= points)."""
        if field_name is None:
            field_name = self.name

        tmp = self.create_mesh(extra_nodes=False)

        aux = self.get_extra_nodes_as_simplices()
        tmp.descs.append( aux[0] )
        tmp.mat_ids.append( aux[1] )
        tmp.conns.append( aux[2] )

##         print tmp
##         pause()
        tmp.write( io = 'auto' )

    ##
    # c: 20.07.2006, r: 15.01.2008
    def get_node_descs( self, region ):
        nds = {}
        for ig, ap in self.aps.iter_aps():
            if ig in region.igs:
                nds[ig] = self.aps.node_desc

        return nds

    ##
    # Modify me for bubble-only approximations to not generate vertex nodes.
    # 12.10.2005, c
    # 26.10.2005
    # 26.05.2006
    # 05.06.2006
    # 25.07.2006
    # 04.09.2006
    def interp_c_vals_to_n_vals( self, vec ):
        """len( vec ) == domain.n_el"""
        n_els = [sub.n_el for sub in self.domain.subs]
        oel = nm.cumsum( [0] + n_els )
        if sum( n_els ) != vec.shape[0]:
            print 'incomatible shape! (%d == %d)' % (sum( n_els ), vec.shape[0])
            raise ValueError

        ##
        # Mesh vertex values. 
        n_vertex = self.domain.n_nod
        nod_vol = nm.zeros( (n_vertex,), nm.float64 )
        dim = vec.shape[1]
        nod_vol_val = nm.zeros( (n_vertex, dim ), nm.float64 )
        for ii, ap in enumerate( self.aps ):
            sub = ap.sub
            ig = sub.iseq
            vg = self.vgs[ii]
            volume = nm.squeeze( vg.variable( 2 ) );

            for ii in range( sub.conn.shape[1] ):
                cc = sub.conn[:,ii]
                nod_vol[cc] += volume
                val = volume[:,nm.newaxis] * vec[oel[ig]:oel[ig+1],:]
                ind2, ind1 = nm.meshgrid( nm.arange( dim ), cc )
                nod_vol_val[ind1,ind2] += val
        
        nod_vol_val = nod_vol_val / nod_vol[:,nm.newaxis]

        ##
        # Field nodes values.
        enod_vol_val = self.interp_v_vals_to_n_vals( nod_vol_val )

        return enod_vol_val
    
    ##
    # 05.06.2006, c
    # 25.07.2006
    # 31.08.2006
    def interp_v_vals_to_n_vals( self, vec ):
        dim = vec.shape[1]
        enod_vol_val = nm.zeros( (self.n_nod, dim), nm.float64 )
        for ig, ap in self.aps.iter_aps():
            group = self.domain.groups[ig]
            offset = group.shape.n_ep
            conn = ap.econn[:,:offset]

            noff = ap.node_offsets.ravel()
            if noff[1] == noff[-1]:
                # Vertex values only...
                ii = nm.unique1d(conn) # Probably wrong?!
                enod_vol_val[ii] = vec[ii]
                continue

            econn = ap.econn

            ginterp = ap.interp.gel.interp
            coors = ap.interp.poly_spaces['v'].node_coors

            bf = ginterp.poly_spaces['v'].eval_base(coors)
            bf = bf[:,0,:].copy()
            
            fea.mu.interp_vertex_data(enod_vol_val, econn, vec, group.conn,
                                      bf, 0)

        return enod_vol_val

    ##
    # 08.08.2006, c
    # 13.02.2007
    def get_coor( self, nods = None, igs = None ):
        """Will igs be ever needed?"""
        if nods is None:
            return self.aps.coors
        else:
            return self.aps.coors[nods]

class SurfaceField(Field):
    """
    A field defined on a surface region.
    """

    def __init__(self, name, dtype, shape, region,
                 space='H1', poly_space_base='lagrange', approx_order=1):
        region.setup_face_indices()

        Field.__init__(self, name, dtype, shape, region,
                       space=space, poly_space_base=poly_space_base,
                       approx_order=approx_order)

    def setup_bases(self):
        """Setup FE bases according to self.approx_order and region cell
        types. Assumes one cell type for the whole region!"""
        gel = self.domain.groups[self.region.igs[0]].gel.surface_facet
        if gel is None:
            raise ValueError('element group has no surface!')

        dim, n_ep = gel.dim, gel.n_vertex

        if n_ep == (dim + 1): # simplex
            kind = 'P'
        else: # tensor product
            kind = 'Q'

        self.gel = gel
        self.base_name = '%d_%d_%s%s' % (dim, n_ep, kind, self.approx_order)

    def setup_approximations(self):
        self.aps = fea.Approximations(self.interp, self.region, is_surface=True)

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type

        if dct != 'surface':
            msg = "dof connectivity type must be 'surface'! (%s)" % dct
            raise ValueError(msg)

        reg = info.get_region()
        reg.select_cells_of_surface(reset=False)

        self.aps.setup_surface_data(reg)

    def setup_dof_conns(self, dof_conns, dpn, dc_type, region):
        """Setup dof connectivities of various kinds as needed by terms."""
        dct = dc_type.type

        if dct != 'surface':
            msg = "dof connectivity type must be 'surface'! (%s)" % dct
            raise ValueError(msg)

        ##
        # Expand nodes into dofs.
        for ig, ap in self.aps.iter_aps(igs=region.igs):
            region_name = region.name # True region name.
            key = (self.name, dpn, region_name, dct, ig)
            if key in dof_conns:
                self.dof_conns[key] = dof_conns[key]
                continue

            sd = ap.surface_data[region_name]
            dc = create_dof_conn(sd.leconn, dpn)
            self.dof_conns[key] = dc

        dof_conns.update(self.dof_conns)
