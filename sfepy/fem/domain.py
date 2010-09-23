from sfepy.base.base import *
from sfepy.base.la import permutations
from geometry_element import GeometryElement
from region import Region, get_dependency_graph, sort_by_dependency, get_parents
from sfepy.fem.parseReg \
     import create_bnf, visit_stack, print_stack, ParseException
import fea
import extmods.meshutils as mu

def _build_orientation_map(n_fp):
    """
    The keys are binary masks of the lexicographical ordering of facet
    vertices. A bit i set to one means `v[i] < v[i+1]`.

    The values are `[original_order, permutation, orientation]`, where
    `permutation` can be used to sort facet vertices lexicographically,
    and `orientation` is the order of the first vertex + 1 times the
    sign. Hence `permuted_facet = facet[permutation]`.
    """
    indices = range(n_fp)

    cmps = [(i1, i2) for i2 in indices for i1 in indices[:i2]]
    powers = [2**ii for ii in range(len(cmps))]

    ori_map = {}
    for indx in permutations(indices):
        key = 0
        sign = 1
        for ip, power in enumerate(powers):
            i1, i2 = cmps[ip]
            less = (indx[i1] < indx[i2])
            key += power * less
            if not less:
                sign *= -1

        isort = nm.argsort(indx)
        ori_map[key] = [indx, isort, sign * (indx[0] + 1)]

    return ori_map, cmps, powers

def _get_signed_orientation(ori, ori_map):
    """
    Transform orientation according to `ori_map`, i.e. from bit mask
    encoding to signed values.
    """
    i_from = nm.array(ori_map.keys())
    i_to = [ii[-1] for ii in ori_map.values()]

    i_map = nm.zeros((i_from.max() + 1,), dtype=nm.int8)
    i_map[i_from] = i_to

    signed_ori = i_map[ori]
    assert_((signed_ori != 0).all())

    return signed_ori

def _permute_facets(facets, ori, ori_map):
    """
    Return a copy of `facets` array with vertices sorted lexicographically.
    """
    assert_((nm.setmember1d(nm.unique(ori), ori_map.keys())).all())

    permuted_facets = facets.copy()

    for key, ori_map in ori_map.iteritems():
        perm = ori_map[1]
        ip = nm.where(ori == key)[0]
        for ic0, ic1 in enumerate(perm):
            permuted_facets[ip,ic0] = facets[ip,ic1]

    return permuted_facets

def _orient_facets(ori, facets, cmps, powers):
    for ip, power in enumerate(powers):
        i1, i2 = cmps[ip]
        iw = nm.where(facets[:,i1] < facets[:,i2])[0]
        ori[iw] += power

class Facets(Struct):

    @staticmethod
    def from_domain(domain, kind):
        groups = domain.groups

        if kind == 'edges':
            n_obj = [group.shape.n_edge_total for group in groups.itervalues()]
            nn = 2
        else:
            n_obj = [group.shape.n_face_total for group in groups.itervalues()]
            nn = 4

        n_all_obj = sum(n_obj)
        indices = nm.zeros((n_all_obj, 3), nm.int32)
        all_facets = nm.empty((n_all_obj, nn), nm.int32)
        all_facets.fill(-1)

        single_facets = {}

        ii = 0
        for ig, group in groups.iteritems():
            conn, gel = group.conn, group.gel
            n_el = group.shape.n_el
            n_all_item = n_obj[ig]

            if kind == 'edges':
                n_item, items = gel.n_edge, gel.edges

            else:
                n_item, items = gel.n_face, gel.faces

            n_fp = items.shape[1]
            single_facets[ig] = items

            io = slice(ii, ii + n_all_item)

            indices[io,0] = ig

            ie = nm.arange(n_el, dtype=nm.int32)
            ie = nm.repeat(ie, n_item)
            indices[io,1] = ie

            iobj = nm.arange(n_item, dtype=nm.int32)
            indices[io,2] = nm.tile(iobj, n_el)

            facets = conn[:, items]
            facets = facets.reshape((n_all_item, n_fp))

            all_facets[io,:n_fp] = facets

            ii += n_all_item

        if (ii != sum( n_obj )):
            msg = 'neighbour_list size mismatch! (%d == %d = sum( %s ))'\
                  % (ii, sum( n_obj ), n_obj)
            raise ValueError( msg )

        obj = Facets('facets', kind, domain, single_facets,
                     n_obj, indices, all_facets)

        return obj

    def __init__(self, name, kind, domain, single_facets,
                 n_obj, indices, facets):
        Struct.__init__(self, name=name, kind=kind, domain=domain,
                        single_facets=single_facets,
                        n_obj=n_obj, indices=indices, facets=facets)
        self.n_all_obj, self.n_col = facets.shape
        self.n_gr = len(self.n_obj)

        self.indx = {}
        ii = 0
        for ig, nn in enumerate(self.n_obj):
            self.indx[ig] = slice(ii, ii+nn)
            ii += nn

        self.n_fps_vec = nm.empty(self.n_all_obj, dtype=nm.int32)
        self.n_fps = {}
        for ig, facet in self.single_facets.iteritems():
            self.n_fps_vec[self.indx[ig]] = facet.shape[1]
            self.n_fps[ig] = facet.shape[1]

    def sort_and_orient(self):
        all_permuted_facets = nm.empty((self.n_all_obj + 2, self.n_col),
                                     dtype=nm.int32)
        all_permuted_facets.fill(-1)

        sentinel = self.domain.shape.n_nod

        aux = nm.repeat(nm.array([sentinel], nm.int32), self.n_col)
        all_permuted_facets[-2] = aux
        all_permuted_facets[-1] = aux + 1
        oris = {}
        ori_maps = {}

        for ig in range(self.n_gr):
            io = self.indx[ig]
            facets = self.facets[io]

            n_fp = self.n_fps[ig]
            ori_map, cmps, powers = _build_orientation_map(n_fp)

            # Determine orientation.
            ori = nm.zeros((facets.shape[0],), dtype=nm.int8)
            _orient_facets(ori, facets, cmps, powers)

            # Permute each facet to have indices in ascending order, so
            # that lexicographic sorting works.
            permuted_facets = _permute_facets(facets, ori, ori_map)
            all_permuted_facets[io] = permuted_facets

            ori = _get_signed_orientation(ori, ori_map)

            oris[ig] = ori
            ori_maps[ig] = ori_map

        self.permuted_facets = all_permuted_facets
        self.oris = oris
        self.ori_maps = ori_maps

    def setup_unique(self):
        """
        `sorted_facets` == `permuted_facets[perm]`
        `permuted_facets` == `sorted_facets[perm_i]`
        `uid` : unique id in order of `sorted_facets`
        `uid_i` : unique id in order of `permuted_facets` or `facets`
        """
        ii = nm.arange(self.permuted_facets.shape[0], dtype=nm.int32)
        aux = nm.concatenate((self.permuted_facets, ii[:,None]), 1).copy()
        sort_cols = nm.arange(self.permuted_facets.shape[1],
                              dtype=nm.int32)

        mu.sort_rows(aux, sort_cols)

        self.perm = perm = aux[:,-1].copy()
        self.sorted_facets = aux[:,:-1].copy()
        aux = nm.arange(perm.shape[0], dtype=nm.int32)
        self.perm_i = nm.zeros_like(self.perm)
        self.perm_i[perm] = aux

        ic = nm.where(nm.abs(nm.diff(self.sorted_facets, axis=0)).sum(1), 0, 1)
        ic = ic.astype(nm.int32)

        self.n_unique = len(ic) - nm.sum(ic) - 1
        self.unique_list = nm.where(ic[:-1] == 0)[0].astype(nm.int32)
        assert_(len(self.unique_list) == self.n_unique)

        ii = nm.cumsum( ic[:-1] == 0, dtype = nm.int32 )
        self.uid = ii.copy()
        self.uid[0], self.uid[1:] = 0, ii[:-1]
        self.uid_i = self.uid[self.perm_i[:-2]]

    def setup_neighbours(self):
        """
        For each unique facet:
           - indices of facets - sparse matrix (n_unique x n_all_obj)
             mtx[i, j] == 1 if facet[j] has uid[i]
           - number of elements it is in
        """
        ones = nm.ones((self.n_all_obj,), dtype=nm.bool)

        self.mtx = sp.coo_matrix((ones, (self.uid, self.perm[:-2])))

        self.n_in_el = self.mtx * ones.astype(nm.int32)

    def mark_surface_facets(self):
        """
        flag: 0 .. inner, 2 .. edge, 3 .. triangle, 4 .. quadrangle
        """
        nn = self.n_in_el[self.uid_i]

        flag = self.n_fps_vec * (nn == 1)

        return flag

##
# 14.06.2006, c
# 15.06.2006
# 19.02.2007
# 02.03.2007
# 02.05.2007
# 30.05.2007
# 05.06.2007
def region_leaf(domain, regions, rdef, functions):
    """
    Create/setup a region instance according to rdef.
    """
    def _region_leaf( level, op ):

        token, details = op['token'], op['orig']

        if token != 'KW_Region':
            parse_def = token + '<' + ' '.join( details ) + '>'
            region = Region('leaf', rdef, domain, parse_def=parse_def)

        if token == 'KW_Region':
            details = details[1][2:]
            aux = regions.find( details )
            if not aux:
                raise ValueError, 'region %s does not exist' % details
            else:
                if rdef[:4] == 'copy':
                    region = aux.copy()
                else:
                    region = aux

        elif token == 'KW_All':
            region.set_vertices( nm.arange( domain.mesh.n_nod,
                                            dtype = nm.int32 ) )
        elif token == 'E_NIR':
            where = details[2]

            if where[0] == '[':
                out = nm.array( eval( where ), dtype = nm.int32 )
                assert_( nm.amin( out ) >= 0 )
                assert_( nm.amax( out ) < domain.mesh.n_nod )
            else:
                coors = domain.get_mesh_coors()
                x = coors[:,0]
                y = coors[:,1]
                if domain.mesh.dim == 3:
                    z = coors[:,2]
                else:
                    z = None
                coor_dict = {'x' : x, 'y' : y, 'z': z}

                out = nm.where( eval( where, {}, coor_dict ) )[0]
            region.set_vertices( out )

        elif token == 'E_NOS':

            if domain.fa: # 3D.
                fa = domain.fa
            else:
                fa = domain.ed

            flag = fa.mark_surface_facets()
            ii = nm.where( flag > 0 )[0]
            aux = nm.unique1d(fa.facets[ii])
            if aux[0] == -1: # Triangular faces have -1 as 4. point.
                aux = aux[1:]

            region.can_cells = False
            region.set_vertices( aux )

        elif token == 'E_NBF':
            where = details[2]

            coors = domain.get_mesh_coors()

            fun = functions[where]
            out = fun(coors, domain=domain)

            region.set_vertices( out )

        elif token == 'E_EBF':
            where = details[2]

            coors = domain.get_mesh_coors()

            fun = functions[where]
            out = fun(coors, domain=domain)

            region.set_cells( out )

        elif token == 'E_EOG':

            group = int( details[3] )

            ig = domain.mat_ids_to_i_gs[group]
            group = domain.groups[ig]
            region.set_from_group( ig, group.vertices, group.shape.n_el )

        elif token == 'E_NOG':

            group = int( details[3] )

            group_nodes = nm.where( domain.mesh.ngroups == group )[0]
            region.set_vertices( group_nodes )

        elif token == 'E_ONIR':
            aux = regions[details[3][2:]]
            region.set_vertices( aux.all_vertices[0:1] )

        elif token == 'E_NI':
            region.set_vertices(nm.array([int(ii) for ii in details[1:]],
                                         dtype=nm.int32))

        elif token == 'E_EI1':
            region.set_cells({0 : nm.array([int(ii) for ii in details[1:]],
                                           dtype=nm.int32)})

        elif token == 'E_EI2':
            num = len(details[1:]) / 2

            cells = {}
            for ii in range(num):
                ig, iel = int(details[1+2*ii]), int(details[2+2*ii])
                cells.setdefault(ig, []).append(iel)

            region.set_cells(cells)

        else:
            output( 'token "%s" unkown - check regions!' % token )
            raise NotImplementedError
        return region

    return _region_leaf

def region_op(level, op, item1, item2):
    token = op['token']
    if token == 'OA_SubN':
        return item1.sub_n( item2 )
    elif token == 'OA_SubE':
        return item1.sub_e( item2 )
    elif token == 'OA_AddN':
        return item1.add_n( item2 )
    elif token == 'OA_AddE':
        return item1.add_e( item2 )
    elif token == 'OA_IntersectN':
        return item1.intersect_n( item2 )
    elif token == 'OA_IntersectE':
        return item1.intersect_e( item2 )
    else:
        raise NotImplementedError, token

##
# 17.07.2006, c
class Domain( Struct ):
    """
    Domain is divided into groups, whose purpose is to have homogeneous
    data shapes."""

    def __init__(self, name, mesh):
        """Create a Domain.

        Parameters
        ----------
        name : str
            Object name.
        mesh : Mesh
            A mesh defining the domain.
        """
        geom_els = {}
        for ig, desc in enumerate(mesh.descs):
            gel = GeometryElement(desc)
            # Create geometry elements of dimension - 1.
            gel.create_surface_facet()

            geom_els[desc] = gel

        interps = {}
        for gel in geom_els.itervalues():
            key = gel.get_interpolation_name()

            gel.interp = interps.setdefault(key,
                                            fea.Interpolant(key, gel))
            gel = gel.surface_facet
            if gel is not None:
                key = gel.get_interpolation_name()
                gel.interp = interps.setdefault(key,
                                                fea.Interpolant(key, gel))

        Struct.__init__(self,
                        name = name,
                        mesh = mesh,
                        geom_els = geom_els,
                        geom_interps = interps)
        self.mat_ids_to_i_gs = {}
        for ig, mat_id in enumerate( mesh.mat_ids ):
            self.mat_ids_to_i_gs[mat_id[0]] = ig

        self.setup_groups()
        self.fix_element_orientation()
        self.setup_neighbour_lists()
        self.reset_regions()

    def setup_groups( self ):

        n_nod, dim = self.mesh.coors.shape
        self.shape = Struct( n_gr = len( self.mesh.conns ), n_el = 0,
                             n_nod = n_nod, dim = dim )
        self.groups = {}
        for ii in range( self.shape.n_gr ):
            gel = self.geom_els[self.mesh.descs[ii]] # Shortcut.
            conn = self.mesh.conns[ii]
            vertices = nm.unique1d( conn )

            n_vertex = vertices.shape[0]
            n_el, n_ep = conn.shape
            n_edge = gel.n_edge           
            n_edge_total = n_edge * n_el

            if gel.dim == 3:
                n_face = gel.n_face
                n_face_total = n_face * n_el
            else:
                n_face = n_face_total = 0

            shape = Struct( n_vertex = n_vertex, n_el = n_el, n_ep = n_ep,
                            n_edge = n_edge, n_edge_total = n_edge_total,
                            n_face = n_face, n_face_total = n_face_total )
            self.groups[ii] = Struct( ig = ii,
                                      vertices = vertices,
                                      conn = conn,
                                      gel = gel,
                                      shape = shape )
            self.shape.n_el += n_el

    ##
    # c: 22.11.2007, r: 25.03.2008
    def iter_groups( self, igs = None ):
        if igs is None:
            for ig in xrange( self.shape.n_gr ): # sorted by ig.
                yield self.groups[ig]
        else:
            for ig in igs:
                yield ig, self.groups[ig]

    ##
    # c: 25.03.2008, r: 25.03.2008
    def get_cell_offsets( self ):
        offs = {}
        off = 0
        for group in self.iter_groups():
            ig = group.ig
            offs[ig] = off
            off += group.shape.n_el
        return offs

    def get_mesh_coors(self):
        """
        Return the coordinates of the underlying mesh vertices.
        """
        return self.mesh.coors

    def get_mesh_bounding_box(self):
        """
        Return the bounding box of the underlying mesh. 

        Returns
        -------
        bbox : ndarray (2, dim)
            The bounding box with min. values in the first row and max. values
            in the second row.
        """
        return self.mesh.get_bounding_box()

    def get_diameter(self):
        """
        Return the diameter of the domain.

        Notes
        -----
        The diameter corresponds to the Friedrichs constant.
        """
        bbox = self.get_mesh_bounding_box()
        return (bbox[1,:] - bbox[0,:]).max()

    def get_conns(self):
        """
        Return the element connectivity groups of the underlying mesh.
        """
        return self.mesh.conns
    
    ##
    # 08.06.2006, c
    # 17.07.2006
    def fix_element_orientation( self ):

        coors = self.mesh.coors
        for ii, group in self.groups.iteritems():

            ori, conn = group.gel.orientation, group.conn

            itry = 0
            while itry < 2:
                flag = -nm.ones( conn.shape[0], dtype = nm.int32 )

                # Changes orientation if it is wrong according to swap*!
                # Changes are indicated by positive flag.
                mu.orient_elements( flag, conn, coors,
                                    ori.roots, ori.vecs,
                                    ori.swap_from, ori.swap_to )
    #            print flag
                if nm.alltrue( flag == 0 ):
                    if itry > 0: output( '...corrected' )
                    itry = -1
                    break

                output('warning: bad element orientation, trying to correct...')
                itry += 1

            if itry == 2 and flag[0] != -1:
                raise RuntimeError, "elements cannot be oriented! (%d, %s)"\
                      % (ii, self.mesh.descs[ii] )
            elif flag[0] == -1:
                output( 'warning: element orienation not checked' )

    ##
    # 19.07.2006, c
    def get_orientation( self, ii, mode = 'edges' ):
        group = self.groups[ii]

        if mode == 'edges':
            oo = nm.reshape(self.ed.oris[ii],
                            (group.shape.n_el, group.gel.n_edge))
            ori = ((1 - oo) / 2).astype(nm.int32)

        elif mode == 'faces':
            output( 'orient faces' )
            raise NotImplementedError
        else:
            output( mode )
            raise ValueError
        
        return ori

    def has_faces( self ):
        return sum( [group.shape.n_face
                     for group in self.iter_groups()] ) > 0
        

    def setup_neighbour_lists(self, create_edge_list=True,
                              create_face_list=True):
        kinds = ['edges', 'faces']

        is_face = self.has_faces()
        flags = [create_edge_list, create_face_list and is_face]

        for ii, kind in enumerate(kinds):
            if flags[ii]:
                output('setting up domain %s...' % kind)

                tt = time.clock()
                obj = Facets.from_domain(self, kind)
                obj.sort_and_orient()
                obj.setup_unique()
                obj.setup_neighbours()

                if kind == 'edges':
                    self.ed = obj

                else:
                    self.fa = obj

                output('...done in %.2f s' % (time.clock() - tt))

        if not is_face:
            self.fa = None

    def get_facets(self, force_faces=False):
        """
        Return edge and face descriptions.
        """
        if force_faces and not self.fa:
            return self.ed, self.ed

        else:
            return self.ed, self.fa

    def reset_regions(self):
        """Reset the list of regions associated with the domain."""
        self.regions = OneTypeList(Region)
        self._region_stack = []
        self._bnf = create_bnf(self._region_stack)

    def create_region(self, name, select, flags=None, check_parents=True,
                      functions=None, add_to_regions=True):
        """Region factory constructor. Append the new region to
        self.regions list."""
        if flags is None:
            flags = {}

        if check_parents:
            parents = get_parents(select)
            for p in parents:
                if p not in [region.name for region in self.regions]:
                    msg = 'parent region %s of %s not found!' % (p, name)
                    raise ValueError(msg)

        stack = self._region_stack
        try:
            self._bnf.parseString(select)
        except ParseException:
            print 'parsing failed:', select
            raise

        region = visit_stack(stack, region_op,
                             region_leaf(self, self.regions, select,
                                         functions))
        region.name = name

        forbid = flags.get('forbid', None)
        if forbid:
            fb = re.compile('^group +\d+(\s+\d+)*$').match(forbid)
            if fb:
                groups = forbid[5:].strip().split()
                forbid = [int( ii ) for ii in groups]
            else:
                raise ValueError('bad forbid! (%s)' % forbid)
            forbidden_igs = [self.mat_ids_to_i_gs[mat_id] for mat_id in forbid]
            region.delete_groups(forbidden_igs)

        region.switch_cells(flags.get('can_cells', True))
            
        region.complete_description(self.ed, self.fa)

        if add_to_regions:
            self.regions.append(region)

        return region
            
    def create_regions(self, region_defs, functions=None):
        output( 'creating regions...' )
        tt = time.clock()

        self.reset_regions()

        ##
        # Sort region definitions by dependencies.
        graph, name_to_sort_name = get_dependency_graph(region_defs)
        sorted_regions = sort_by_dependency(graph)
##         print sorted_regions
        
        ##
        # Define regions.
        for name in sorted_regions:
            sort_name = name_to_sort_name[name]
            rdef = region_defs[sort_name]

            region = self.create_region(name, rdef.select,
                                        flags=rdef,
                                        check_parents=False,
                                        functions=functions)
            output(' ', region.name)

        output( '...done in %.2f s' % (time.clock() - tt) )

        return self.regions

    ##
    # 31.07.2007, c
    def get_element_diameters( self, ig, cells, vg, mode, square = True ):
        group = self.groups[ig]
        diameters = nm.empty( (len( cells ), 1, 1, 1), dtype = nm.float64 )
        if vg is None:
            diameters.fill( 1.0 )
        else:
            vg.get_element_diameters( diameters, group.gel.edges,
                                    self.get_mesh_coors().copy(), group.conn,
                                    cells, mode )
        if square:
            out = diameters.squeeze()
        else:
            out = nm.sqrt( diameters.squeeze() )

        return out
            
    ##
    # 29.08.2007, re-c from 00.01.18, r: 26.03.2008
    def surface_faces( self ):

        if not self.fa:
            print "no faces defined!"
            raise ValueError

        fa = self.fa
        flag = fa.mark_surface_facets()

        surf_faces = []
        itri = nm.where( flag == 1 )[0]
        if itri.size:
            surf_faces.append( fa.facets[itri,:3] )
        itet = nm.where( flag == 2 )[0]
        if itet.size:
            surf_faces.append( fa.facets[itet,:4] )

        isurf = nm.where( flag >= 1 )[0]
        if isurf.size:
            lst = fa.indices[isurf]

        return lst, surf_faces
