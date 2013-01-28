import time
import re

import numpy as nm

from sfepy.base.base import output, assert_, OneTypeList, Struct
from sfepy.fem.facets import Facets
from geometry_element import GeometryElement
from region import Region, get_dependency_graph, sort_by_dependency, get_parents
from sfepy.fem.parseReg import create_bnf, visit_stack, ParseException
from sfepy.fem.refine import refine_2_3, refine_2_4, refine_3_4, refine_3_8
from sfepy.fem.fe_surface import FESurface
import fea

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
            aux = nm.unique(fa.facets[ii])
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
            group = int(details[3])

            ig = domain.mat_ids_to_i_gs[group]
            group = domain.groups[ig]
            region.set_from_group(ig, group.vertices, group.shape.n_el)

        elif token == 'E_EOSET':
            raise NotImplementedError('element sets not implemented!')

        elif token == 'E_NOG':
            group = int(details[3])
            group_nodes = nm.where(domain.mesh.ngroups == group)[0]

            region.set_vertices(group_nodes)

        elif token == 'E_NOSET':
            try:
                set_nodes = domain.mesh.nodal_bcs[details[3]]

            except KeyError:
                msg = 'undefined nodal set! (%s)' % details[3]
                raise ValueError(msg)

            region.set_vertices(set_nodes)

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

    def __init__(self, name, mesh, verbose=False):
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
        self.setup_facets(verbose=verbose)
        self.reset_regions()
        self.clear_surface_groups()

    def setup_groups( self ):

        n_nod, dim = self.mesh.coors.shape
        self.shape = Struct( n_gr = len( self.mesh.conns ), n_el = 0,
                             n_nod = n_nod, dim = dim )
        self.groups = {}
        for ii in range( self.shape.n_gr ):
            gel = self.geom_els[self.mesh.descs[ii]] # Shortcut.
            conn = self.mesh.conns[ii]
            vertices = nm.unique( conn )

            n_vertex = vertices.shape[0]
            n_el, n_ep = conn.shape
            n_edge = gel.n_edge           
            n_edge_total = n_edge * n_el

            if gel.dim == 3:
                n_face = gel.n_face
                n_face_total = n_face * n_el
            else:
                n_face = n_face_total = 0

            shape = Struct(n_vertex=n_vertex, n_el=n_el, n_ep=n_ep,
                           n_edge=n_edge, n_edge_total=n_edge_total,
                           n_face=n_face, n_face_total=n_face_total,
                           dim=self.mesh.dims[ii])

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

    def get_mesh_coors(self, actual=False):
        """
        Return the coordinates of the underlying mesh vertices.
        """
        if actual and hasattr(self.mesh, 'coors_act'):
            return self.mesh.coors_act
        else:
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

    def fix_element_orientation(self):
        """
        Ensure element nodes ordering giving positive element volume.

        The groups with elements of lower dimension than the space dimension
        are skipped.
        """
        from extmods.mesh import orient_elements

        coors = self.mesh.coors
        for ii, group in self.groups.iteritems():
            if group.shape.dim < self.shape.dim: continue

            ori, conn = group.gel.orientation, group.conn

            itry = 0
            while itry < 2:
                flag = -nm.ones(conn.shape[0], dtype=nm.int32)

                # Changes orientation if it is wrong according to swap*!
                # Changes are indicated by positive flag.
                orient_elements(flag, conn, coors,
                                ori.roots, ori.vecs,
                                ori.swap_from, ori.swap_to)

                if nm.alltrue( flag == 0 ):
                    if itry > 0: output('...corrected')
                    itry = -1
                    break

                output('warning: bad element orientation, trying to correct...')
                itry += 1

            if itry == 2 and flag[0] != -1:
                raise RuntimeError('elements cannot be oriented! (%d, %s)'
                                   % (ii, self.mesh.descs[ii]))
            elif flag[0] == -1:
                output('warning: element orienation not checked')

    def has_faces( self ):
        return sum( [group.shape.n_face
                     for group in self.iter_groups()] ) > 0

    def setup_facets(self, create_edges=True, create_faces=True,
                     verbose=False):
        """
        Setup the edges and faces (in 3D) of domain elements.
        """
        kinds = ['edges', 'faces']

        is_face = self.has_faces()
        create = [create_edges, create_faces and is_face]

        for ii, kind in enumerate(kinds):
            if create[ii]:
                if verbose:
                    output('setting up domain %s...' % kind)

                tt = time.clock()
                obj = Facets.from_domain(self, kind)
                obj.sort_and_orient()
                obj.setup_unique()
                obj.setup_neighbours()

                # 'ed' or 'fa'
                setattr(self, kind[:2], obj)

                if verbose:
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

    def save_regions(self, filename, region_names=None):
        """
        Save regions as individual meshes.

        Parameters
        ----------
        filename : str
            The output filename.
        region_names : list, optional
            If given, only the listed regions are saved.
        """
        import os

        if region_names is None:
            region_names = self.regions.get_names()

        trunk, suffix = os.path.splitext(filename)

        output('saving regions...')
        for name in region_names:
            region = self.regions[name]
            output(name)
            aux = self.mesh.from_region(region, self.mesh)
            aux.write('%s_%s%s' % (trunk, region.name, suffix),
                      io='auto')
        output('...done')

    def save_regions_as_groups(self, filename, region_names=None):
        """
        Save regions in a single mesh but mark them by using different
        element/node group numbers.

        If regions overlap, the result is undetermined, with exception of the
        whole domain region, which is marked by group id 0.

        Region masks are also saved as scalar point data for output formats
        that support this.

        Parameters
        ----------
        filename : str
            The output filename.
        region_names : list, optional
            If given, only the listed regions are saved.
        """

        output('saving regions as groups...')
        aux = self.mesh.copy()
        n_ig = c_ig = 0
        n_nod = self.shape.n_nod

        # The whole domain region should go first.
        names = (region_names if region_names is not None
                 else self.regions.get_names())
        for name in names:
            region = self.regions[name]
            if region.all_vertices.shape[0] == n_nod:
                names.remove(region.name)
                names = [region.name] + names
                break

        out = {}
        for name in names:
            region = self.regions[name]
            output(region.name)

            aux.ngroups[region.all_vertices] = n_ig
            n_ig += 1

            mask = nm.zeros((n_nod, 1), dtype=nm.float64)
            mask[region.all_vertices] = 1.0
            out[name] = Struct(name='region', mode='vertex', data=mask,
                               var_name=name, dofs=None)

            if region.has_cells():
                for ig in region.igs:
                    ii = region.get_cells(ig)
                    aux.mat_ids[ig][ii] = c_ig
                    c_ig += 1

        aux.write(filename, io='auto', out=out)
        output('...done')

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
        itri = nm.where(flag == 3)[0]
        if itri.size:
            surf_faces.append( fa.facets[itri,:3] )
        itet = nm.where(flag == 4)[0]
        if itet.size:
            surf_faces.append( fa.facets[itet,:4] )

        isurf = nm.where( flag >= 1 )[0]
        if isurf.size:
            lst = fa.indices[isurf]

        return lst, surf_faces

    def clear_surface_groups(self):
        """
        Remove surface group data.
        """
        self.surface_groups = {}

    def create_surface_group(self, region):
        """
        Create a new surface group corresponding to `region` if it does
        not exist yet.

        Notes
        -----
        Surface groups define surface facet connectivity that is needed
        for :class:`sfepy.fem.mappings.SurfaceMapping`.
        """
        for ig in region.igs:
            groups = self.surface_groups.setdefault(ig, {})
            if region.name not in groups:
                region.select_cells_of_surface(reset=False)

                group = self.groups[ig]
                gel_faces = group.gel.get_surface_entities()

                name = 'surface_group_%s_%d' % (region.name, ig)
                surface_group = FESurface(name, region, gel_faces,
                                          group.conn, ig)

                groups[region.name] = surface_group

    def refine(self):
        """
        Uniformly refine the domain mesh.

        Returns
        -------
        domain : Domain instance
            The new domain with the refined mesh.

        Notes
        -----
        Works only for meshes with single element type! Does not
        preserve node groups!
        """

        names = set()
        for group in self.groups.itervalues():
            names.add(group.gel.name)

        if len(names) != 1:
            msg = 'refine() works only for meshes with single element type!'
            raise NotImplementedError(msg)

        el_type = names.pop()
        if el_type == '2_3':
            mesh = refine_2_3(self.mesh, self.ed)

        elif el_type == '2_4':
            mesh = refine_2_4(self.mesh, self.ed)

        elif el_type == '3_4':
            mesh = refine_3_4(self.mesh, self.ed)

        elif el_type == '3_8':
            mesh = refine_3_8(self.mesh, self.ed, self.fa)

        else:
            msg = 'unsupported element type! (%s)' % el_type
            raise NotImplementedError(msg)

        domain = Domain(self.name + '_r', mesh)

        return domain
