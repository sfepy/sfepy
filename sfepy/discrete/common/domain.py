
import numpy as nm

from sfepy.base.base import output, assert_, OneTypeList, Struct
from sfepy.base.timing import Timer
from sfepy.discrete.common.region import (Region, get_dependency_graph,
                                          sort_by_dependency, get_parents)
from sfepy.discrete.parse_regions import create_bnf, visit_stack, ParseException


def region_leaf(domain, regions, rdef, functions, tdim):
    """
    Create/setup a region instance according to rdef.
    """
    n_coor = domain.shape.n_nod
    dim = domain.shape.dim
    cmesh = domain.cmesh_tdim[tdim]

    def _region_leaf(level, op):
        token, details = op['token'], op['orig']

        if token != 'KW_Region':
            parse_def = token + '<' + ' '.join(details) + '>'
            if token != 'E_COG':
                region = Region('leaf', rdef, domain, parse_def=parse_def,
                                tdim=tdim)

        if token == 'KW_Region':
            details = details[1][2:]
            aux = regions.find(details)
            if not aux:
                raise ValueError('region %s does not exist' % details)
            else:
                if rdef[:4] == 'copy':
                    region = aux.copy()
                else:
                    region = aux

        elif token == 'KW_All':
            region.vertices = nm.arange(n_coor, dtype=nm.uint32)

        elif token == 'E_VIR':
            where = details[2]

            if where[0] == '[':
                vertices = nm.array(eval(where), dtype=nm.uint32)
                assert_(nm.amin(vertices) >= 0)
                assert_(nm.amax(vertices) < n_coor)
            else:
                coors = cmesh.coors
                y = z = None
                x = coors[:,0]

                if dim > 1:
                    y = coors[:,1]

                if dim > 2:
                    z = coors[:,2]

                coor_dict = {'x' : x, 'y' : y, 'z': z}

                vertices = nm.where(eval(where, {}, coor_dict))[0]

            region.vertices = vertices

        elif token == 'E_VOS':
            facets = cmesh.get_surface_facets()

            region.set_kind('facet')
            region.facets = facets

        elif token == 'E_VBF':
            where = details[2]

            coors = cmesh.coors

            fun = functions[where]
            vertices = fun(coors, domain=domain)

            region.vertices = vertices

        elif token == 'E_CBF':
            where = details[2]

            coors = domain.get_centroids(dim)

            fun = functions[where]
            cells = fun(coors, domain=domain)

            region.cells = cells

        elif token == 'E_COG':
            group = int(details[3])
            td = 0
            for k in range(4):
                if domain.cmesh_tdim[k] is not None:
                    cg = domain.cmesh_tdim[k].cell_groups
                    if nm.any(cg == group):
                        td = k
                        break

            if td == 0:
                raise ValueError('cell group %s does not exist' % group)

            region = Region('leaf', rdef, domain, parse_def=parse_def, tdim=td)
            region.cells = \
                nm.where(domain.cmesh_tdim[td].cell_groups == group)[0]

        elif token == 'E_COSET':
            raise NotImplementedError('element sets not implemented!')

        elif token == 'E_VOG':
            group = int(details[3])

            region.vertices = nm.where(cmesh.vertex_groups == group)[0]

        elif token == 'E_VOSET':
            try:
                vertices = domain.vertex_set_bcs[details[3]]

            except KeyError:
                msg = 'undefined vertex set! (%s)' % details[3]
                raise ValueError(msg)

            region.vertices = vertices

        elif token == 'E_OVIR':
            aux = regions[details[3][2:]]
            region.vertices = aux.vertices[0:1]

        elif token == 'E_VI':
            region.vertices = nm.array([int(ii) for ii in details[1:]],
                                       dtype=nm.uint32)

        elif token == 'E_CI':
            region.cells = nm.array([int(ii) for ii in details[1:]],
                                    dtype=nm.uint32)

        else:
            output('token "%s" unkown - check regions!' % token)
            raise NotImplementedError
        return region

    return _region_leaf

def region_op(level, op_code, item1, item2):
    token = op_code['token']
    op = {'S' : '-', 'A' : '+', 'I' : '*'}[token[3]]

    if token[-1] == 'V':
        return item1.eval_op_vertices(item2, op)

    elif token[-1] == 'E':
        return item1.eval_op_edges(item2, op)

    elif token[-1] == 'F':
        return item1.eval_op_faces(item2, op)

    elif token[-1] == 'S':
        return item1.eval_op_facets(item2, op)

    elif token[-1] == 'C':
        return item1.eval_op_cells(item2, op)

    else:
        raise ValueError('unknown region operator token! (%s)' % token)

class Domain(Struct):

    def __init__(self, name, mesh=None, nurbs=None, bmesh=None, regions=None,
                 verbose=False):
        Struct.__init__(self, name=name, mesh=mesh, nurbs=nurbs, bmesh=bmesh,
                        regions=regions, verbose=verbose)

    def get_centroids(self, dim):
        """
        Return the coordinates of centroids of mesh entities with dimension
        `dim`.
        """
        return self.cmesh.get_centroids(dim)

    def has_faces(self):
        return self.shape.tdim == 3

    def reset_regions(self):
        """
        Reset the list of regions associated with the domain.
        """
        self.regions = OneTypeList(Region)
        self._region_stack = []
        self._bnf = create_bnf(self._region_stack)

    def create_extra_tdim_region(self, region, functions, tdim):
        from sfepy.discrete.fem.geometry_element import (GeometryElement,
            create_geometry_elements)
        from sfepy.discrete import PolySpace
        """
        Create a new region which has its own cmesh with
        topological dimension tdim.
        """
        mesh = self.mesh
        if mesh.cmesh_tdim[tdim] is not None:
            raise ValueError(f'cmesh of dimension {tdim} already exists!')

        aux = mesh.from_region(region, mesh, tdim=tdim)
        cmesh = aux.cmesh
        new_mat_id = nm.max([nm.max(k.cell_groups) for k in mesh.cmesh_tdim
                             if k is not None]) + 1
        cmesh.cell_groups[:] = new_mat_id
        mesh.cmesh_tdim[tdim] = cmesh
        mesh.descs += aux.descs
        mesh.dims += aux.dims
        mesh.n_el += aux.n_el
        cmesh.set_local_entities(create_geometry_elements())
        cmesh.setup_entities()

        gel = GeometryElement(aux.descs[0])
        if gel.dim > 0:
            gel.create_surface_facet()

        new_gel_entry = {aux.descs[0]: gel}
        self.geom_els.update(new_gel_entry)

        self.fix_element_orientation(geom_els=new_gel_entry, force_check=True)

        key = gel.get_interpolation_name()

        gel.poly_space = PolySpace.any_from_args(key, gel, 1)
        gel = gel.surface_facet
        if gel is not None:
            key = gel.get_interpolation_name()
            gel.poly_space = PolySpace.any_from_args(key, gel, 1)

        select = f'cells of group {new_mat_id}'
        self._bnf.parseString(select)
        region = visit_stack(self._region_stack, region_op,
                             region_leaf(self, self.regions, select,
                                         functions, tdim))
        region.field_dim = tdim

        return region

    def create_region(self, name, select, kind='cell', parent=None,
                      check_parents=True, extra_options=None, functions=None,
                      add_to_regions=True, allow_empty=False):
        """
        Region factory constructor. Append the new region to
        self.regions list.
        """
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
            print('parsing failed:', select)
            raise

        eopts = extra_options

        tdim = self.shape.tdim if parent is None else self.regions[parent].tdim
        if eopts is not None:
            tdim = eopts.get('cell_tdim', tdim)

        region = visit_stack(stack, region_op,
                             region_leaf(self, self.regions, select,
                                         functions, tdim))

        finalize = True
        if eopts is not None:
            if 'mesh_dim' in eopts:
                region = self.create_extra_tdim_region(region, functions,
                                                       eopts['mesh_dim'])
            if not(eopts.get('finalize', True)):
                finalize = False

            if 'vertices_from' in eopts:
                vreg = eopts['vertices_from']
                region.entities[0] = self.regions[vreg].vertices.copy()
                finalize = False

        region.name = name
        region.definition = select
        region.set_kind(kind)
        if finalize:
            region.finalize(allow_empty=allow_empty)
        region.parent = parent
        region.extra_options = extra_options
        region.update_shape()

        if add_to_regions:
            self.regions.append(region)

        return region

    def create_regions(self, region_defs, functions=None, allow_empty=False):
        output('creating regions...')
        timer = Timer(start=True)

        self.reset_regions()

        ##
        # Sort region definitions by dependencies.
        graph, name_to_sort_name = get_dependency_graph(region_defs)
        sorted_regions = sort_by_dependency(graph)

        ##
        # Define regions.
        for name in sorted_regions:
            sort_name = name_to_sort_name[name]
            rdef = region_defs[sort_name]

            region = self.create_region(name, rdef.select,
                                        kind=rdef.get('kind', 'cell'),
                                        parent=rdef.get('parent', None),
                                        check_parents=False,
                                        extra_options=rdef.get('extra_options',
                                                               None),
                                        functions=functions,
                                        allow_empty=allow_empty)
            output(' ', region.name)

        output('...done in %.2f s' % timer.stop())

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
            dim = region.tdim
            is_surf = not region.can[dim] and region.can[dim - 1]
            aux = self.mesh.from_region(region, self.mesh, is_surface=is_surf)
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
            if region.vertices.shape[0] == n_nod:
                names.remove(region.name)
                names = [region.name] + names
                break

        out = {}
        for name in names:
            region = self.regions[name]
            output(region.name)

            aux.cmesh.vertex_groups[region.vertices] = n_ig
            n_ig += 1

            mask = nm.zeros((n_nod, 1), dtype=nm.float64)
            mask[region.vertices] = 1.0
            out[name] = Struct(name='region', mode='vertex', data=mask,
                               var_name=name, dofs=None)

            if region.has_cells():
                ii = region.get_cells()
                aux.cmesh.cell_groups[ii] = c_ig
                c_ig += 1

        aux.write(filename, io='auto', out=out)
        output('...done')
