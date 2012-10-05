"""
Notes
-----

Important attributes of continuous (order > 0) :class:`Field` and
:class:`SurfaceField` instances:

- `vertex_remap` : `econn[:, :n_vertex] = vertex_remap[conn]`
- `vertex_remap_i` : `conn = vertex_remap_i[econn[:, :n_vertex]]`

where `conn` is the mesh vertex connectivity, `econn` is the
region-local field connectivity.
"""
import time
import numpy as nm

from sfepy.base.base import output, iter_dict_of_lists, get_default, assert_
from sfepy.base.base import Struct, basestr
import fea
from sfepy.fem.mesh import Mesh
from sfepy.fem.meshio import convert_complex_output
from sfepy.fem.utils import extend_cell_data, prepare_remap, invert_remap
from sfepy.fem.fe_surface import FESurface
from sfepy.fem.integrals import Integral
from sfepy.fem.linearizer import get_eval_dofs, get_eval_coors, create_output

def parse_approx_order(approx_order):
    """
    Parse the uniform approximation order value (str or int).
    """
    ao_msg = 'unsupported approximation order! (%s)'
    force_bubble = False
    discontinuous = False

    try:
        ao = int(approx_order)
    except ValueError:
        mode = approx_order[-1].lower()
        if mode == 'b':
            ao = int(approx_order[:-1])
            force_bubble = True

        elif mode == 'd':
            ao = int(approx_order[:-1])
            discontinuous = True

        else:
            raise ValueError(ao_msg % approx_order)

    if ao < 0:
        raise ValueError(ao_msg % approx_order)

    elif ao == 0:
        discontinuous = True

    return ao, force_bubble, discontinuous

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
            iloc = n_ep * idof + inod # Hack: For DBD order.
            dc[:,iloc] = dpn * conn[:,inod] + idof

    return dc

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

        if info.primary is not None:
            var = info.primary
            field = var.get_field()
            field.setup_extra_data(info.ps_tg, info, info.is_trace)
            field.setup_dof_conns(dof_conns, var.n_components,
                                  info.dc_type, info.get_region(),
                                  info.is_trace)

        if info.has_virtual and not info.is_trace:
            # This is needed regardless make_virtual.
            var = info.virtual
            field = var.get_field()
            field.setup_extra_data(info.v_tg, info, False)
            field.setup_dof_conns(dof_conns, var.n_components,
                                  info.dc_type,
                                  info.get_region(can_trace=False))

    if verbose:
        output('...done in %.2f s' % (time.clock() - tt))

    return dof_conns

def get_eval_expression(expression, ig,
                        fields, materials, variables,
                        functions=None, mode='eval', term_mode=None,
                        extra_args=None, verbose=True, kwargs=None):
    """
    Get the function for evaluating an expression given a list of elements,
    and reference element coordinates.
    """
    from sfepy.fem.evaluate import eval_in_els_and_qp

    def _eval(iels, coors):
        val = eval_in_els_and_qp(expression, ig, iels, coors,
                                 fields, materials, variables,
                                 functions=functions, mode=mode,
                                 term_mode=term_mode,
                                 extra_args=extra_args, verbose=verbose,
                                 kwargs=kwargs)
        return val[..., 0]

    return _eval

def create_expression_output(expression, name, primary_field_name,
                             fields, materials, variables,
                             functions=None, mode='eval', term_mode=None,
                             extra_args=None, verbose=True, kwargs=None,
                             min_level=0, max_level=1, eps=1e-4):
    """
    Create output mesh and data for the expression using the adaptive
    linearizer.

    Parameters
    ----------
    expression : str
        The expression to evaluate.
    name : str
        The name of the data.
    primary_field_name : str
        The name of field that defines the element groups and polynomial
        spaces.
    fields : dict
        The dictionary of fields used in `variables`.
    materials : Materials instance
        The materials used in the expression.
    variables : Variables instance
        The variables used in the expression.
    functions : Functions instance, optional
        The user functions for materials etc.
    mode : one of 'eval', 'el_avg', 'qp'
        The evaluation mode - 'qp' requests the values in quadrature points,
        'el_avg' element averages and 'eval' means integration over
        each term region.
    term_mode : str
        The term call mode - some terms support different call modes
        and depending on the call mode different values are
        returned.
    extra_args : dict, optional
        Extra arguments to be passed to terms in the expression.
    verbose : bool
        If False, reduce verbosity.
    kwargs : dict, optional
        The variables (dictionary of (variable name) : (Variable
        instance)) to be used in the expression.
    min_level : int
        The minimum required level of mesh refinement.
    max_level : int
        The maximum level of mesh refinement.
    eps : float
        The relative tolerance parameter of mesh adaptivity.

    Returns
    -------
    out : dict
        The output dictionary.
    """
    field = fields[primary_field_name]
    vertex_coors = field.coors[:field.n_vertex_dof, :]

    coors = []
    vdofs = []
    conns = []
    mat_ids = []
    levels = []
    offset = 0
    for ig, ap in field.aps.iteritems():
        ps = ap.interp.poly_spaces['v']
        gps = ap.interp.gel.interp.poly_spaces['v']
        group = field.domain.groups[ig]
        vertex_conn = ap.econn[:, :group.shape.n_ep]

        eval_dofs = get_eval_expression(expression, ig,
                                        fields, materials, variables,
                                        functions=functions,
                                        mode=mode, extra_args=extra_args,
                                        verbose=verbose, kwargs=kwargs)
        eval_coors = get_eval_coors(vertex_coors, vertex_conn, gps)

        (level, _coors, conn,
         _vdofs, _mat_ids) = create_output(eval_dofs, eval_coors,
                                           group.shape.n_el, ps,
                                           min_level=min_level,
                                           max_level=max_level, eps=eps)

        _mat_ids[:] = field.domain.mesh.mat_ids[ig][0]

        coors.append(_coors)
        vdofs.append(_vdofs)
        conns.append(conn + offset)
        mat_ids.append(_mat_ids)
        levels.append(level)

        offset += _coors.shape[0]

    coors = nm.concatenate(coors, axis=0)
    vdofs = nm.concatenate(vdofs, axis=0)
    mesh = Mesh.from_data('linearized_mesh', coors, None, conns, mat_ids,
                          field.domain.mesh.descs)

    out = {}
    out[name] = Struct(name='output_data', mode='vertex',
                       data=vdofs, var_name=name, dofs=None,
                       mesh=mesh, levels=levels)

    out = convert_complex_output(out)

    return out

class Field(Struct):
    """
    Base class for finite element fields.

    Notes
    -----
    - Region can span over several groups -> different Aproximation
      instances
    - interps and hence node_descs are per region (must have single
      geometry!)
    - no two interps can be in a same group -> no two aps (with
      different regions) can be in a same group -> aps can be uniquely
      indexed with ig
    """
    _all = None

    @staticmethod
    def from_args(name, dtype, shape, region, approx_order=1,
                  space='H1', poly_space_base='lagrange'):
        """
        Create a Field subclass instance corresponding to a given space.

        Parameters
        ----------
        name : str
            The field name.
        dtype : numpy.dtype
            The field data type: float64 or complex128.
        shape : int/tuple/str
            The field shape: 1 or (1,) or 'scalar', space dimension (2, or
            (2,) or 3 or (3,)) or 'vector'. The field shape determines
            the shape of the FE base functions and can be different from
            a FieldVariable instance shape. (TODO)
        region : Region
            The region where the field is defined.
        approx_order : int/str
            The FE approximation order, e.g. 0, 1, 2, '1B' (1 with bubble).
        space : str
            The function space name.
        poly_space_base : str
            The name of polynomial space base.

        Notes
        -----
        Assumes one cell type for the whole region!
        """
        conf = Struct(name=name, dtype=dtype, shape=shape, region=region.name,
                      approx_order=approx_order, space=space,
                      poly_space_base=poly_space_base)
        return Field.from_conf(conf, {region.name : region})

    @staticmethod
    def from_conf(conf, regions):
        """
        Create a Field subclass instance based on the configuration.
        """
        if Field._all is None:
            import sfepy
            from sfepy.base.base import load_classes

            field_files = [ii for ii in sfepy.get_paths('sfepy/fem/fields*.py')
                           if 'fields_base.py' not in ii]
            Field._all = load_classes(field_files, [Field], ignore_errors=True,
                                      name_attr='family_name')
        table = Field._all

        space = conf.get_default_attr('space', 'H1')
        poly_space_base = conf.get_default_attr('poly_space_base', 'lagrange')

        key = space + '_' + poly_space_base

        approx_order = parse_approx_order(conf.approx_order)
        ao, force_bubble, discontinuous = approx_order

        if isinstance(conf.region, tuple):
            # Surface fields.
            region_name, kind = conf.region
            region = regions[region_name]

            cls = table[kind + '_' + key]
            obj = cls(conf.name, conf.dtype, conf.shape, region,
                      approx_order=approx_order[:2])


        else:
            # Volume fields.
            kind = 'volume'

            if discontinuous:
                cls = table[kind + '_' + key + '_discontinuous']

            else:
                cls = table[kind + '_' + key]

            obj = cls(conf.name, conf.dtype, conf.shape, regions[conf.region],
                      approx_order=approx_order[:2])

        return obj

    def __init__(self, name, dtype, shape, region, approx_order=1):
        """
        Create a Field.

        name : str
            The field name.
        dtype : numpy.dtype
            The field data type: float64 or complex128.
        shape : int/tuple/str
            The field shape: 1 or (1,) or 'scalar', space dimension (2, or
            (2,) or 3 or (3,)) or 'vector'. The field shape determines
            the shape of the FE base functions and can be different from
            a FieldVariable instance shape. (TODO)
        region : Region
            The region where the field is defined.
        approx_order : int/str
            The FE approximation order, e.g. 0, 1, 2, '1B' (1 with bubble).

        Notes
        -----
        Assumes one cell type for the whole region!
        """
        if isinstance(shape, basestr):
            try:
                shape = {'scalar' : (1,),
                         'vector' : (region.domain.shape.dim,)}[shape]
            except KeyError:
                raise ValueError('unsupported field shape! (%s)', shape)

        elif isinstance(shape, int):
            shape = (shape,)

        if not self._check_region(region):
            raise ValueError('unsuitable region for field %s! (%s)' %
                             (name, region.name))

        Struct.__init__(self,
                        name=name,
                        dtype=dtype,
                        shape=shape,
                        region=region)
        self.domain = self.region.domain
        self.igs = self.region.igs

        self.clear_dof_conns()

        self._set_approx_order(approx_order)
        self._setup_geometry()
        self._setup_kind()

        self._create_interpolant()
        self._setup_approximations()
        self._setup_global_base()
        self.setup_coors()
        self.clear_mappings(clear_all=True)

    def _set_approx_order(self, approx_order):
        """
        Set a uniform approximation order.
        """
        if isinstance(approx_order, tuple):
            self.approx_order = approx_order[0]
            self.force_bubble = approx_order[1]

        else:
            self.approx_order = approx_order
            self.force_bubble = False

    def _setup_kind(self):
        name = self.get_default_attr('family_name', None,
                                     'An abstract Field method called!')
        aux = name.split('_')
        self.space = aux[1]
        self.poly_space_base = aux[2]

    def _create_interpolant(self):
        name = '%s_%s_%s_%d%s' % (self.gel.name, self.space,
                                  self.poly_space_base, self.approx_order,
                                  'B' * self.force_bubble)
        self.interp = fea.Interpolant(name, self.gel, self.space,
                                      self.poly_space_base, self.approx_order,
                                      self.force_bubble)

    def _setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.Approximation(name, self.interp, self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def get_true_order(self):
        """
        Get the true approximation order depending on the reference
        element geometry.

        For example, for P1 (linear) approximation the true order is 1,
        while for Q1 (bilinear) approximation in 2D the true order is 2.
        """
        gel = self.gel
        if (gel.dim + 1) == gel.n_vertex:
            order = self.approx_order

        else:
            order = gel.dim * self.approx_order

        if self.force_bubble:
            bubble_order = gel.dim + 1
            order = max(order, bubble_order)

        return order

    def is_higher_order(self):
        """
        Return True, if the field's approximation order is greater than one.
        """
        return self.force_bubble or (self.approx_order > 1)

    def _setup_global_base(self):
        """
        Setup global DOF/base functions, their indices and connectivity of the
        field. Called methods implemented in subclasses.
        """
        self._setup_facet_orientations()

        self._init_econn()

        self.n_vertex_dof, self.vertex_remap = self._setup_vertex_dofs()
        self.vertex_remap_i = invert_remap(self.vertex_remap)

        aux = self._setup_edge_dofs()
        self.n_edge_dof, self.edge_dofs, self.edge_remap = aux

        aux = self._setup_face_dofs()
        self.n_face_dof, self.face_dofs, self.face_remap = aux

        aux = self._setup_bubble_dofs()
        self.n_bubble_dof, self.bubble_dofs, self.bubble_remaps = aux

        self.n_nod = self.n_vertex_dof + self.n_edge_dof \
                     + self.n_face_dof + self.n_bubble_dof

        self._setup_esurface()

    def _setup_esurface(self):
        """
        Setup extended surface entities (edges in 2D, faces in 3D),
        i.e. indices of surface entities into the extended connectivity.
        """
        node_desc = self.node_desc

        for ig, ap in self.aps.iteritems():
            gel = ap.interp.gel
            ap.efaces = gel.get_surface_entities().copy()

            nd = node_desc.edge
            if nd is not None:
                efs = []
                for eof in gel.get_edges_per_face():
                    efs.append(nm.concatenate([nd[ie] for ie in eof]))
                efs = nm.array(efs).squeeze()

                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
                ap.efaces = nm.hstack((ap.efaces, efs))

            efs = node_desc.face
            if efs is not None:
                efs = nm.array(efs).squeeze()

                if efs.ndim < 2:
                    efs = efs[:,nm.newaxis]
                ap.efaces = nm.hstack((ap.efaces, efs))

    def setup_coors(self, coors=None):
        """
        Setup coordinates of field nodes.
        """
        mesh = self.domain.mesh
        self.coors = nm.empty((self.n_nod, mesh.dim), nm.float64)

        if coors is None:
            coors = mesh.coors

        # Mesh vertex nodes.
        if self.n_vertex_dof:
            indx = self.vertex_remap_i
            self.coors[:self.n_vertex_dof] = nm.take(coors, indx, axis=0)

        for ig, ap in self.aps.iteritems():
            ap.eval_extra_coor(self.coors, coors)

    def get_vertices(self):
        """
        Return indices of vertices belonging to the field region.
        """
        return self.vertex_remap_i

    def get_dofs_in_region(self, region, merge=False, clean=False,
                           warn=False, igs=None):
        """
        Return indices of DOFs that belong to the given region.
        """
        if igs is None:
            igs = region.igs

        nods = []
        for ig in self.igs:
            if not ig in igs:
                nods.append(None)
                continue

            nn = self.get_dofs_in_region_group(region, ig)
            nods.append(nn)

        if merge:
            nods = [nn for nn in nods if nn is not None]
            nods = nm.unique(nm.hstack(nods))

        elif clean:
            for nn in nods[:]:
                if nn is None:
                    nods.remove(nn)
                    if warn is not None:
                        output(warn + ('%s' % region.name))

        return nods

    def _get_facet_dofs(self, facets, get_facets, remap, dofs, ig):
        ii = get_facets(ig)
        g_uid = facets.uid_i[ii]
        uid = remap[g_uid]

        return dofs[uid[uid >= 0]].ravel()

    def get_dofs_in_region_group(self, region, ig, merge=True):
        """
        Return indices of DOFs that belong to the given region and group.
        """
        node_desc = self.node_desc

        dofs = []

        vdofs = nm.empty((0,), dtype=nm.int32)
        if node_desc.vertex is not None:
            ii = region.get_vertices(ig)
            vdofs = self.vertex_remap[ii]
            vdofs = vdofs[vdofs >= 0]
        dofs.append(vdofs)

        edofs = nm.empty((0,), dtype=nm.int32)
        if node_desc.edge is not None:
            edofs = self._get_facet_dofs(self.domain.ed,
                                         region.get_edges,
                                         self.edge_remap,
                                         self.edge_dofs, ig)
        dofs.append(edofs)

        fdofs = nm.empty((0,), dtype=nm.int32)
        if node_desc.face is not None:
            fdofs = self._get_facet_dofs(self.domain.fa,
                                         region.get_faces,
                                         self.face_remap,
                                         self.face_dofs, ig)
        dofs.append(fdofs)

        bdofs = nm.empty((0,), dtype=nm.int32)
        if ((node_desc.bubble is not None)
            and region.can_cells and region.true_cells[ig]):
            ii = region.get_cells(ig)
            group_els = self.bubble_remaps[ig][ii]
            bdofs = self.bubble_dofs[ig][group_els[group_els >= 0]].ravel()
        dofs.append(bdofs)

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs

    def extend_dofs(self, dofs, fill_value=None):
        """
        Extend DOFs to the whole domain using the `fill_value`, or the
        smallest value in `dofs` if `fill_value` is None.
        """
        if fill_value is None:
            if dofs.shape[1] > 1: # Vector.
                fill_value = nm.amin(nm.abs(dofs))
            else: # Scalar.
                fill_value = nm.amin(dofs)

        if self.approx_order != 0:
            indx = self.get_vertices()

            n_nod = self.domain.shape.n_nod
            new_dofs = nm.empty((n_nod, dofs.shape[1]), dtype=self.dtype)
            new_dofs.fill(fill_value)
            new_dofs[indx] = dofs[:indx.size]

        else:
            new_dofs = extend_cell_data(dofs, self.domain, self.region,
                                        val=fill_value)

        return new_dofs

    def remove_extra_dofs(self, dofs):
        """
        Remove DOFs defined in higher order nodes (order > 1).
        """
        if self.approx_order != 0:
            new_dofs = dofs[:self.n_vertex_dof]

        else:
            new_dofs = dofs

        return new_dofs

    def linearize(self, dofs, min_level=0, max_level=1, eps=1e-4):
        """
        Linearize the solution for post-processing.

        Parameters
        ----------
        dofs : array, shape (n_nod, n_component)
            The array of DOFs reshaped so that each column corresponds
            to one component.
        min_level : int
            The minimum required level of mesh refinement.
        max_level : int
            The maximum level of mesh refinement.
        eps : float
            The relative tolerance parameter of mesh adaptivity.

        Returns
        -------
        mesh : Mesh instance
            The adapted, nonconforming, mesh.
        vdofs : array
            The DOFs defined in vertices of `mesh`.
        levels : array of ints
            The refinement level used for each element group.
        """
        assert_(dofs.ndim == 2)

        n_nod, dpn = dofs.shape

        assert_(n_nod == self.n_nod)
        assert_(dpn == self.shape[0])

        vertex_coors = self.coors[:self.n_vertex_dof, :]

        coors = []
        vdofs = []
        conns = []
        mat_ids = []
        levels = []
        offset = 0
        for ig, ap in self.aps.iteritems():
            ps = ap.interp.poly_spaces['v']
            gps = ap.interp.gel.interp.poly_spaces['v']
            group = self.domain.groups[ig]
            vertex_conn = ap.econn[:, :group.shape.n_ep]

            eval_dofs = get_eval_dofs(dofs, ap.econn, ps, ori=ap.ori)
            eval_coors = get_eval_coors(vertex_coors, vertex_conn, gps)

            (level, _coors, conn,
             _vdofs, _mat_ids) = create_output(eval_dofs, eval_coors,
                                               group.shape.n_el, ps,
                                               min_level=min_level,
                                               max_level=max_level, eps=eps)

            _mat_ids[:] = self.domain.mesh.mat_ids[ig][0]

            coors.append(_coors)
            vdofs.append(_vdofs)
            conns.append(conn + offset)
            mat_ids.append(_mat_ids)
            levels.append(level)

            offset += _coors.shape[0]

        coors = nm.concatenate(coors, axis=0)
        vdofs = nm.concatenate(vdofs, axis=0)
        mesh = Mesh.from_data('linearized_mesh', coors, None, conns, mat_ids,
                              self.domain.mesh.descs)

        return mesh, vdofs, levels

    def get_output_approx_order(self):
        """
        Get the approximation order used in the output file.
        """
        return min(self.approx_order, 1)

    def create_output(self, dofs, var_name, dof_names=None,
                      key=None, extend=True, fill_value=None,
                      linearization=None):
        """
        Convert the DOFs corresponding to the field to a dictionary of
        output data usable by Mesh.write().

        Parameters
        ----------
        dofs : array, shape (n_nod, n_component)
            The array of DOFs reshaped so that each column corresponds
            to one component.
        var_name : str
            The variable name corresponding to `dofs`.
        dof_names : tuple of str
            The names of DOF components.
        key : str, optional
            The key to be used in the output dictionary instead of the
            variable name.
        extend : bool
            Extend the DOF values to cover the whole domain.
        fill_value : float or complex
           The value used to fill the missing DOF values if `extend` is True.
        linearization : Struct or None
            The linearization configuration for higher order approximations.

        Returns
        -------
        out : dict
            The output dictionary.
        """
        linearization = get_default(linearization, Struct(kind='strip'))

        out = {}
        if linearization.kind is None:
            out[key] = Struct(name='output_data', mode='full',
                              data=dofs, var_name=var_name,
                              dofs=dof_names, field_name=self.name)

        elif ((not self.is_higher_order())
            or (linearization.kind == 'strip')):
            if extend:
                ext = self.extend_dofs(dofs, fill_value)

            else:
                ext = self.remove_extra_dofs(dofs)

            if ext is not None:
                approx_order = self.get_output_approx_order()

                if approx_order != 0:
                    # Has vertex data.
                    out[key] = Struct(name='output_data', mode='vertex',
                                      data=ext, var_name=var_name,
                                      dofs=dof_names)

                else:
                    ext.shape = (ext.shape[0], 1, ext.shape[1], 1)
                    out[key] = Struct(name='output_data', mode='cell',
                                      data=ext, var_name=var_name,
                                      dofs=dof_names)

        else:
            mesh, vdofs, levels = self.linearize(dofs,
                                                 linearization.min_level,
                                                 linearization.max_level,
                                                 linearization.eps)
            out[key] = Struct(name='output_data', mode='vertex',
                              data=vdofs, var_name=var_name, dofs=dof_names,
                              mesh=mesh, levels=levels)

        out = convert_complex_output(out)

        return out

    def clear_dof_conns(self):
        self.dof_conns = {}

    def setup_dof_conns(self, dof_conns, dpn, dc_type, region, is_trace=False):
        """Setup dof connectivities of various kinds as needed by terms."""
        if isinstance(dc_type, Struct):
            dct = dc_type.type

        else:
            dct = dc_type

        ##
        # Expand nodes into dofs.
        can_point = True
        for ig, ap in self.aps.iteritems():
            if ig not in region.igs: continue

            region_name = region.name # True region name.
            key = (self.name, dpn, region_name, dct, ig, is_trace)
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
                conn = sd.get_connectivity(is_trace=is_trace)
                dc = create_dof_conn(conn, dpn)
                self.dof_conns[key] = dc
                if is_trace:
                    trkey = key[:-1] + (False,)
                    if trkey in dof_conns:
                        self.dof_conns[trkey] = dof_conns[trkey]
                        continue
                    conn = sd.get_connectivity(is_trace=False)
                    dc = create_dof_conn(conn, dpn)
                    self.dof_conns[trkey] = dc

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

    def create_mesh(self, extra_nodes=True):
        """
        Create a mesh from the field region, optionally including the field
        extra nodes.
        """
        mesh = self.domain.mesh

        if self.approx_order != 0:
            conns, mat_ids, descs = [], [], []
            for ig, ap in self.aps.iteritems():
                group = self.domain.groups[ig]
                if extra_nodes:
                    conn = ap.econn
                else:
                    offset = group.shape.n_ep
                    conn = ap.econn[:,:offset]
                conns.append(conn)
                mat_ids.append(mesh.mat_ids[ig])
                descs.append(mesh.descs[ig])

            if extra_nodes:
                coors = self.coors

            else:
                coors = self.coors[:self.n_vertex_dof]

            mesh = Mesh.from_data(self.name, coors, None, conns,
                                  mat_ids, descs)

        return mesh

    def interp_to_qp(self, dofs):
        """
        Interpolate DOFs into quadrature points.

        The quadrature order is given by the field approximation order.

        Parameters
        ----------
        dofs : array
            The array of DOF values of shape `(n_nod, n_component)`.

        Returns
        -------
        data_qp : array
            The values interpolated into the quadrature points.
        integral : Integral
            The corresponding integral defining the quadrature points.
        """
        integral = Integral('i', order=self.approx_order)

        data_qp = []
        for ig, ap in self.aps.iteritems():
            bf = ap.get_base('v', False, integral)
            bf = bf[:,0,:].copy()

            vals = nm.dot(bf, dofs[ap.econn])
            vals = nm.swapaxes(vals, 0, 1)
            vals.shape = vals.shape + (1,)

            data_qp.append(vals)

        data_qp = nm.concatenate(data_qp, axis=0)

        return data_qp, integral

    def get_coor(self, nods=None):
        """
        Get coordinates of the field nodes.

        Parameters
        ----------
        nods : array, optional
           The indices of the required nodes. If not given, the
           coordinates of all the nodes are returned.
        """
        if nods is None:
            return self.coors
        else:
            return self.coors[nods]

    def clear_mappings(self, clear_all=False):
        """
        Clear current reference mappings.
        """
        self.mappings = {}
        if clear_all:
            self.mappings0 = {}

    def save_mappings(self):
        """
        Save current reference mappings to `mappings0` attribute.
        """
        self.mappings0 = self.mappings.copy()

    def create_mapping(self, ig, region, integral, integration):
        """
        Create a new reference mapping.
        """
        ap = self.aps[ig]

        out = ap.describe_geometry(self, integration, region, integral)
        return out

    def get_mapping(self, ig, region, integral, integration,
                    get_saved=False, return_key=False):
        """
        For given region, integral and integration type, get a reference
        mapping, i.e. jacobians, element volumes and base function
        derivatives for Volume-type geometries, and jacobians, normals
        and base function derivatives for Surface-type geometries
        corresponding to the field approximation.

        The mappings are cached in the field instance in `mappings`
        attribute. The mappings can be saved to `mappings0` using
        `Field.save_mappings`. The saved mapping can be retrieved by
        passing `get_saved=True`. If the required (saved) mapping
        is not in cache, a new one is created.

        Returns
        -------
        geo : VolumeGeometry or SurfaceGeometry instance
            The geometry object that describes the mapping.
        mapping : VolumeMapping or SurfaceMapping instance
            The mapping.
        key : tuple
            The key of the mapping in `mappings` or `mappings0`.
        """
        ap = self.aps[ig]
        # Share full group mappings.
        shape = self.domain.groups[ig].shape
        if ((region.shape[ig].n_vertex == shape.n_vertex)
            and (region.shape[ig].n_cell == shape.n_el)):
            region_name = ig

        else:
            region_name = region.name

        key = (integral.get_key(), region_name, ig, integration)

        # out is (geo, mapping) tuple.
        if get_saved:
            out = self.mappings0.get(key, None)

        else:
            out = self.mappings.get(key, None)

        if out is None:
            out = ap.describe_geometry(self, integration, region, integral,
                                       return_mapping=True)
            self.mappings[key] = out

        if return_key:
            out = out + (key,)

        return out

class VolumeField(Field):
    """
    Finite element field base class over volume elements (element dimension
    equals space dimension).
    """

    def _check_region(self, region):
        """
        Check whether the `region` can be used for the
        field. Non-surface fields require the region to span whole
        element groups.

        Returns
        -------
        ok : bool
            True if the region is usable for the field.
        """
        ok = True
        for ig in region.igs:
            shape = region.domain.groups[ig].shape
            if region.shape[ig].n_vertex < shape.n_vertex:
                ok = False
                break

        return ok

    def _setup_geometry(self):
        """
        Setup the field region geometry.
        """
        self.gel = self.domain.groups[self.region.igs[0]].gel

    def _create_interpolant(self):
        name = '%s_%s_%s_%d%s' % (self.gel.name, self.space,
                                  self.poly_space_base, self.approx_order,
                                  'B' * self.force_bubble)
        self.interp = fea.Interpolant(name, self.gel, self.space,
                                      self.poly_space_base, self.approx_order,
                                      self.force_bubble)

    def _setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.Approximation(name, self.interp, self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def _init_econn(self):
        """
        Initialize the extended DOF connectivity.
        """
        for ig, ap in self.aps.iteritems():
            n_ep = ap.n_ep['v']
            n_cell = self.region.get_n_cells(ig)
            ap.econn = nm.zeros((n_cell, n_ep), nm.int32)

    def _setup_vertex_dofs(self):
        """
        Setup vertex DOF connectivity.
        """
        if self.node_desc.vertex is None:
            return 0, None

        region = self.region

        all_vertices = region.get_vertices_of_cells()
        remap = prepare_remap(all_vertices, region.n_v_max)
        n_dof = all_vertices.shape[0]

        ##
        # Remap vertex node connectivity to field-local numbering.
        for ig, ap in self.aps.iteritems():
            group = self.domain.groups[ig]
            offset = group.shape.n_ep
            cells = region.get_cells(ig)
            ap.econn[:,:offset] = nm.take(remap,
                                          nm.take(group.conn, cells, axis=0))

        return n_dof, remap

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type

        if geometry != None:
            geometry_flag = 'surface' in geometry
        else:
            geometry_flag = False

        if (dct == 'surface') or (geometry_flag):
            reg = info.get_region()
            # Calls reg.select_cells_of_surface(reset=False)...
            self.domain.create_surface_group(reg)
            self._setup_surface_data(reg, is_trace)

        elif dct == 'edge':
            raise NotImplementedError('dof connectivity type %s' % dct)

        elif dct == 'point':
            self._setup_point_data(self, info.region)

        elif dct not in ('volume', 'scalar'):
            raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def _setup_surface_data(self, region, is_trace=False):
        for ig, ap in self.aps.iteritems():
            if ig not in region.igs: continue
            if region.name not in ap.surface_data:
                ap.setup_surface_data(region)

        for ig, ap in self.aps.iteritems():
            if region.name in ap.surface_data and is_trace:
                sd = ap.surface_data[region.name]
                sd.setup_mirror_connectivity(region)

    def _setup_point_data(self, field, region):
        # Point data only in the first group to avoid multiple
        # assembling of nodes on group boundaries.
        ap = self.aps[self.igs[0]]
        if region.name not in ap.point_data:
            ap.setup_point_data(field, region)

    def setup_dof_conns(self, dof_conns, dpn, dc_type, region, is_trace=False):
        """Setup dof connectivities of various kinds as needed by terms."""
        if isinstance(dc_type, Struct):
            dct = dc_type.type

        else:
            dct = dc_type

        ##
        # Expand nodes into dofs.
        can_point = True
        for ig, ap in self.aps.iteritems():
            if ig not in region.igs: continue

            region_name = region.name # True region name.
            key = (self.name, dpn, region_name, dct, ig, is_trace)
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
                conn = sd.get_connectivity(is_trace=is_trace)
                dc = create_dof_conn(conn, dpn)
                self.dof_conns[key] = dc
                if is_trace:
                    trkey = key[:-1] + (False,)
                    if trkey in dof_conns:
                        self.dof_conns[trkey] = dof_conns[trkey]
                        continue
                    conn = sd.get_connectivity(is_trace=False)
                    dc = create_dof_conn(conn, dpn)
                    self.dof_conns[trkey] = dc

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

    def average_qp_to_vertices(self, data_qp, integral):
        """
        Average data given in quadrature points in region elements into
        region vertices.

        .. math::
           u_n = \sum_e (u_{e,avg} * volume_e) / \sum_e volume_e
               = \sum_e \int_{volume_e} u / \sum volume_e
        """
        domain = self.domain
        if domain.shape.n_el != data_qp.shape[0]:
            msg = 'incomatible shape! (%d == %d)' % (domain.shape.n_el,
                                                     data_qp.shape[0])
            raise ValueError(msg)

        n_vertex = domain.shape.n_nod
        dim = data_qp.shape[2]

        nod_vol = nm.zeros((n_vertex,), dtype=nm.float64)
        data_vertex = nm.zeros((n_vertex, dim), dtype=nm.float64)
        for ig, ap in self.aps.iteritems():
            vg = ap.describe_geometry(self, 'volume', ap.region, integral)

            volume = nm.squeeze(vg.volume)
            iels = ap.region.cells[ig]

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

        return data_vertex

class SurfaceField(Field):
    """
    Finite element field base class over surface (element dimension is one
    less than space dimension).
    """

    def _check_region(self, region):
        """
        Check whether the `region` can be used for the
        field.

        Returns
        -------
        ok : bool
            True if the region is usable for the field.
        """
        region.setup_face_indices()

        ok = True
        for ig in region.igs:
            n_cell = region.get_n_cells(ig, True)
            if n_cell == 0:
                ok = False
                break

        return ok

    def _setup_geometry(self):
        """
        Setup the field region geometry.
        """
        self.gel = self.domain.groups[self.region.igs[0]].gel.surface_facet
        if self.gel is None:
            raise ValueError('element group has no surface!')

    def _create_interpolant(self):
        name = '%s_%s_%s_%d%s' % (self.gel.name, self.space,
                                  self.poly_space_base, self.approx_order,
                                  'B' * self.force_bubble)
        self.interp = fea.SurfaceInterpolant(name, self.gel, self.space,
                                             self.poly_space_base,
                                             self.approx_order,
                                             self.force_bubble)

    def _setup_approximations(self):
        self.aps = {}
        self.aps_by_name = {}
        for ig in self.igs:
            name = self.interp.name + '_%s_ig%d' % (self.region.name, ig)
            ap = fea.SurfaceApproximation(name, self.interp, self.region, ig)
            self.aps[ig] = ap
            self.aps_by_name[ap.name] = ap

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type

        if dct != 'surface':
            msg = "dof connectivity type must be 'surface'! (%s)" % dct
            raise ValueError(msg)

        reg = info.get_region()
        reg.select_cells_of_surface(reset=False)

        for ig, ap in self.aps.iteritems():
            if ig not in reg.igs: continue

            if reg.name not in ap.surface_data:
                # Defined in setup_vertex_dofs()
                msg = 'no surface data of surface field! (%s)' % reg.name
                raise ValueError(msg)

        for ig, ap in self.aps.iteritems():
            if reg.name in ap.surface_data and is_trace:
                sd = ap.surface_data[reg.name]
                sd.setup_mirror_connectivity(reg)

    def _init_econn(self):
        """
        Initialize the extended DOF connectivity.
        """
        for ig, ap in self.aps.iteritems():
            n_ep = ap.n_ep['v']
            n_cell = self.region.get_n_cells(ig, True)
            ap.econn = nm.zeros((n_cell, n_ep), nm.int32)

    def _setup_vertex_dofs(self):
        """
        Setup vertex DOF connectivity.
        """
        if self.node_desc.vertex is None:
            return 0, None

        region = self.region

        remap = prepare_remap(region.all_vertices, region.n_v_max)
        n_dof = region.all_vertices.shape[0]

        ##
        # Remap vertex node connectivity to field-local numbering.
        for ig, ap in self.aps.iteritems():
            group = self.domain.groups[ig]
            faces = group.gel.get_surface_entities()
            aux = FESurface('aux', region, faces, group.conn, ig)
            ap.econn[:,:aux.n_fp] = aux.leconn
            ap.surface_data[region.name] = aux

        return n_dof, remap

    def _setup_bubble_dofs(self):
        """
        Setup bubble DOF connectivity.
        """
        return 0, None, None

    def setup_dof_conns(self, dof_conns, dpn, dc_type, region, is_trace=False):
        """Setup dof connectivities of various kinds as needed by terms."""
        dct = dc_type.type

        if dct != 'surface':
            msg = "dof connectivity type must be 'surface'! (%s)" % dct
            raise ValueError(msg)

        ##
        # Expand nodes into dofs.
        for ig, ap in self.aps.iteritems():
            if ig not in region.igs: continue

            region_name = region.name # True region name.
            key = (self.name, dpn, region_name, dct, ig, is_trace)
            if key in dof_conns:
                self.dof_conns[key] = dof_conns[key]
                continue

            sd = ap.surface_data[region_name]
            conn = sd.get_connectivity(local=True, is_trace=is_trace)
            dc = create_dof_conn(conn, dpn)
            self.dof_conns[key] = dc

        dof_conns.update(self.dof_conns)

    def average_qp_to_vertices(self, data_qp, integral):
        """
        Average data given in quadrature points in region elements into
        region vertices.

        .. math::
           u_n = \sum_e (u_{e,avg} * area_e) / \sum_e area_e
               = \sum_e \int_{area_e} u / \sum area_e
        """
        region = self.region

        n_cells = region.get_n_cells(None, True)
        if n_cells != data_qp.shape[0]:
            msg = 'incomatible shape! (%d == %d)' % (n_cells,
                                                     data_qp.shape[0])
            raise ValueError(msg)

        n_vertex = len(region.all_vertices)
        nc = data_qp.shape[2]

        nod_vol = nm.zeros((n_vertex,), dtype=nm.float64)
        data_vertex = nm.zeros((n_vertex, nc), dtype=nm.float64)
        offset = 0
        for ig, ap in self.aps.iteritems():
            sg = ap.describe_geometry(self, 'surface', ap.region, integral)

            area = nm.squeeze(sg.area)
            n_cells = region.get_n_cells(ig, True)
            iels = offset + nm.arange(n_cells, dtype=nm.int32)
            offset += n_cells

            data_e = nm.zeros((area.shape[0], 1, nc, 1), dtype=nm.float64)
            sg.integrate(data_e, data_qp[iels])

            ir = nm.arange(nc, dtype=nm.int32)

            sd = self.domain.surface_groups[ig][region.name]
            # Should be vertex connectivity!
            conn = sd.get_connectivity(local=True)
            for ii, cc in enumerate(conn):
                # Assumes unique nodes in cc!
                ind2, ind1 = nm.meshgrid(ir, cc)
                data_vertex[ind1,ind2] += data_e[iels[ii],0,:,0]
                nod_vol[cc] += area[ii]
        data_vertex /= nod_vol[:,nm.newaxis]

        return data_vertex
