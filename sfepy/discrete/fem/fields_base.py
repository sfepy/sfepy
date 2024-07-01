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
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import output, get_default, assert_
from sfepy.base.base import Struct
from sfepy.base.timing import Timer
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete import PolySpace
from sfepy.discrete.fem.mesh import Mesh
from sfepy.discrete.fem.meshio import convert_complex_output
from sfepy.discrete.fem.utils import (extend_cell_data, prepare_remap,
                                      invert_remap, get_min_value)
from sfepy.discrete.fem.mappings import FEMapping
from sfepy.discrete.fem.fe_surface import FESurface, FEPhantomSurface
from sfepy.discrete.integrals import Integral
from sfepy.discrete.fem.linearizer import (get_eval_dofs, get_eval_coors,
                                           create_output)

def _find_geometry(region):
    cmesh = region.cmesh
    if region.kind == 'cell':
        ct = cmesh.cell_types
        for _gel in region.domain.geom_els.values():
            if (ct[region.cells] == cmesh.key_to_index[_gel.name]).all():
                gel = _gel
                break
        else:
            raise ValueError(f'region {region.name} of contains multiple'
                             ' reference geometries!')

        is_surface = False

    elif region.kind == 'facet':
        for _gel in region.domain.geom_els.values():
            gel = _gel.surface_facet
            break

        if gel is None:
            raise ValueError('cells with no surface!')

        is_surface = True

    else:
        raise ValueError('cannot find geometry element for region'
                         f' {region.name} of kind {region.kind}, '
                         f'"region.kind" must be "cell" or "facet"!')

    return gel, is_surface

def set_mesh_coors(domain, fields, coors, update_fields=False, actual=False,
                   clear_all=True, extra_dofs=False):
    if actual:
        if not hasattr(domain.mesh, 'coors_act'):
            domain.mesh.coors_act = nm.zeros_like(domain.mesh.coors)
        domain.mesh.coors_act[:] = coors[:domain.mesh.n_nod]
    else:
        domain.cmesh.coors[:] = coors[:domain.mesh.n_nod]

    if update_fields:
        for field in fields.values():
            field.set_coors(coors, extra_dofs=extra_dofs)
            field.clear_mappings(clear_all=clear_all)


def eval_nodal_coors(coors, mesh_coors, region, poly_space, geom_poly_space,
                     econn, only_extra=True):
    """
    Compute coordinates of nodes corresponding to `poly_space`, given
    mesh coordinates and `geom_poly_space`.
    """
    if only_extra:
        iex = (poly_space.nts[:, 0] > 0).nonzero()[0]
        if iex.shape[0] == 0:
            return

        qp_coors = poly_space.node_coors[iex, :]
        econn = econn[:, iex].copy()

    else:
        qp_coors = poly_space.node_coors

    ##
    # Evaluate geometry interpolation basis functions in (extra) nodes.
    bf = geom_poly_space.eval_basis(qp_coors)
    bf = bf[:, 0, :].copy()

    ##
    # Evaluate extra coordinates with 'bf'.
    cmesh = region.cmesh
    conn = cmesh.get_incident(0, region.cells, region.tdim)
    conn.shape = (econn.shape[0], -1)

    ecoors = nm.dot(bf, mesh_coors[conn])
    coors[econn] = nm.swapaxes(ecoors, 0, 1)


def _interp_to_faces(vertex_vals, bfs, faces):
    dim = vertex_vals.shape[1]
    n_face = faces.shape[0]
    n_qp = bfs.shape[0]

    faces_vals = nm.zeros((n_face, n_qp, dim), nm.float64)
    for ii, face in enumerate(faces):
        vals = vertex_vals[face, :dim]
        faces_vals[ii, :, :] = nm.dot(bfs[:, 0, :], vals)

    return(faces_vals)


def get_eval_expression(expression,
                        fields, materials, variables,
                        functions=None, mode='eval', term_mode=None,
                        extra_args=None, verbose=True, kwargs=None):
    """
    Get the function for evaluating an expression given a list of elements,
    and reference element coordinates.
    """
    from sfepy.discrete.evaluate import eval_in_els_and_qp

    def _eval(iels, coors):
        val = eval_in_els_and_qp(expression, iels, coors,
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

    ps = field.poly_space
    gps = field.gel.poly_space
    vertex_conn = field.econn[:, :field.gel.n_vertex]

    eval_dofs = get_eval_expression(expression,
                                    fields, materials, variables,
                                    functions=functions,
                                    mode=mode, extra_args=extra_args,
                                    verbose=verbose, kwargs=kwargs)
    eval_coors = get_eval_coors(vertex_coors, vertex_conn, gps)

    (level, coors, conn,
     vdofs, mat_ids) = create_output(eval_dofs, eval_coors,
                                     vertex_conn.shape[0], ps,
                                     min_level=min_level,
                                     max_level=max_level, eps=eps)

    mesh = Mesh.from_data('linearized_mesh', coors, None, [conn], [mat_ids],
                          field.domain.mesh.descs)

    out = {}
    out[name] = Struct(name='output_data', mode='vertex',
                       data=vdofs, var_name=name, dofs=None,
                       mesh=mesh, level=level)

    out = convert_complex_output(out)

    return out


class FEField(Field):
    """
    Base class for finite element fields.

    Notes
    -----
    - interps and hence node_descs are per region (must have single
      geometry!)

    Field shape information:

    - ``shape`` - the shape of the basis functions in a point
    - ``n_components`` - the number of DOFs per FE node
    - ``val_shape`` - the shape of field value (the product of DOFs and
      basis functions) in a point
    """

    def __init__(self, name, dtype, shape, region, approx_order=1):
        """
        Create a finite element field.

        Parameters
        ----------
        name : str
            The field name.
        dtype : numpy.dtype
            The field data type: float64 or complex128.
        shape : int/tuple/str
            The field shape: 1 or (1,) or 'scalar', space dimension (2, or (2,)
            or 3 or (3,)) or 'vector', or a tuple. The field shape determines
            the shape of the FE basis functions and is related to the number of
            components of variables and to the DOF per node count, depending
            on the field kind.
        region : Region
            The region where the field is defined.
        approx_order : int or tuple
            The FE approximation order. The tuple form is (order, has_bubble),
            e.g. (1, True) means order 1 with a bubble function.

        Notes
        -----
        Assumes one cell type for the whole region!
        """
        field_dim = region.field_dim if hasattr(region, 'field_dim')\
            else region.domain.shape.dim
        shape = parse_shape(shape, field_dim)
        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region)
        self.domain = self.region.domain
        self.cmesh = self.region.cmesh

        self._set_approx_order(approx_order)
        self.gel, self.is_surface = _find_geometry(self.region)
        self._setup_kind()
        self._setup_shape()

        self.extra_data = {}
        self.ori = None
        self._create_interpolant()
        self._setup_global_basis()
        self.setup_coors()
        self.clear_mappings(clear_all=True)
        self.clear_qp_basis()
        self.basis_transform = None
        self.econn0 = None
        self.unused_dofs = None
        self.stored_subs = None

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

    def _create_interpolant(self):
        name = '%s_%s_%s_%d%s' % (self.gel.name, self.space,
                                  self.poly_space_basis, self.approx_order,
                                  'B' * self.force_bubble)
        ps = PolySpace.any_from_args(name, self.gel, self.approx_order,
                                     basis=self.poly_space_basis,
                                     force_bubble=self.force_bubble)
        self.poly_space = ps

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

    def _setup_global_basis(self):
        """
        Setup global DOF/basis functions, their indices and connectivity of the
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
        self.n_bubble_dof, self.bubble_dofs, self.bubble_remap = aux

        self.n_nod = (self.n_vertex_dof + self.n_edge_dof
                      + self.n_face_dof + self.n_bubble_dof)

        self._setup_esurface()

    def _init_econn(self):
        """
        Initialize the extended DOF connectivity.
        """
        n_ep = self.poly_space.n_nod
        n_cell = self.region.get_n_cells(is_surface=self.is_surface)
        self.econn = nm.zeros((n_cell, n_ep), nm.int32)

    def _setup_vertex_dofs(self):
        """
        Setup vertex DOF connectivity.
        """
        if self.node_desc.vertex is None:
            return 0, None

        region = self.region

        remap = prepare_remap(region.vertices, region.n_v_max)
        n_dof = region.vertices.shape[0]

        # Remap vertex node connectivity to field-local numbering.
        if self.is_surface:
            aux = FESurface.from_region('aux', region)
            self.econn[:, :aux.n_fp] = aux.leconn
            self.extra_data[f'sd_{region.name}'] = aux
        else:
            conn = self.domain.get_conn(tdim=region.tdim, cells=region.cells)
            self.econn[:, :conn.shape[1]] = nm.take(remap, conn)

        return n_dof, remap

    def _setup_esurface(self):
        """
        Setup extended surface entities (edges in 2D, faces in 3D),
        i.e. indices of surface entities into the extended connectivity.
        """
        node_desc = self.node_desc

        gel = self.gel
        self.efaces = gel.get_surface_entities().copy()

        nd = node_desc.edge
        if nd is not None:
            efs = []
            for eof in gel.get_edges_per_face():
                efs.append(nm.concatenate([nd[ie] for ie in eof]))
            efs = nm.array(efs).squeeze()

            if efs.ndim < 2:
                efs = efs[:, nm.newaxis]
            self.efaces = nm.hstack((self.efaces, efs))

        efs = node_desc.face
        if efs is not None:
            efs = nm.array(efs).squeeze()

            if efs.ndim < 2:
                efs = efs[:, nm.newaxis]
            self.efaces = nm.hstack((self.efaces, efs))

        if gel.dim == 3:
            self.eedges = gel.edges.copy()
            efs = node_desc.edge
            if efs is not None:
                efs = nm.array(efs).squeeze()

                if efs.ndim < 2:
                    efs = efs[:, nm.newaxis]
                self.eedges = nm.hstack((self.eedges, efs))

    def set_coors(self, coors, extra_dofs=False):
        """
        Set coordinates of field nodes.
        """
        # Mesh vertex nodes.
        if self.n_vertex_dof:
            indx = self.vertex_remap_i
            self.coors[:self.n_vertex_dof] = nm.take(coors,
                                                     indx.astype(nm.int32),
                                                     axis=0)

        n_ex_dof = self.n_bubble_dof + self.n_edge_dof + self.n_face_dof

        # extra nodes
        if n_ex_dof:
            if extra_dofs:
                if self.n_nod != coors.shape[0]:
                    raise NotImplementedError
                self.coors[:] = coors
            else:
                gps = self.gel.poly_space
                ps = self.poly_space
                eval_nodal_coors(self.coors, coors, self.region,
                                 ps, gps, self.econn)

    def setup_coors(self):
        """
        Setup coordinates of field nodes.
        """
        mesh = self.domain.mesh
        self.coors = nm.empty((self.n_nod, mesh.dim), nm.float64)
        self.set_coors(mesh.coors)

    def get_vertices(self):
        """
        Return indices of vertices belonging to the field region.
        """
        return self.vertex_remap_i

    def _get_facet_dofs(self, rfacets, remap, dofs):
        facets = remap[rfacets]

        return dofs[facets[facets >= 0]].ravel()

    def get_data_shape(self, integral, integration='cell', region_name=None):
        """
        Get element data dimensions.

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
        data_shape : 4 ints
            The `(n_el, n_qp, dim, n_en)` for volume shape kind,
            `(n_fa, n_qp, dim, n_fn)` for surface shape kind and
            `(n_nod, 0, 0, 1)` for point shape kind.

        Notes
        -----
        Integration modes:
        - 'cell': integration over cells/elements
        - 'facet': integration over cell facets (faces, edges)
        - 'facet_extra': same as 'facet' but also the normal derivatives
          are evaluated
        - 'point': point integration
        - 'custom': user defined integration

        Dimensions:
        - `n_el`, `n_fa` = number of elements/facets
        - `n_qp` = number of quadrature points per element/facet
        - `dim` = spatial dimension
        - `n_en`, `n_fn` = number of element/facet nodes
        - `n_nod` = number of element nodes
        """
        region = self.domain.regions[region_name]
        shape = region.shape
        dim = region.field_dim if hasattr(region, 'field_dim') else region.dim

        if integration is None:
            integration == region.kind

        if 'facet' in integration:
            name = f'sd_{region_name}'
            if name not in self.extra_data:
                reg = self.domain.regions[region_name]
                self.domain.create_surface_group(reg)
                self.setup_surface_data(reg, None)

            sd = self.extra_data[name]

            # This works also for surface fields.
            key = sd.face_type
            weights = self.get_qp(key, integral).weights
            n_qp = weights.shape[0]

            if integration == 'facet_extra':
                data_shape = (sd.n_fa, n_qp, dim, self.econn.shape[1])
            else:
                data_shape = (sd.n_fa, n_qp, dim, sd.n_fp)

        elif (integration == 'cell' and self.region.tdim > 1 and
              region.tdim == 1):
            data_shape = (shape.n_cell, 0, dim, 2)  # bar elements

        elif integration in ('cell', 'custom'):
            _, weights = integral.get_qp(self.gel.name)
            n_qp = weights.shape[0]

            data_shape = (shape.n_cell, n_qp, dim, self.econn.shape[1])

        elif integration == 'point':
            dofs = self.get_dofs_in_region(region, merge=True)
            data_shape = (dofs.shape[0], 0, 0, 1)

        else:
            raise NotImplementedError('unsupported integration type! (%s)'
                                      % integration)

        return data_shape

    def get_dofs_in_region(self, region, merge=True):
        """
        Return indices of DOFs that belong to the given region.
        """
        node_desc = self.node_desc

        dofs = []

        vdofs = nm.empty((0,), dtype=nm.int32)
        if node_desc.vertex is not None:
            vdofs = self.vertex_remap[region.vertices]
            vdofs = vdofs[vdofs >= 0]
        dofs.append(vdofs)

        edofs = nm.empty((0,), dtype=nm.int32)
        if node_desc.edge is not None:
            edofs = self._get_facet_dofs(region.edges,
                                         self.edge_remap,
                                         self.edge_dofs)
        dofs.append(edofs)

        fdofs = nm.empty((0,), dtype=nm.int32)
        if node_desc.face is not None:
            fdofs = self._get_facet_dofs(region.faces,
                                         self.face_remap,
                                         self.face_dofs)
        dofs.append(fdofs)

        bdofs = nm.empty((0,), dtype=nm.int32)
        if (node_desc.bubble is not None) and region.has_cells():
            els = self.bubble_remap[region.cells]
            bdofs = self.bubble_dofs[els[els >= 0]].ravel()
        dofs.append(bdofs)

        if merge:
            dofs = nm.concatenate(dofs)

        return dofs

    def clear_qp_basis(self):
        """
        Remove cached quadrature points and basis functions.
        """
        self.qp_coors = {}
        self.bf = {}

    def get_qp(self, key, integral):
        """
        Get quadrature points and weights corresponding to the given key and
        integral. The key is 'v', 's#' or 'b#', where # is the number of face
        vertices. For 'b#', the quadrature must already be created by calling
        :func:`FEField.create_bqp()`, usually through
        :func:`FEField.create_mapping()`.
        """
        qpkey = (integral.order, key)

        if qpkey not in self.qp_coors:
            if key[0] == 'b':
                raise ValueError(f'the quadrature "{qpkey}" does not exist!')

            if (key[0] == 's') and not self.is_surface:
                dim = self.gel.dim - 1
                if isinstance(self.gel.surface_facet, dict):
                    n_fp = int(key[1:])
                else:
                    n_fp = self.gel.surface_facet.n_vertex
                geometry = '%d_%d' % (dim, n_fp)

            else:
                geometry = self.gel.name

            vals, weights = integral.get_qp(geometry)
            self.qp_coors[qpkey] = Struct(vals=vals, weights=weights)

        return self.qp_coors[qpkey]

    def substitute_dofs(self, subs, restore=False):
        """
        Perform facet DOF substitutions according to `subs`.

        Modifies `self.econn` in-place and sets `self.econn0`,
        `self.unused_dofs` and `self.basis_transform`.
        """
        if restore and (self.stored_subs is not None):
            self.econn0 = self.econn
            self.econn, self.unused_dofs, basis_transform = self.stored_subs

        else:
            if subs is None:
                self.econn0 = self.econn
                return

            else:
                self.econn0 = self.econn.copy()

            self._substitute_dofs(subs)

            self.unused_dofs = nm.setdiff1d(self.econn0, self.econn)

            basis_transform = self._eval_basis_transform(subs)

        self.set_basis_transform(basis_transform)

    def restore_dofs(self, store=False):
        """
        Undoes the effect of :func:`FEField.substitute_dofs()`.
        """
        if self.econn0 is None:
            raise ValueError('no original DOFs to restore!')

        if store:
            self.stored_subs = (self.econn,
                                self.unused_dofs,
                                self.basis_transform)

        else:
            self.stored_subs = None

        self.econn = self.econn0
        self.econn0 = None
        self.unused_dofs = None
        self.basis_transform = None

    def set_basis_transform(self, transform):
        """
        Set local element basis transformation.

        The basis transformation is applied in :func:`FEField.eval_basis()` and
        :func:`FEField.create_mapping()`.

        Parameters
        ----------
        transform : array, shape `(n_cell, n_ep, n_ep)`
            The array with `(n_ep, n_ep)` transformation matrices for each cell
            in the field's region, where `n_ep` is the number of element DOFs.
        """
        self.basis_transform = transform

    def restore_substituted(self, vec):
        """
        Restore values of the unused DOFs using the transpose of the applied
        basis transformation.
        """
        if (self.econn0 is None) or (self.basis_transform is None):
            raise ValueError('no original DOF values to restore!!')

        vec = vec.reshape((self.n_nod, self.n_components)).copy()
        evec = vec[self.econn]

        vec[self.econn0] = nm.einsum('cji,cjk->cik', self.basis_transform,
                                     evec, optimize=True)

        return vec.ravel()

    def eval_basis(self, key, derivative, integral, iels=None,
                   from_geometry=False, basis_only=True):
        qp = self.get_qp(key, integral)

        if from_geometry:
            ps = self.gel.poly_space

        else:
            ps = self.poly_space

        _key = key if not from_geometry else 'g' + key
        bf_key = (integral.order, _key, derivative)

        if bf_key not in self.bf:
            ori = self.ori
            self.bf[bf_key] = ps.eval_basis(qp.vals, diff=derivative, ori=ori,
                                            transform=self.basis_transform)

        bf = self.bf[bf_key]
        if iels is not None and bf.ndim == 4:
            bf = bf[iels]

        if basis_only:
            return bf

        else:
            return bf, qp.weights

    def create_bqp(self, region_name, integral):
        gel = self.gel

        sd = self.extra_data[f'sd_{region_name}']
        bqpkey = (integral.order, sd.bkey)
        if bqpkey not in self.qp_coors:
            qp = self.get_qp(sd.face_type, integral)

            ps_s = self.gel.surface_facet.poly_space
            bf_s = ps_s.eval_basis(qp.vals)

            coors, faces = gel.coors, gel.get_surface_entities()

            vals = _interp_to_faces(coors, bf_s, faces)
            self.qp_coors[bqpkey] = Struct(name='BQP_%s' % sd.bkey,
                                           vals=vals, weights=qp.weights)

    def create_bqp_key(self, integral, bkey):
        gel = self.gel
        sd_bkey, face_type = f'b{bkey}', f's{bkey}'

        bqpkey = (integral.order, sd_bkey)
        if bqpkey not in self.qp_coors:
            qp = self.get_qp(face_type, integral)

            ps_s = self.gel.surface_facet[bkey].poly_space
            bf_s = ps_s.eval_basis(qp.vals)

            coors, faces = gel.coors, gel.get_surface_entities()

            vals = _interp_to_faces(coors, bf_s,
                                    faces[:, :bf_s.shape[-1]])
            self.qp_coors[bqpkey] = Struct(name=f'BQP_{sd_bkey}',
                                           vals=vals, weights=qp.weights)

    def extend_dofs(self, dofs, fill_value=None):
        """
        Extend DOFs to the whole domain using the `fill_value`, or the
        smallest value in `dofs` if `fill_value` is None.
        """
        if fill_value is None:
            if nm.isrealobj(dofs):
                fill_value = get_min_value(dofs)

            else:
                # Complex values - treat real and imaginary parts separately.
                fill_value = get_min_value(dofs.real)
                fill_value += 1j * get_min_value(dofs.imag)

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

        ps = self.poly_space
        gps = self.gel.poly_space

        vertex_conn = self.econn[:, :self.gel.n_vertex]

        eval_dofs = get_eval_dofs(dofs, self.econn, ps, ori=self.ori)
        eval_coors = get_eval_coors(vertex_coors, vertex_conn, gps)

        (level, coors, conn,
         vdofs, mat_ids) = create_output(eval_dofs, eval_coors,
                                         vertex_conn.shape[0], ps,
                                         min_level=min_level,
                                         max_level=max_level, eps=eps)

        mesh = Mesh.from_data('linearized_mesh', coors, None, [conn],
                              [mat_ids], self.domain.mesh.descs)

        return mesh, vdofs, level

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
        reg_name = self.region.name
        if linearization.kind is None:
            out[key] = Struct(name='output_data', mode='full',
                              data=dofs, var_name=var_name,
                              region_name=reg_name,
                              dofs=dof_names, field_name=self.name)

        elif linearization.kind == 'strip':
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
                                      dofs=dof_names, region_name=reg_name)

                else:
                    ext.shape = (ext.shape[0], 1, ext.shape[1], 1)
                    out[key] = Struct(name='output_data', mode='cell',
                                      data=ext, var_name=var_name,
                                      dofs=dof_names, region_name=reg_name)

        else:
            mesh, vdofs, levels = self.linearize(dofs,
                                                 linearization.min_level,
                                                 linearization.max_level,
                                                 linearization.eps)
            out[key] = Struct(name='output_data', mode='vertex',
                              data=vdofs, var_name=var_name, dofs=dof_names,
                              mesh=mesh, levels=levels, region_name=reg_name)

        out = convert_complex_output(out)

        return out

    def create_mesh(self, extra_nodes=True):
        """
        Create a mesh from the field region, optionally including the field
        extra nodes.
        """
        mesh = self.domain.mesh

        if self.approx_order != 0:
            if extra_nodes:
                conn = self.econn

            else:
                conn = self.econn[:, :self.gel.n_vertex]

            conns = [conn]
            mat_ids = [self.cmesh.cell_groups]
            tdim = self.cmesh.tdim
            descs = f'{tdim}_{conn.shape[1]}'
            if descs not in mesh.descs:
                msg = f'element type {descs} not in mesh! ({mesh.descs})'
                raise ValueError(msg)

            if extra_nodes:
                coors = self.coors

            else:
                coors = self.coors[:self.n_vertex_dof]

            mesh = Mesh.from_data(self.name, coors[:, :tdim], None, conns,
                                  mat_ids, [descs])

        return mesh

    def get_evaluate_cache(self, cache=None, share_geometry=False,
                           verbose=False):
        """
        Get the evaluate cache for :func:`Variable.evaluate_at()
        <sfepy.discrete.variables.Variable.evaluate_at()>`.

        Parameters
        ----------
        cache : Struct instance, optional
            Optionally, use the provided instance to store the cache data.
        share_geometry : bool
            Set to True to indicate that all the evaluations will work on the
            same region. Certain data are then computed only for the first
            probe and cached.
        verbose : bool
            If False, reduce verbosity.

        Returns
        -------
        cache : Struct instance
            The evaluate cache.
        """
        try:
            from scipy.spatial import cKDTree as KDTree
        except ImportError:
            from scipy.spatial import KDTree

        from sfepy.discrete.fem.geometry_element import \
            create_geometry_elements

        if cache is None:
            cache = Struct(name='evaluate_cache')

        timer = Timer(start=True)
        if (cache.get('cmesh', None) is None) or not share_geometry:
            mesh = self.create_mesh(extra_nodes=False)
            cache.cmesh = cmesh = self.cmesh

            gels = create_geometry_elements()

            cmesh.set_local_entities(gels)
            cmesh.setup_entities()

            cache.centroids = cmesh.get_centroids(cmesh.tdim)

            if self.gel.name != '3_8':
                cache.normals0 = cmesh.get_facet_normals()
                cache.normals1 = None

            else:
                cache.normals0 = cmesh.get_facet_normals(0)
                cache.normals1 = cmesh.get_facet_normals(1)

        output('cmesh setup: %f s' % timer.stop(), verbose=verbose)

        timer.start()
        if (cache.get('kdtree', None) is None) or not share_geometry:
            cache.kdtree = KDTree(cmesh.coors)

        output('kdtree: %f s' % timer.stop(), verbose=verbose)

        return cache

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

        bf = self.eval_basis('v', False, integral)
        bf = bf[:, 0, :].copy()

        data_qp = nm.dot(bf, dofs[self.econn])
        data_qp = nm.swapaxes(data_qp, 0, 1)
        data_qp.shape = data_qp.shape + (1,)

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

    def get_econn(self, conn_type, region, trace_region=None, local=False):
        """
        Get extended connectivity of the given type in the given region.

        Parameters
        ----------
        conn_type: tuple or string
            DOF connectivity type, eg. ('cell', 3) or 'cell'.
            If the topological dimension not specified, it is taken from
            region.tdim.
        region: sfepy.discrete.common.region.Region
            The region for which the connectivity is required.
        trace_region: None or string
            If not None, return mirror connectivity according to `local`.
        local: bool
            If True, return local connectivity w.r.t. facet nodes,
            otherwise return global connectivity w.r.t. all mesh nodes.

        Returns
        -------
        econn: numpy.ndarray
            The extended connectivity array.
        """
        if isinstance(conn_type, tuple):
            integration, tdim = conn_type
        else:
            integration, tdim = conn_type, region.tdim

        if integration == 'cell' and tdim == 1 and self.region.tdim > 1:
            # bar elements
            conn = self.extra_data[f'bars_{region.name}']

        elif (integration in ('cell', 'custom')) and (trace_region is None):
            if region.name == self.region.name:
                conn = self.econn
            else:
                tco = region.kind == 'cell'
                cells = region.get_cells(true_cells_only=tco)
                ii = self.region.get_cell_indices(cells, true_cells_only=tco)
                conn = nm.take(self.econn, ii, axis=0)

        elif integration == 'cell' and trace_region is not None:
            name = f'sd_{region.name}'
            sd = self.extra_data[name]  # FEPhantomSurface
            conn = sd.get_connectivity(local=local, trace_region=trace_region)

        elif integration == 'facet':
            name = f'sd_{region.name}'
            if name not in self.extra_data:
                self.domain.create_surface_group(region)
                self.setup_surface_data(region)

            if self.is_surface:
                local = True
            sd = self.extra_data[name]
            conn = sd.get_connectivity(local=local, trace_region=trace_region)

        elif integration == 'point':
            conn = self.extra_data[f'pd_{region.name}']

        else:
            raise ValueError(f'unknown integration type! ({integration})')

        return conn

    def setup_extra_data(self, info):
        for dct, tdim in set(info.dof_conn_types.values()):
            if dct == 'facet':
                reg = info.get_region()
                mreg_name = info.get_region_name(can_trace=False)
                mreg_name = None if reg.name == mreg_name else mreg_name
                self.domain.create_surface_group(reg)
                self.setup_surface_data(reg, mreg_name)

            elif dct == 'edge':
                raise NotImplementedError('dof connectivity type %s' % dct)

            elif dct == 'point':
                self.setup_point_data(self, info.region)

            elif dct == 'cell' and tdim == 1 and self.region.tdim > 1:
                # bar elements
                self.setup_bar_data(self, info.region)

            elif dct not in ('cell', 'custom'):
                raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def setup_surface_data(self, region, trace_region=None):
        """nodes[leconn] == econn"""
        """nodes are sorted by node number -> same order as region.vertices"""
        name = f'sd_{region.name}'
        if name not in self.extra_data:
            if trace_region is not None and region.tdim == (region.dim - 1):
                sd = FEPhantomSurface(name, region, self.econn)
            else:
                sd = FESurface(name, region, self.efaces, self.econn,
                               self.region)
            self.extra_data[name] = sd

        if name in self.extra_data and trace_region is not None:
            sd = self.extra_data[name]
            sd.setup_mirror_connectivity(region, trace_region)

    def setup_point_data(self, field, region):
        name = f'pd_{region.name}'
        if name not in self.extra_data:
            conn = field.get_dofs_in_region(region, merge=True)
            conn.shape += (1,)
            self.extra_data[name] = conn

    def setup_bar_data(self, field, region):
        name = f'bars_{region.name}'
        if name not in self.extra_data:
            conn = region.domain.get_conn(tdim=1)[region.cells]
            self.extra_data[name] = conn

    def create_mapping(self, region, integral, integration,
                       return_mapping=True):
        """
        Create a new reference mapping.

        Compute jacobians, element volumes and basis function derivatives
        for Volume-type geometries (volume mappings), and jacobians,
        normals and basis function derivatives for Surface-type
        geometries (surface mappings).

        Notes
        -----
        - surface mappings are defined on the surface region
        - surface mappings require field order to be > 0
        """
        domain = self.domain
        coors = domain.get_mesh_coors(actual=True)

        iels = region.get_cells(true_cells_only=(region.kind == 'cell'))
        transform = (self.basis_transform[iels] if self.basis_transform
                     is not None else None)

        geo_ps = self.gel.poly_space
        ps = self.poly_space

        if region.kind == 'cell':
            qp = self.get_qp('v', integral)
            bf = self.eval_basis('v', 0, integral, iels=iels)

            dconn = domain.get_conn(tdim=region.tdim, cells=iels)
            mapping = FEMapping(coors, dconn, poly_space=geo_ps)
            out = mapping.get_mapping(qp.vals, qp.weights, bf, poly_space=ps,
                                      ori=self.ori, transform=transform)

        elif region.kind == 'facet':
            assert_(self.approx_order > 0)

            if self.ori is not None:
                msg = 'surface integrals do not work yet with the' \
                      ' hierarchical basis!'
                raise ValueError(msg)

            if self.basis_transform is not None:
                msg = 'surface integrals do not work with the' \
                      ' basis transform!'
                raise ValueError(msg)

            sd = domain.surface_groups[region.name]
            esd = self.extra_data[f'sd_{region.name}']

            face_indices = region.get_facet_indices()
            cells = face_indices[:, 0]
            dconn = domain.get_conn(tdim=region.tdim, cells=cells)
            conn = sd.get_connectivity()
            mapping = FEMapping(coors, conn, poly_space=geo_ps)

            if not self.is_surface:
                if isinstance(ps.geometry.surface_facet, dict):
                    if integration == 'facet_extra':
                        msg = ('facet integration not supported for '
                               f'element type {self.gel.name}!')
                        raise ValueError(msg)

                    nfc, dim = self.gel.n_face, self.gel.coors.shape[1]
                    bkeys = list(ps.geometry.surface_facet.keys())
                    for bkey in bkeys:
                        self.create_bqp_key(integral, bkey)

                    nqp = nm.max([integral.qps[f'{dim - 1}_{bkey}'].n_point
                                  for bkey in bkeys])

                    flag = ''.join(str(k) for k in bkeys)
                    qp = Struct(
                        name=f'BQP_b{flag}',
                        vals=nm.zeros((nfc, nqp, dim), dtype=nm.float64),
                        weights=nm.zeros((nfc, nqp), dtype=nm.float64)
                    )

                    efc_map = nm.count_nonzero(
                        nm.diff(nm.sort(self.gel.faces)), axis=1) + 1
                    indxs = {}
                    fcaxes = {}
                    ffcidxs = []
                    fc_map = efc_map[sd.fis[:, 1]] * (-1)
                    for ikey, bkey in enumerate(bkeys):
                        idx = efc_map == bkey
                        qp0 = self.qp_coors[(integral.order, f'b{bkey}')]
                        nqp0 = qp0.vals.shape[1]
                        qp.vals[idx, :nqp0, :] = qp0.vals[idx]
                        qp.weights[idx, :nqp0] = qp0.weights
                        if qp0.weights.shape[0] < qp.weights.shape[1]:
                            qp.weights[idx, qp0.weights.shape[0]:] = 0

                        ffcidx = nm.where(efc_map == bkey)[0][0]
                        ffcidxs.append(ffcidx)
                        indx = self.efaces[ffcidx]
                        indx = nm.roll(indx[:bkey], -1)[::-1]
                        indxs[ikey] = indx

                        fcco = ps.geometry.coors[indx]
                        fcax = [nm.where((fcco[1] - fcco[0]) == 1)[0][0],
                                nm.where((fcco[-1] - fcco[0]) == 1)[0][0]]
                        fcaxes[ikey] = fcax

                        fc_map[fc_map == -bkey] = ikey

                    abf = ps.eval_basis(qp.vals, transform=transform)
                    bf = nm.zeros((nfc, nqp, 1, max(bkeys)), dtype=nm.float64)
                    for ifc, efc in enumerate(self.efaces):
                        bkey = efc_map[ifc]
                        bf[ifc, ..., :bkey] = abf[ifc, ...][..., efc[:bkey]]

                    mapping.set_basis_indices(indxs)
                    weights = nm.ascontiguousarray(qp.weights[ffcidxs][fc_map])
                    bf = nm.ascontiguousarray(bf[ffcidxs][fc_map])
                    out = mapping.get_mapping(qp.vals[ffcidxs], weights, bf,
                                              extra=(None, None, None),
                                              is_face=True,
                                              fc_bf_map=(fc_map, fcaxes))
                else:
                    self.create_bqp(region.name, integral)
                    qp = self.qp_coors[(integral.order, esd.bkey)]

                    abf = ps.eval_basis(qp.vals[0], transform=transform)
                    bf = abf[..., self.efaces[0]]

                    indx = self.gel.get_surface_entities()[0]
                    # Fix geometry element's 1st facet orientation for gradients.
                    indx = nm.roll(indx, -1)[::-1]
                    mapping.set_basis_indices(indx)

                    if integration == 'facet_extra':
                        se_bf_bg = geo_ps.eval_basis(qp.vals, diff=True)
                        se_bf_bg = se_bf_bg[sd.fis[:, 1]]
                        se_ebf_bg = self.eval_basis(esd.bkey, 1, integral)
                        se_ebf_bg = se_ebf_bg[sd.fis[:, 1]]
                        remap = prepare_remap(cells, cells.max() + 1)
                        se_conn = dconn[remap[sd.fis[:, 0]], :]
                    else:
                        se_bf_bg, se_ebf_bg, se_conn = None, None, None

                    out = mapping.get_mapping(qp.vals[0], qp.weights, bf,
                        extra=(se_conn, se_bf_bg, se_ebf_bg), is_face=True)

            else:
                # Do not use BQP for surface fields.
                qp = self.get_qp(sd.face_type, integral)
                bf = ps.eval_basis(qp.vals, transform=transform)

                out = mapping.get_mapping(qp.vals, qp.weights, bf,
                                          is_face=True)

        else:
            out = mapping = None

        if out is not None:
            # Store the integral used.
            out.integral = integral
            out.qp = qp
            out.ps = ps

        if return_mapping:
            out = (out, mapping)

        return out

    def average_qp_to_vertices(self, data_qp, integral):
        r"""
        Average data given in quadrature points in region elements into
        region vertices.

        .. math::
           u_n = \sum_e (u_{e,avg} * area_e) / \sum_e area_e
               = \sum_e \int_{area_e} u / \sum area_e
        """
        region = self.region

        n_cells = region.get_n_cells(is_surface=self.is_surface)
        if n_cells != data_qp.shape[0]:
            msg = 'incomatible shape! (%d == %d)' % (n_cells,
                                                     data_qp.shape[0])
            raise ValueError(msg)

        n_vertex = self.n_vertex_dof
        nc = data_qp.shape[2]

        nod_vol = nm.zeros((n_vertex,), dtype=nm.float64)
        data_vertex = nm.zeros((n_vertex, nc), dtype=nm.float64)

        rg = self.get_mapping(self.region, integral, region.kind)[0]

        area = nm.squeeze(rg.volume)
        iels = self.region.get_cells()

        data_e = nm.zeros((area.shape[0], 1, nc, 1), dtype=nm.float64)
        rg.integrate(data_e, data_qp[iels])

        ir = nm.arange(nc, dtype=nm.int32)

        if region.kind == 'cell':
            conn = self.econn[:, :self.gel.n_vertex]
        elif region.kind == 'facet':
            sd = self.domain.surface_groups[region.name]
            # Should be vertex connectivity!
            conn = sd.get_connectivity(local=True)

        for ii, cc in enumerate(conn):
            # Assumes unique nodes in cc!
            ind2, ind1 = nm.meshgrid(ir, cc)
            data_vertex[ind1, ind2] += data_e[iels[ii], 0, :, 0]
            nod_vol[cc] += area[ii]

        data_vertex /= nod_vol[:, nm.newaxis]

        return data_vertex

    def _setup_bubble_dofs(self):
        """
        Setup bubble DOF connectivity for surface field.
        """
        if self.is_surface:
            return 0, None, None


class H1Mixin(Struct):
    """
    Methods of fields specific to H1 space.
    """

    def _setup_shape(self):
        """
        Setup the field's shape-related attributes, see :class:`Field`.
        """
        self.n_components = nm.prod(self.shape)
        self.val_shape = self.shape
