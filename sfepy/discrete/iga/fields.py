"""
Fields for isogeometric analysis.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import assert_, basestr, Struct
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete.iga.mappings import IGMapping
from sfepy.discrete.iga.iga import get_bezier_element_entities
from six.moves import range

def parse_approx_order(approx_order):
    if (approx_order is None): return 0

    aux = approx_order.split('+')
    if len(aux) == 2:
        return int(aux[1])

    else:
        return 0

class IGField(Field):
    """
    Bezier extraction based NURBS field for isogeometric analysis.

    Notes
    -----
    The field has to cover the whole IGA domain. The field's NURBS basis can
    have higher degree than the domain NURBS basis.
    """
    family_name = 'volume_H1_iga'

    def __init__(self, name, dtype, shape, region, approx_order=None,
                 **kwargs):
        """
        Create a Bezier element isogeometric analysis field.

        Parameters
        ----------
        name : str
            The field name.
        dtype : numpy.dtype
            The field data type: float64 or complex128.
        shape : int/tuple/str
            The field shape: 1 or (1,) or 'scalar', space dimension (2, or (2,)
            or 3 or (3,)) or 'vector', or a tuple. The field shape determines
            the shape of the FE base functions and is related to the number of
            components of variables and to the DOF per node count, depending
            on the field kind.
        region : Region
            The region where the field is defined.
        approx_order : str or tuple, optional
            The field approximation order string or tuple with the first
            component in the form 'iga+<nonnegative int>'. Other components are
            ignored. The nonnegative int corresponds to the number of times the
            degree is elevated by one w.r.t. the domain NURBS description.
        **kwargs : dict
            Additional keyword arguments.
        """
        shape = parse_shape(shape, region.domain.shape.dim)

        if approx_order is None:
            elevate_times = 0

        else:
            if isinstance(approx_order, basestr): approx_order = (approx_order,)
            elevate_times = parse_approx_order(approx_order[0])

        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region, elevate_times=elevate_times)

        self.domain = self.region.domain
        self.nurbs = self.domain.nurbs.elevate(elevate_times)

        self._setup_kind()

        self.n_components = nm.prod(self.shape)
        self.val_shape = self.shape
        self.n_nod = self.nurbs.weights.shape[0]
        self.n_efun = nm.prod(self.nurbs.degrees + 1)
        self.approx_order = self.nurbs.degrees.max()

        self.mappings = {}

        self.is_surface = False

    def _get_facets(self, tdim):
        aux = get_bezier_element_entities(self.nurbs.degrees)
        return aux[2 - tdim]

    def get_true_order(self):
        return nm.prod(self.nurbs.degrees)

    def is_higher_order(self):
        """
        Return True, if the field's approximation order is greater than one.
        """
        return (self.nurbs.degrees > 1).any()

    def get_econn(self, conn_type, region, is_trace=False, integration=None):
        """
        Get DOF connectivity of the given type in the given region.
        """
        ct = conn_type.type if isinstance(conn_type, Struct) else conn_type

        nconn = self.nurbs.conn

        if ct == 'volume':
            if region.name == self.region.name:
                conn = nconn

            else:
                cells = region.get_cells(true_cells_only=True)
                conn = nm.take(nconn, cells.astype(nm.int32), axis=0)

        elif ct == 'surface':
            fis = region.get_facet_indices()
            tdim = region.kind_tdim
            facets = self._get_facets(tdim)

            conn = []
            for ii, fi in enumerate(fis):
                conn.append(nconn[fi[0], facets[fi[1]]])

        else:
            raise ValueError('unsupported connectivity type! (%s)' % ct)

        return conn

    def get_data_shape(self, integral, integration='volume', region_name=None):
        """
        Get element data dimensions.

        Parameters
        ----------
        integral : Integral instance
            The integral describing used numerical quadrature.
        integration : 'volume'
            The term integration type. Only 'volume' type is implemented.
        region_name : str
            The name of the region of the integral.

        Returns
        -------
        data_shape : 4 ints
            The `(n_el, n_qp, dim, n_en)` for volume shape kind.

        Notes
        -----
        - `n_el` = number of elements
        - `n_qp` = number of quadrature points per element/facet
        - `dim` = spatial dimension
        - `n_en` = number of element nodes
        """
        region = self.domain.regions[region_name]

        _, weights = integral.get_qp(self.domain.gel.name)
        n_qp = weights.shape[0]

        data_shape = (region.shape.n_cell, n_qp, region.dim, self.n_efun)

        return data_shape

    def get_dofs_in_region(self, region, merge=True):
        """
        Return indices of DOFs that belong to the given region and group.

        Notes
        -----
        `merge` is not used.
        """
        idim = region.kind_tdim
        if idim < (self.domain.shape.tdim - 1):
            raise ValueError('region "%s" has no facets or cells!'
                             % region.name)

        if idim == region.tdim: # Cells.
            rconn = self.get_econn('volume', region)
            dofs = nm.unique(rconn)

        else: # Facets.
            rconn = self.get_econn('surface', region)
            dofs = nm.unique(nm.concatenate(rconn))

        if not merge:
            dofs = [dofs]

        return dofs

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        """
        Set the values of DOFs given by the `region` using a function of space
        coordinates or value `fun`.

        If `fun` is a function, the l2 projection that is global for all region
        facets is used to set the DOFs.

        If `dpn > 1`, and `fun` is a function, it has to return the values
        DOF-by-DOF, i.e. a single one-dimensional vector with all values of the
        first component, then of the second one etc. concatenated
        together.

        Parameters
        ----------
        fun : float or array of length dpn or callable
            The DOF values.
        region : Region
            The region containing the DOFs.
        dpn : int, optional
            The DOF-per-node count. If not given, the number of field
            components is used.
        warn : str, optional
            The warning message printed when the region selects no DOFs.

        Returns
        -------
        nods : array, shape (n_dof,)
            The field DOFs (or node indices) given by the region.
        vals : array, shape (dpn, n_dof)
            The values of the DOFs, DOF-by-DOF when raveled in C (row-major)
            order.
        """
        if region is None:
            region = self.region

        if dpn is None:
            dpn = self.n_components

        nods = []
        vals = []

        aux = self.get_dofs_in_region(region)
        nods = nm.unique(nm.hstack(aux))

        if nm.isscalar(fun):
            vals = nm.repeat([fun], nods.shape[0] * dpn)

        elif isinstance(fun, nm.ndarray):
            assert_(len(fun) == dpn)
            vals = nm.repeat(fun, nods.shape[0])

        elif callable(fun):
            import scipy.sparse as sps
            from sfepy.solvers.ls import solve
            from sfepy.discrete.integrals import Integral
            from sfepy.discrete.fem.utils import prepare_remap
            import sfepy.discrete.iga as iga
            from sfepy.discrete.iga.extmods.igac import eval_mapping_data_in_qp

            nurbs = self.nurbs
            facets = self._get_facets(region.kind_tdim)

            # Region facet connectivity.
            rconn = self.get_econn('surface', region)

            # Local connectivity.
            remap = prepare_remap(nods, nods.max() + 1)
            lconn = [remap[ii] for ii in rconn]

            # Cell and face(cell) ids for each facet.
            fis = region.get_facet_indices()

            # Integral given by max. NURBS surface degree.
            fdegrees = iga.get_surface_degrees(nurbs.degrees)
            order = fdegrees.max()
            integral = Integral('i', order=2*order)
            vals, weights = integral.get_qp(self.domain.gel.surface_facet_name)

            # Boundary QP - use tensor product structure.
            bvals = iga.create_boundary_qp(vals, region.tdim)

            # Compute facet basis, jacobians and physical BQP.
            n_dof = len(nods)
            rhs = nm.zeros((dpn, n_dof), dtype=nm.float64)
            rows, cols, mvals = [], [], []
            all_qp = []
            all_fbfs = []
            all_dets = []
            for ii, (ie, ifa) in enumerate(fis):
                qp_coors = bvals[ifa]

                bfs, _, dets = eval_mapping_data_in_qp(qp_coors, nurbs.cps,
                                                       nurbs.weights,
                                                       nurbs.degrees,
                                                       nurbs.cs,
                                                       nurbs.conn,
                                                       nm.array([ie]))
                # Facet basis.
                fbfs = bfs[..., facets[ifa]][0, :, 0, :]

                # Weight Jacobians by quadrature point weights.
                dets = nm.abs(dets) * weights[None, :, None, None]
                dets = dets[0, :, 0, :]

                # Physical BQP.
                fcps = nurbs.cps[nurbs.conn[ie, facets[ifa]]]
                qp = nm.dot(fbfs, fcps)

                all_qp.append(qp)
                all_fbfs.append(fbfs)
                all_dets.append(dets)

            # DOF values in the physical BQP.
            qps = nm.concatenate(all_qp)
            vals = nm.asarray(fun(qps))
            vals.shape = (dpn, qps.shape[0])

            n_qp_face = len(bvals[0])

            # Assemble l2 projection system.
            for ii, (ie, ifa) in enumerate(fis):
                # Assembling indices.
                elc = lconn[ii]

                fvals = vals[:, n_qp_face * ii : n_qp_face * (ii + 1)]

                fbfs = all_fbfs[ii]
                dets = all_dets[ii]

                # Local projection system.
                for idof in range(dpn):
                    lrhs = (fbfs * (fvals[idof, :, None] * dets)).sum(0)
                    rhs[idof, elc] += lrhs

                lmtx = ((fbfs[..., None] * fbfs[:, None, :])
                        * dets[..., None]).sum(0)

                er, ec = nm.meshgrid(elc, elc)
                rows.append(er.ravel())
                cols.append(ec.ravel())
                mvals.append(lmtx.ravel())

            rows = nm.concatenate(rows)
            cols = nm.concatenate(cols)
            mvals = nm.concatenate(mvals)
            mtx = sps.coo_matrix((mvals, (rows, cols)), shape=(n_dof, n_dof))

            vals = nm.zeros((n_dof, dpn), dtype=nm.float64)

            # Solve l2 projection system.
            for idof in range(dpn):
                dofs = solve(mtx, rhs[idof, :])
                vals[remap[nods], idof] = dofs

        else:
            raise ValueError('unknown function/value type! (%s)' % type(fun))

        vals.shape = (len(nods), -1)

        return nods, vals

    def setup_extra_data(self, geometry, info, is_trace):
        dct = info.dc_type.type

        if dct != 'volume':
            raise ValueError('unknown dof connectivity type! (%s)' % dct)

    def create_mapping(self, region, integral, integration):
        """
        Create a new reference mapping.
        """
        vals, weights = integral.get_qp(self.domain.gel.name)
        mapping = IGMapping(self.domain, region.cells, nurbs=self.nurbs)
        cmap = mapping.get_mapping(vals, weights)

        return cmap, mapping

    def create_mesh(self, extra_nodes=True):
        """
        Create a mesh corresponding to the field region. For IGA fields, this
        is directly the topological mesh. The `extra_nodes` argument is
        ignored.
        """
        return self.domain.mesh

    def create_eval_mesh(self):
        """
        Create a mesh with the original NURBS connectivity for evaluating the
        field. The mesh coordinates are the NURBS control points.
        """
        return self.domain.eval_mesh

    def create_basis_context(self):
        """
        Create the context required for evaluating the field basis.
        """
        from sfepy.discrete.iga.extmods.igac import CNURBSContext

        nurbs = self.nurbs

        ctx = CNURBSContext(nurbs.cps, nurbs.weights, nurbs.degrees,
                            nurbs.cs, nurbs.conn, i_max=100)

        return ctx

    def create_output(self, dofs, var_name, dof_names=None,
                      key=None, **kwargs):
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

        Returns
        -------
        out : dict
            The output dictionary.
        """
        from sfepy.discrete.iga.utils import create_mesh_and_output

        key = key if key is not None else var_name

        num = 25 if self.region.dim == 3 else 101
        pars = (nm.linspace(0, 1, num),) * self.region.dim
        mesh, out = create_mesh_and_output(self.nurbs, pars, **{key : dofs})
        out[key].var_name = var_name
        out[key].dofs = dof_names
        out['__mesh__'] = mesh

        return out
