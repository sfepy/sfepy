"""
Reference mappings for isogeometric analysis.
"""
import numpy as nm

from sfepy.discrete.common.mappings import Mapping, PyCMapping
import sfepy.discrete.iga.extmods.igac as iga


class IGMapping(Mapping):
    """
    Reference mapping for isogeometric analysis based on Bezier extraction.

    Parameters
    ----------
    domain : IGDomain instance
        The mapping domain.
    cells : array
        The mapping region cells. (All domain cells required.)
    nurbs : NurbsPatch instance, optional
        If given, the `nurbs` is used instead of `domain.nurbs`. The `nurbs`
        has to be obtained by degree elevation of `domain.nurbs`.
    """

    def __init__(self, domain, cells, nurbs=None):
        self.domain = domain
        self.cells = cells
        self.nurbs = domain.nurbs if nurbs is None else nurbs
        self.v_shape = (len(cells), -1, self.domain.shape.dim)
        self.s_shape = (len(cells), -1, 1)

    def get_geometry(self):
        """
        Return reference element geometry as a GeometryElement instance.
        """
        return self.domain.gel

    def get_physical_qps(self, qp_coors):
        """
        Get physical quadrature points corresponding to given reference
        Bezier element quadrature points.

        Returns
        -------
        qps : array
            The physical quadrature points ordered element by element,
            i.e. with shape (n_el, n_qp, dim).
        """
        nurbs = self.nurbs
        variable = nm.ones((nurbs.weights.shape[0], 1), dtype=nm.float64)
        qps, _, _ = iga.eval_variable_in_qp(variable, qp_coors, nurbs.cps,
                                            nurbs.weights, nurbs.degrees,
                                            nurbs.cs, nurbs.conn, self.cells)
        qps = qps.reshape(self.v_shape)

        return qps

    def get_mapping(self, qp_coors, weights):
        """
        Get the mapping for given quadrature points and weights.

        Returns
        -------
        cmap : CMapping instance
            The reference mapping.

        Notes
        -----
        Does not set total volume of the C mapping structure!
        """
        nurbs = self.nurbs
        bfs, bfgs, dets = iga.eval_mapping_data_in_qp(qp_coors, nurbs.cps,
                                                      nurbs.weights,
                                                      nurbs.degrees, nurbs.cs,
                                                      nurbs.conn, self.cells)
        # Weight Jacobians by quadrature point weights.
        dets = nm.abs(dets) * weights[None, :, None, None]

        # Cell volumes.
        volumes = dets.sum(axis=1)[..., None]

        pycmap = PyCMapping(bfs, dets, volumes, bfgs, None, self.v_shape[2])

        return pycmap
