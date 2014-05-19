"""
Reference mappings for isogeometric analysis.
"""
import numpy as nm

from sfepy.discrete.common.mappings import Mapping
import sfepy.discrete.iga.iga as iga
from sfepy.discrete.fem.extmods.mappings import CMapping

class IGMapping(Mapping):
    """
    Reference mapping for isogeometric analysis based on Bezier extraction.
    """

    def __init__(self, domain, cells):
        self.domain = domain
        self.cells = cells
        self.v_shape = (self.domain.shape.n_el, -1, self.domain.shape.dim)
        self.s_shape = (self.domain.shape.n_el, -1, 1)

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
        nurbs = self.domain.nurbs
        variable = nm.ones((nurbs.weights.shape[0], 1), dtype=nm.float64)
        qps, _, _ = iga.eval_variable_in_qp(variable, qp_coors, nurbs.cps,
                                            nurbs.weights, nurbs.degrees,
                                            nurbs.cs, nurbs.conn)
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
        nurbs = self.domain.nurbs
        bfs, bfgs, dets = iga.eval_mapping_data_in_qp(qp_coors, nurbs.cps,
                                                      nurbs.weights,
                                                      nurbs.degrees, nurbs.cs,
                                                      nurbs.conn)
        # Weight Jacobians by quadrature point weights.
        dets = nm.abs(dets) * weights[None, :, None, None]

        # Cell volumes.
        volumes = dets.sum(axis=1)[..., None]

        cmap = CMapping(self.v_shape[0], qp_coors.shape[0], self.v_shape[2],
                        bfs.shape[3], mode='volume', flag=1)

        cmap.bf[:] = bfs
        cmap.bfg[:] = bfgs
        cmap.det[:] = dets
        cmap.volume[:] = volumes

        return cmap
