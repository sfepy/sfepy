"""
Finite element reference mappings for structural elements.
"""
import numpy as nm

from sfepy.linalg import dot_sequences as ddot
from sfepy.discrete.common.mappings import Mapping, PhysicalQPs
import sfepy.mechanics.shell10x as shell10x

class Shell10XMapping(Mapping):
    """
    The reference mapping for the shell10x element.
    """

    def __init__(self, region, field):
        self.region = region
        self.field = field

    def get_physical_qps(self, qp_coors):
        """
        Get physical quadrature points corresponding the given reference
        element quadrature points.

        Returns
        -------
        qps : array
            The physical quadrature points ordered element by element,
            i.e. with shape (n_el, n_qp, dim).
        """
        phys_qps = PhysicalQPs()

        bf = self.bfu.squeeze()

        qps_loc = nm.einsum('qi,cqij->cqj', bf, self.coors_loc_3d)

        qps = nm.einsum('cqi,cji->cqj', qps_loc, self.mtx_t)
        qps += self.e_centres[:, None, :]

        n_el, n_qp = qps.shape[0], qps.shape[1]

        phys_qps.num = n_el * n_qp
        phys_qps.shape = qps.shape

        qps.shape = (phys_qps.num, qps.shape[2])
        phys_qps.values = qps

        return phys_qps

    def get_mapping(self, qp_coors, weights):
        """
        Get the mapping for given quadrature points and weights.
        """
        domain = self.region.domain
        mesh = domain.mesh

        iels = self.region.get_cells()
        conn = nm.take(domain.get_conn(tdim=self.region.tdim),
                       iels.astype(nm.int32), axis=0)

        e_coors = mesh.coors[conn]
        mtx_t = shell10x.create_transformation_matrix(e_coors)

        e_centres = self.region.cmesh.get_centroids(2)[iels]

        coors_loc = ddot((e_coors - e_centres[:, None, :]), mtx_t)

        ebs = shell10x.create_local_bases(coors_loc)
        rops = shell10x.create_rotation_ops(ebs)

        ps = self.field.poly_space

        qp_coors0 = nm.array([[0.5, 0.5, 0.5]])
        # Should be thickness, but not used anywhere.
        qp_weights0 = nm.array([0.0])

        dxidx0, det0 = shell10x.get_mapping_data(ebs, rops, ps, coors_loc,
                                                 qp_coors0, qp_weights0,
                                                 special_dx3=True)

        aux = shell10x.get_mapping_data(ebs, rops, ps, coors_loc,
                                        qp_coors, weights)
        coors_loc_3d, bfu, bfgm, dxidx, det = aux

        self.coors_loc = coors_loc
        self.e_centres = e_centres
        self.mtx_t = mtx_t
        self.ebs = ebs
        self.rops = rops

        self.dxidx0 = dxidx0
        self.det0 = det0

        self.coors_loc_3d = coors_loc_3d
        self.bfu = bfu
        self.bfgm = bfgm
        self.dxidx = dxidx
        self.det = det
        self.volume = (det * weights).sum(1)

        return self
