import numpy as nm

from sfepy.base.base import Struct, assert_
from sfepy.discrete.fem.mappings import VolumeMapping, SurfaceMapping
from poly_spaces import PolySpace
from fe_surface import FESurface

class SurfaceInterpolant(Interpolant):
    """
    Like Interpolant, but for use with SurfaceField and
    SurfaceApproximation.
    """

    def __init__(self, name, gel, space='H1', base='lagrange',
                 approx_order=1, force_bubble=False):
        Interpolant.__init__(self, name, gel, space=space, base=base,
                             approx_order=approx_order,
                             force_bubble=force_bubble)

        # Make alias 'v' <-> 's#'.
        ps = self.poly_spaces['v']
        self.poly_spaces['s%d' % ps.n_nod] = ps

    def get_geom_poly_space(self, key):
        assert_(key[0] == 's')

        ps = self.gel.poly_space

        return ps

class Approximation(Struct):

    def __init__(self, name, interp, region, is_surface=False):
        """interp, region are borrowed."""

class DiscontinuousApproximation(Approximation):

    def eval_extra_coor(self, coors, mesh_coors):
        """
        Compute coordinates of extra nodes. For discontinuous
        approximations, all nodes are treated as extra.
        """
        gps = self.interp.gel.poly_space
        ps = self.interp.poly_spaces['v']

        eval_nodal_coors(coors, mesh_coors, self.region, ps, gps,
                         self.econn, only_extra=False)

class SurfaceApproximation(Approximation):

    def __init__(self, name, interp, region):
        Approximation.__init__(self, name, interp, region, is_surface=True)

    def get_qp(self, key, integral):
        """
        Get quadrature points and weights corresponding to the given key
        and integral. The key is 's#', where # is the number of
        face vertices.
        """
        assert_(key[0] == 's')
        qpkey = (integral.order, key)

        if not self.qp_coors.has_key(qpkey):
            interp = self.interp
            geometry = interp.gel.name

            vals, weights = integral.get_qp(geometry)
            self.qp_coors[qpkey] = Struct(vals=vals, weights=weights)

        return self.qp_coors[qpkey]
