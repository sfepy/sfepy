"""
Fields for Discontinous Galerkin
"""
import numpy as nm

# sfepy imports
from sfepy.base.base import assert_, basestr, Struct
from sfepy.discrete.common.fields import parse_shape, Field
from sfepy.discrete import Integral
from sfepy.discrete.iga.mappings import IGMapping
from sfepy.discrete.iga.iga import get_bezier_element_entities
from six.moves import range

# local imports
from dg_basis import LegendrePolySpace

class DGField(Field):
    family_name = 'volume_DG'

    def __init__(self, name, dtype, shape, region, space="H1",
                 poly_space_base=None,  approx_order=None):
        # TODO what else will be needed in field fo DG?

        if poly_space_base is None:
            # TODO where does polyspace object come from in other filed implementations
            # TODO how to get geometry for it
            poly_space_base = LegendrePolySpace("legendre", region.geometry, approx_order)

        self.integrate = lambda fun: Integral("dg_fi", order=approx_order).integrate(fun, order=approx_order)

        Struct.__init__(self, name=name, dtype=dtype, shape=shape,
                        region=region, space=space,
                        poly_space_base=poly_space_base,
                        approx_order=approx_order)

    def create_mapping(self, region, integral, integration):
        # TODO create a ne reference mapping, maybe steal this from FE
        raise NotImplemented

    def set_dofs(self, fun=0.0, region=None, dpn=None, warn=None):
        # TODO used setup DOF directly, basically sampleIC from TSs

        mesh = self.region.mesh
        sic = nm.zeros((2, self.mesh.n_el, 1), dtype=nm.float64)

        c = (mesh.coors[1:] + mesh.coors[:-1]) / 2  # center
        s = (mesh.coors[1:] - mesh.coors[:-1]) / 2  # scale
        sic[0, :] = self.integrate(lambda t: fun(c + t * s)) / 2
        sic[1, :] = 3 * self.integrate(lambda t: t * fun(c + t * s)) / 2

        vals = sic
        nods = c

        return nods, vals

# _get_facets
# create_basis_context
# create_eval_mesh
# create_mesh
# create_output
# get_data_shape
# get_dofs_in_region
# get_econn
# get_true_order
# is_higher_order
# setup_extra_data

