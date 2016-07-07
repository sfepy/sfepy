import numpy as nm

from sfepy.base.base import Struct
from sfepy.terms.terms import Term

class ContactInfo(Struct):
    """
    Various contact-related data of contact terms.
    """
    pass

class ContactTerm(Term):
    r"""
    """
    name = 'dw_contact'
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : 'str', 'virtual' : ('D', 'state'), 'state' : 'D'}
    integration = 'surface'

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.ci = None

    @staticmethod
    def function(out, geo, fmode):
        out_qp = nm.zeros((out.shape[0], geo.n_qp) + out.shape[2:],
                          dtype=out.dtype)
        status = geo.integrate(out, nm.ascontiguousarray(out_qp))

        return status

    def get_fargs(self, material, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(virtual)

        region0 = self.region
        region1 = self.region.domain.regions[material]

        if self.ci is None:
            self.ci = ContactInfo()

        print(region0.name)
        print(region0.shape)
        print(region0.facets)
        print(region0.get_facet_indices())
        print(state.field.surface_data[region0.name].fis)

        print(region1.name)
        print(region1.shape)
        print(region1.facets)
        print(region1.get_facet_indices())
        print(state.field.surface_data[region1.name].fis)

        print(geo)
        print(geo.normal)

        mesh_coors = self.region.domain.mesh.coors
        print(mesh_coors[region0.vertices])

        # Uses field connectivity (higher order nodes).
        sd0 = state.field.surface_data[region0.name]
        sd1 = state.field.surface_data[region1.name]

        # Uses mesh connectivity.
        sdg0 = self.region.domain.surface_groups[region0.name]
        sdg1 = self.region.domain.surface_groups[region1.name]

        print(mesh_coors[sdg0.econn])
        print(mesh_coors[sdg1.econn])

        qps = self.get_physical_qps()
        qp_coors = qps.values
        u_qp = self.get(state, 'val').reshape(qp_coors.shape)

        # Deformed QP coordinates.
        coors = u_qp + qp_coors

        print(coors)

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return geo, fmode
