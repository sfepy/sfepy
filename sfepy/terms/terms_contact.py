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
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D'}
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

    def get_fargs(self, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(virtual)

        if self.ci is None:
            self.ci = ContactInfo()

        print self.region.name
        print self.region.shape
        print self.region.facets
        print self.region.get_facet_indices()
        print geo
        print geo.normal

        qps = self.get_physical_qps()
        qp_coors = qps.values
        u_qp = self.get(state, 'val').reshape(qp_coors.shape)

        # Deformed QP coordinates.
        coors = u_qp + qp_coors

        print coors

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return geo, fmode
