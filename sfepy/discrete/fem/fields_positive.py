from sfepy.discrete.fem.fields_base import (VolumeField, SurfaceField,
                                            H1Mixin)
from sfepy.discrete.fem.fields_nodal import GlobalNodalLikeBasis

class H1BernsteinVolumeField(H1Mixin, GlobalNodalLikeBasis, VolumeField):
    family_name = 'volume_H1_bernstein'

    def create_basis_context(self):
        """
        Create the context required for evaluating the field basis.
        """
        # Hack for tests to pass - the reference coordinates are determined
        # from vertices only - we can use the Lagrange basis context for the
        # moment. The true context for Field.evaluate_at() is not implemented.
        gps = self.gel.poly_space
        mesh = self.create_mesh(extra_nodes=False)

        ctx = geo_ctx = gps.create_context(mesh.cmesh, 0, 1e-15, 100, 1e-8)
        ctx.geo_ctx = geo_ctx

        return ctx

class H1BernsteinSurfaceField(H1Mixin, GlobalNodalLikeBasis, SurfaceField):
    family_name = 'surface_H1_bernstein'
