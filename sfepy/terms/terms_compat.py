from sfepy.terms.terms_dot import DotProductTerm
from sfepy.terms.terms_basic import IntegrateTerm, IntegrateOperatorTerm,\
    VolumeTerm, IntegrateMatTerm
from sfepy.terms.terms_elastic import CauchyStrainTerm
from sfepy.terms.terms_navier_stokes import GradTerm, DivTerm


# deprecated names, will be removed in future

class DotVolumeProductTerm(DotProductTerm):
    name = 'dw_volume_dot'


class DotSurfaceProductTerm(DotProductTerm):
    name = 'dw_surface_dot'


class IntegrateVolumeTerm(IntegrateTerm):
    name = 'ev_volume_integrate'


class IntegrateSurfaceTerm(IntegrateTerm):
    name = 'ev_surface_integrate'


class IntegrateVolumeOperatorTerm(IntegrateOperatorTerm):
    name = 'dw_volume_integrate'


class IntegrateSurfaceOperatorTerm(IntegrateOperatorTerm):
    name = 'dw_surface_integrate'


class VolumeXTerm(VolumeTerm):
    name = 'd_volume'


class SurfaceTerm(VolumeTerm):
    name = 'd_surface'


class IntegrateVolumeMatTerm(IntegrateMatTerm):
    name = 'ev_volume_integrate_mat'


class IntegrateSurfaceMatTerm(IntegrateMatTerm):
    name = 'ev_surface_integrate_mat'


class CauchyStrainSTerm(CauchyStrainTerm):
    name = 'ev_cauchy_strain_s'


class SurfaceGradTerm(GradTerm):
    name = 'ev_surface_grad'


class SurfaceDivTerm(DivTerm):
    name = 'ev_surface_div'