from sfepy.terms.terms import *
from sfepy.terms.terms_base import VectorVector
from sfepy.terms.terms_hyperelastic_tl import HyperElasticTLBase
from sfepy.homogenization.utils import iter_sym

class FibresActiveTLTerm(VectorVector, HyperElasticTLBase):
    r"""
    :Description:
    Hyperelastic active fibres term. Effective stress
    :math:`S_{ij} = A f_{\rm max} \exp{\left\{-(\frac{\epsilon -
    \varepsilon_{\rm opt}}{s})^2\right\}} d_i d_j`,
    where :math:`\epsilon = E_{ij} d_i d_j` is the Green strain
    :math:`\ull{E}` projected to the fibre direction :math:`\ul{d}`.

    :Definition:
    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        material_1 : :math:`f_{\rm max}`,
        material_2 : :math:`\varepsilon_{\rm opt}`,
        material_3 : :math:`s`,
        material_4 : :math:`\ul{d}`,
        material_5 : :math:`A`,
        virtual    : :math:`\ul{v}`,
        state      : :math:`\ul{u}`

    """
    name = 'dw_tl_fib_a'
    arg_types = ('material_1', 'material_2', 'material_3',
                 'material_4', 'material_5', 'virtual', 'state')
    geometry = [(Volume, 'virtual')]
    family_data_names = ['E']
    
    def compute_crt_data( self, family_data, mode, **kwargs ):
        pars = self.get_args(['material_1', 'material_2', 'material_3',
                              'material_4', 'material_5'], **kwargs)
        fmax, eps_opt, s, fdir, act = pars
        strainE = family_data[0]

        eps = nm.zeros_like(fmax)
        omega = nm.empty_like(strainE)
        for ii, (ir, ic) in enumerate(iter_sym(fdir.shape[2])):
            omega[...,ii,0] = fdir[...,ir,0] * fdir[...,ic,0]
            eps[...,0,0] += omega[...,ii,0] * strainE[...,ii,0]

        tau = act * fmax * nm.exp(-((eps - eps_opt) / s)**2.0)
        
        if mode == 0:
            out = omega * tau

        else:
            shape = list(strainE.shape)
            shape[-1] = shape[-2]
            out = nm.empty(shape, dtype=nm.float64)

            for ir in range(omega.shape[2]):
                for ic in range(omega.shape[2]):
                    out[...,ir,ic] = omega[...,ir,0] * omega[...,ic,0]

            out[:] *= -2.0 * ((eps - eps_opt) / (s**2.0)) * tau

        return out
 
