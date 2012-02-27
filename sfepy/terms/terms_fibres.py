import numpy as nm

from sfepy.terms.terms_hyperelastic_tl import HyperElasticTLBase
from sfepy.homogenization.utils import iter_sym

def fibre_function(out, pars, green_strain, fmode):
    """
    Depending on `fmode`, compute fibre stress (0) or tangent modulus (!= 0).
    """
    fmax, eps_opt, s, fdir, act = pars

    eps = nm.zeros_like(fmax)
    omega = nm.empty_like(green_strain)
    for ii, (ir, ic) in enumerate(iter_sym(fdir.shape[2])):
        omega[..., ii, 0] = fdir[..., ir, 0] * fdir[..., ic, 0]
        eps[..., 0, 0] += omega[..., ii, 0] * green_strain[..., ii, 0]

    tau = act * fmax * nm.exp(-((eps - eps_opt) / s)**2.0)

    if fmode == 0:
        out[:] = omega * tau

    else:
        for ir in range(omega.shape[2]):
            for ic in range(omega.shape[2]):
                out[..., ir, ic] = omega[..., ir, 0] * omega[..., ic, 0]

        out[:] *= -2.0 * ((eps - eps_opt) / (s**2.0)) * tau

    return out

class FibresActiveTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic active fibres term. Effective stress
    :math:`S_{ij} = A f_{\rm max} \exp{\left\{-(\frac{\epsilon -
    \varepsilon_{\rm opt}}{s})^2\right\}} d_i d_j`,
    where :math:`\epsilon = E_{ij} d_i d_j` is the Green strain
    :math:`\ull{E}` projected to the fibre direction :math:`\ul{d}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material_1 : :math:`f_{\rm max}`
        - material_2 : :math:`\varepsilon_{\rm opt}`
        - material_3 : :math:`s`
        - material_4 : :math:`\ul{d}`
        - material_5 : :math:`A`
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_tl_fib_a'
    arg_types = ('material_1', 'material_2', 'material_3',
                 'material_4', 'material_5', 'virtual', 'state')
    family_data_names = ['green_strain']

    def get_fargs(self, mat1, mat2, mat3, mat4, mat5, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        fargs = HyperElasticTLBase.get_fargs(self,
                                             (mat1, mat2, mat3, mat4, mat5),
                                             virtual, state,
                                             mode, term_mode, diff_var,
                                             **kwargs)
        return fargs

    @staticmethod
    def stress_function(out, pars, green_strain):
        fibre_function(out, pars, green_strain, 0)

    @staticmethod
    def tan_mod_function(out, pars, green_strain):
        fibre_function(out, pars, green_strain, 1)
