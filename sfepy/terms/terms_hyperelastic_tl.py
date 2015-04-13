import numpy as nm

from sfepy.base.base import assert_, Struct
from sfepy.terms.terms import terms
from sfepy.terms.terms_hyperelastic_base import HyperElasticBase

class HyperElasticTLBase(HyperElasticBase):
    """
    Base class for all hyperelastic terms in TL formulation family.

    The subclasses should have the following static method attributes:
    - `stress_function()` (the stress)
    - `tan_mod_function()` (the tangent modulus)

    The common (family) data are cached in the evaluate cache of state
    variable.
    """
    family_function = staticmethod(terms.dq_finite_strain_tl)
    weak_function = staticmethod(terms.dw_he_rtm)
    fd_cache_name = 'tl_common'
    hyperelastic_mode = 0

    def compute_family_data(self, state):
        ap, vg = self.get_approximation(state)

        vec = self.get_vector(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = dim * (dim + 1) / 2

        shapes = {
            'mtx_f' : (n_el, n_qp, dim, dim),
            'det_f' : (n_el, n_qp, 1, 1),
            'sym_c' : (n_el, n_qp, sym, 1),
            'tr_c' : (n_el, n_qp, 1, 1),
            'in2_c' : (n_el, n_qp, 1, 1),
            'sym_inv_c' : (n_el, n_qp, sym, 1),
            'green_strain' : (n_el, n_qp, sym, 1),
        }
        data = Struct(name='tl_family_data')
        for key, shape in shapes.iteritems():
            setattr(data, key, nm.zeros(shape, dtype=nm.float64))

        self.family_function(data.mtx_f,
                             data.det_f,
                             data.sym_c,
                             data.tr_c,
                             data.in2_c,
                             data.sym_inv_c,
                             data.green_strain,
                             vec, vg, ap.econn)
        return data

class NeoHookeanTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic neo-Hookean term. Effective stress
    :math:`S_{ij} = \mu J^{-\frac{2}{3}}(\delta_{ij} -
    \frac{1}{3}C_{kk}C_{ij}^{-1})`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`\mu`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_tl_he_neohook'
    family_data_names = ['det_f', 'tr_c', 'sym_inv_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_neohook)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_neohook)

class MooneyRivlinTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic Mooney-Rivlin term. Effective stress
    :math:`S_{ij} = \kappa J^{-\frac{4}{3}} (C_{kk} \delta_{ij} - C_{ij}
    - \frac{2}{3 } I_2 C_{ij}^{-1})`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`\kappa`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_tl_he_mooney_rivlin'
    family_data_names = ['det_f', 'tr_c', 'sym_inv_c', 'sym_c', 'in2_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_mooney_rivlin)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_mooney_rivlin)

class BulkPenaltyTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic bulk penalty term. Stress
    :math:`S_{ij} = K(J-1)\; J C_{ij}^{-1}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`K`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """

    name = 'dw_tl_bulk_penalty'
    family_data_names = ['det_f', 'sym_inv_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_bulk)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_bulk)

class BulkActiveTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic bulk active term. Stress :math:`S_{ij} = A J C_{ij}^{-1}`,
    where :math:`A` is the activation in :math:`[0, F_{\rm max}]`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(\ul{u}) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - material : :math:`A`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """

    name = 'dw_tl_bulk_active'
    family_data_names = ['det_f', 'sym_inv_c']

    stress_function = staticmethod(terms.dq_tl_he_stress_bulk_active)
    tan_mod_function = staticmethod(terms.dq_tl_he_tan_mod_bulk_active)

class BulkPressureTLTerm(HyperElasticTLBase):
    r"""
    Hyperelastic bulk pressure term. Stress
    :math:`S_{ij} = -p J C_{ij}^{-1}`.

    :Definition:

    .. math::
        \int_{\Omega} S_{ij}(p) \delta E_{ij}(\ul{u};\ul{v})

    :Arguments:
        - virtual : :math:`\ul{v}`
        - state   : :math:`\ul{u}`
        - state_p : :math:`p`
    """
    name = 'dw_tl_bulk_pressure'
    arg_types = ('virtual', 'state', 'state_p')
    arg_shapes = {'virtual' : ('D', 'state'), 'state' : 'D', 'state_p' : 1}
    family_data_names = ['det_f', 'sym_inv_c']

    family_function = staticmethod(terms.dq_finite_strain_tl)
    weak_function = staticmethod(terms.dw_he_rtm)
    weak_dp_function = staticmethod(terms.dw_tl_volume)

    stress_function = staticmethod(terms.dq_tl_stress_bulk_pressure)
    tan_mod_u_function = staticmethod(terms.dq_tl_tan_mod_bulk_pressure_u)

    def compute_data(self, family_data, mode, **kwargs):
        det_f, sym_inv_c = family_data.det_f, family_data.sym_inv_c
        p_qp = family_data.p_qp

        if mode == 0:
            out = nm.empty_like(sym_inv_c)
            fun = self.stress_function

        elif mode == 1:
            shape = list(sym_inv_c.shape)
            shape[-1] = shape[-2]
            out = nm.empty(shape, dtype=nm.float64)
            fun = self.tan_mod_u_function

        else:
            raise ValueError('bad mode! (%d)' % mode)

        fun(out, p_qp, det_f, sym_inv_c)

        return out

    def get_fargs(self, virtual, state, state_p,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgv, _ = self.get_mapping(state)

        fd = self.get_family_data(state, 'tl_common', self.family_data_names)
        fd.p_qp = self.get(state_p, 'val')

        if mode == 'weak':
            if diff_var != state_p.name:
                if diff_var is None:
                    stress = self.compute_data(fd, 0, **kwargs)
                    self.stress_cache = stress
                    tan_mod = nm.array([0], ndmin=4, dtype=nm.float64)

                    fmode = 0

                else:
                    stress = self.stress_cache
                    if stress is None:
                        stress = self.compute_data(fd, 0, **kwargs)

                    tan_mod = self.compute_data(fd, 1, **kwargs)
                    fmode = 1

                fargs = (self.weak_function,
                         stress, tan_mod, fd.mtx_f, fd.det_f, vgv, fmode, 0)

            else:
                vgs, _ = self.get_mapping(state_p)

                fargs =  (self.weak_dp_function,
                          fd.mtx_f, fd.sym_inv_c, fd.det_f, vgs, vgv, 1, -1)

            return fargs

        elif mode == 'el_avg':
            if term_mode == 'strain':
                out_qp = fd.green_strain

            elif term_mode == 'stress':
                out_qp = self.compute_data(fd, 0, **kwargs)

            else:
                raise ValueError('unsupported term mode in %s! (%s)'
                                 % (self.name, term_mode))

            return self.integrate, out_qp, vgv, 1

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, virtual, state, state_p,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = dim * (dim + 1) / 2

        return (n_el, 1, sym, 1), state.dtype

class VolumeTLTerm(HyperElasticTLBase):
    r"""
    Volume term (weak form) in the total Lagrangian formulation.

    :Definition:

    .. math::
         \begin{array}{l}
         \int_{\Omega} q J(\ul{u}) \\
         \mbox{volume mode: vector for } K \from \Ical_h: \int_{T_K}
         J(\ul{u}) \\
         \mbox{rel\_volume mode: vector for } K \from \Ical_h:
         \int_{T_K} J(\ul{u}) / \int_{T_K} 1
         \end{array}

    :Arguments:
        - virtual : :math:`q`
        - state   : :math:`\ul{u}`
    """
    name = 'dw_tl_volume'
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual' : (1, None), 'state' : 'D'}
    family_data_names = ['mtx_f', 'det_f', 'sym_inv_c']

    function = staticmethod(terms.dw_tl_volume)

    def get_fargs(self, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgs, _ = self.get_mapping(virtual)
        vgv, _ = self.get_mapping(state)

        fd = self.get_family_data(state, 'tl_common', self.family_data_names)

        if mode == 'weak':
            if diff_var is None:
                fmode = 0

            else:
                fmode = 1

        elif (mode == 'eval') or (mode == 'el_avg'):
            if term_mode == 'volume':
                fmode = 2

            elif term_mode == 'rel_volume':
                fmode = 3

            else:
                raise ValueError('unsupported term evaluation mode in %s! (%s)'
                                 % (self.name, term_mode))

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

        return fd.mtx_f, fd.sym_inv_c, fd.det_f, vgs, vgv, 0, fmode

    def get_eval_shape(self, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

class DiffusionTLTerm(HyperElasticTLBase):
    r"""
    Diffusion term in the total Lagrangian formulation with
    linearized deformation-dependent permeability
    :math:`\ull{K}(\ul{u}) = J \ull{F}^{-1} \ull{k} f(J) \ull{F}^{-T}`,
    where :math:`\ul{u}` relates to the previous time step :math:`(n-1)`
    and
    :math:`f(J) = \max\left(0, \left(1 + \frac{(J - 1)}{N_f}\right)\right)^2`
    expresses the dependence on volume compression/expansion.

    :Definition:

    .. math::
        \int_{\Omega} \ull{K}(\ul{u}^{(n-1)}) : \pdiff{q}{\ul{X}}
        \pdiff{p}{\ul{X}}

    :Arguments:
        - material_1 : :math:`\ull{k}`
        - material_2 : :math:`N_f`
        - virtual    : :math:`q`
        - state      : :math:`p`
        - parameter  : :math:`\ul{u}^{(n-1)}`
    """
    name = 'dw_tl_diffusion'
    arg_types = ('material_1', 'material_2', 'virtual', 'state', 'parameter')
    arg_shapes = {'material_1' : 'D, D', 'material_2' : '1, 1',
                  'virtual' : (1, 'state'), 'state' : 1, 'parameter' : 'D'}
    family_data_names = ['mtx_f', 'det_f']

    function = staticmethod(terms.dw_tl_diffusion)

    def get_fargs(self, perm, ref_porosity, virtual, state, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vgv, _ = self.get_mapping(parameter)

        fd = self.get_family_data(parameter, 'tl_common',
                                  self.family_data_names)
        grad = self.get(state, 'grad')

        if mode == 'weak':
            if diff_var is None:
                fmode = 0

            else:
                fmode = 1

        elif mode == 'el_avg':
            if term_mode == 'diffusion_velocity':
                fmode = 2

            else:
                raise ValueError('unsupported term evaluation mode in %s! (%s)'
                                 % (self.name, term_mode))

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

        return grad, perm, ref_porosity, fd.mtx_f, fd.det_f, vgv, fmode

    def get_eval_shape(self, perm, ref_porosity, virtual, state, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, dim, 1), state.dtype

class HyperElasticSurfaceTLBase(HyperElasticBase):
    """
    Base class for all hyperelastic surface terms in TL formulation family.
    """
    family_function = staticmethod(terms.dq_tl_finite_strain_surface)
    fd_cache_name = 'tl_surface_common'

    def compute_family_data(self, state):
        ap, sg = self.get_approximation(state)
        sd = ap.surface_data[self.region.name]

        vec = self.get_vector(state)

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        shapes = {
            'mtx_f' : (n_el, n_qp, dim, dim),
            'det_f' : (n_el, n_qp, 1, 1),
            'inv_f' : (n_el, n_qp, dim, dim),
        }
        data = Struct(name='tl_surface_family_data')
        for key, shape in shapes.iteritems():
            setattr(data, key, nm.zeros(shape, dtype=nm.float64))

        self.family_function(data.mtx_f,
                             data.det_f,
                             data.inv_f,
                             vec, sg, sd.fis, ap.econn)
        return data

class SurfaceFluxTLTerm(HyperElasticSurfaceTLBase):
    r"""
    Surface flux term in the total Lagrangian formulation, consistent with
    :class:`DiffusionTLTerm`.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{\nu} \cdot \ull{K}(\ul{u}^{(n-1)}) \pdiff{p}{\ul{X}}

    :Arguments:
        - material_1 : :math:`\ull{k}`
        - material_2 : :math:`N_f`
        - parameter_1 : :math:`p`
        - parameter_2 : :math:`\ul{u}^{(n-1)}`
    """
    name = 'd_tl_surface_flux'
    arg_types = ('material_1', 'material_2', 'parameter_1', 'parameter_2')
    arg_shapes = {'material_1' : 'D, D', 'material_2' : '1, 1',
                  'parameter_1' : 1, 'parameter_2' : 'D'}
    family_data_names = ['det_f', 'inv_f']
    integration = 'surface_extra'

    function = staticmethod(terms.d_tl_surface_flux)

    def get_fargs(self, perm, ref_porosity, pressure, displacement,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, sg = self.get_approximation(displacement)
        fd = self.get_family_data(displacement, 'tl_surface_common',
                                  self.family_data_names)
        grad = self.get(pressure, 'grad')

        fmode = {'eval' : 0, 'el_avg' : 1, 'el' : 0}.get(mode, 0)

        return grad, perm, ref_porosity, fd.inv_f, fd.det_f, sg, fmode

    def get_eval_shape(self, perm, ref_porosity, pressure, displacement,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_en, n_c = self.get_data_shape(displacement)

        return (n_fa, 1, 1, 1), pressure.dtype

class SurfaceTractionTLTerm(HyperElasticSurfaceTLBase):
    r"""
    Surface traction term in the total Lagrangian formulation, expressed
    using :math:`\ul{\nu}`, the outward unit normal vector w.r.t. the
    undeformed surface, :math:`\ull{F}(\ul{u})`, the deformation gradient,
    :math:`J = \det(\ull{F})`, and :math:`\ull{\sigma}` a given traction,
    often equal to a given pressure, i.e.
    :math:`\ull{\sigma} = \pi \ull{I}`.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{\nu} \cdot \ull{F}^{-1} \cdot \ull{\sigma} \cdot
        \ul{v} J

    :Arguments:
        - material : :math:`\ull{\sigma}`
        - virtual  : :math:`\ul{v}`
        - state    : :math:`\ul{u}`
    """
    name = 'dw_tl_surface_traction'
    arg_types = ('opt_material', 'virtual', 'state')
    arg_shapes = [{'opt_material' : 'D, D', 'virtual' : ('D', 'state'),
                   'state' : 'D'},
                  {'opt_material' : None}]
    family_data_names = ['det_f', 'inv_f']
    integration = 'surface_extra'

    function = staticmethod(terms.dw_tl_surface_traction)

    def check_shapes(self, mat, virtual, state):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        if mat is not None:
            assert_(mat.shape == (n_el, n_qp, dim, dim))

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, sg = self.get_approximation(virtual)
        sd = ap.surface_data[self.region.name]
        bf = ap.get_base(sd.bkey, 0, self.integral)

        fd = self.get_family_data(state, 'tl_surface_common',
                                  self.family_data_names)

        if mat is None:
            eye = nm.eye(sg.dim, dtype=nm.float64)
            mat = nm.tile(eye, ((1, sg.n_qp, 1, 1)))

        if diff_var is None:
            fmode = 0

        else:
            fmode = 1

        return mat, fd.det_f, fd.inv_f, bf, sg, sd.fis, fmode

class VolumeSurfaceTLTerm(HyperElasticSurfaceTLBase):
    r"""
    Volume of a :math:`D`-dimensional domain, using a surface integral in the
    total Lagrangian formulation, expressed using :math:`\ul{\nu}`, the outward
    unit normal vector w.r.t. the undeformed surface, :math:`\ull{F}(\ul{u})`,
    the deformation gradient, and :math:`J = \det(\ull{F})`. Uses the
    approximation of :math:`\ul{u}` for the deformed surface coordinates
    :math:`\ul{x}`.

    :Definition:

    .. math::
        1 / D \int_{\Gamma} \ul{\nu} \cdot \ull{F}^{-1} \cdot \ul{x} J

    :Arguments:
        - parameter : :math:`\ul{u}`
    """
    name = 'd_tl_volume_surface'
    arg_types = ('parameter',)
    arg_shapes = {'parameter' : 'D'}
    family_data_names = ['det_f', 'inv_f']
    integration = 'surface_extra'

    function = staticmethod(terms.d_tl_volume_surface)

    def check_shapes(self, parameter):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        assert_(dim == n_c)

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, sg = self.get_approximation(parameter)
        sd = ap.surface_data[self.region.name]
        bf = ap.get_base(sd.bkey, 0, self.integral)

        fd = self.get_family_data(parameter, 'tl_surface_common',
                                  self.family_data_names)

        asc = nm.ascontiguousarray

        coors0 = parameter.field.get_coor()
        coors = asc(coors0 + parameter().reshape(coors0.shape))

        return coors, fd.det_f, fd.inv_f, bf, sg, asc(sd.econn)

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, 1, 1, 1), parameter.dtype
