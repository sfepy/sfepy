import numpy as nm

from sfepy.base.base import assert_
from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms

class ZeroTerm(Term):
    r"""
    A do-nothing term useful for introducing additional variables into the
    equations.

    :Definition:

    .. math::
        0

    :Arguments:
        - virtual  : :math:`q` or :math:`\ul{v}`
        - state : :math:`p` or :math:`\ul{u}`
    """
    name = 'dw_zero'
    arg_types = ('virtual', 'state')
    arg_shapes = {'virtual' : ('N', None), 'state' : 'N'}

    @staticmethod
    def function(out):
        out.fill(0.0)

        return 0

    def get_fargs(self, vvar, svar,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        return ()


class IntegrateTerm(Term):
    r"""
    Evaluate (weighted) variable in a region.

    Depending on evaluation mode, integrate a variable over a region
    ('eval'), average it in elements ('el_avg') or interpolate it into
    quadrature points ('qp'). For a surface region and vector variables,
    setting `term_mode` to `'flux'` leads to computing corresponding fluxes
    for the three modes instead.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\cal{D}} y \mbox{ , } \int_{\cal{D}} \ul{y}
        \mbox{ , } \int_\Gamma \ul{y} \cdot \ul{n}\\
        \int_{\cal{D}} c y \mbox{ , } \int_{\cal{D}} c \ul{y}
        \mbox{ , } \int_\Gamma c \ul{y} \cdot \ul{n} \mbox{ flux }

    :Arguments:
        - material : :math:`c` (optional)
        - parameter : :math:`y` or :math:`\ul{y}`
    """
    name = 'ev_integrate'
    arg_types = ('opt_material', 'parameter')
    arg_shapes = [{'opt_material' : '1, 1', 'parameter' : 'N'},
                  {'opt_material' : None}]
    integration = ('cell', 'facet')

    @staticmethod
    def function(out, val_qp, vg, fmode):
        if fmode == 2:
            out[:] = val_qp
            status = 0

        elif fmode == 5:
            normal = vg.normal
            out[:] = dot_sequences(val_qp, normal)
            status = 0

        else:
            status = vg.integrate(out, val_qp, fmode)

        return status

    def get_fargs(self, material, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        val_qp = self.get(parameter, 'val')
        if material is not None:
            val_qp *= material

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)
        if term_mode == 'flux':
            n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(parameter)
            if n_c == dim:
                fmode += 3

        return val_qp, vg, fmode

    def get_eval_shape(self, material, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        if term_mode == 'flux':
            n_c = 1

        return (n_el, n_qp, n_c, 1), parameter.dtype


class IntegrateOperatorTerm(Term):
    r"""
    Integral of a test function weighted by a scalar function
    :math:`c`.

    :Definition:

    .. math::
        \int_{\cal{D}} q \mbox{ or } \int_{\cal{D}} c q

    :Arguments:
        - material : :math:`c` (optional)
        - virtual  : :math:`q`
    """
    name = 'dw_integrate'
    arg_types = ('opt_material', 'virtual')
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, None)},
                  {'opt_material' : None}]
    integration = ('cell', 'facet')

    @staticmethod
    def function(out, material, bf, geo):
        bf_t = nm.tile(bf.transpose((0, 1, 3, 2)), (out.shape[0], 1, 1, 1))
        bf_t = nm.ascontiguousarray(bf_t)
        if material is not None:
            status = geo.integrate(out, material * bf_t)
        else:
            status = geo.integrate(out, bf_t)
        return status

    def get_fargs(self, material, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        assert_(virtual.n_components == 1)
        geo, _ = self.get_mapping(virtual)

        return material, geo.bf, geo


class VolumeTerm(Term):
    r"""
    Volume or surface of a domain. Uses approximation of the parameter variable.

    :Definition:

    .. math::
        \int_{\cal{D}} 1

    :Arguments:
        - parameter : any variable
    """
    name = 'ev_volume'
    arg_types = ('parameter',)
    arg_shapes = [{'parameter' : 'N'}]
    integration = ('cell', 'facet')

    @staticmethod
    def function(out, geo):
        out[:] = geo.volume

        return 0

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(parameter)

        return geo,

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_cell, n_qp, dim, n_n, n_c = self.get_data_shape(parameter)

        return (n_cell, 1, 1, 1), parameter.dtype


class VolumeSurfaceTerm(Term):
    r"""
    Volume of a :math:`D`-dimensional domain, using a surface integral. Uses
    approximation of the parameter variable.

    :Definition:

    .. math::
        1 / D \int_\Gamma \ul{x} \cdot \ul{n}

    :Arguments:
        - parameter : any variable
    """
    name = 'ev_volume_surface'
    arg_types = ('parameter',)
    arg_shapes = {'parameter' : 'N'}
    integration = 'facet'

    function = staticmethod(terms.d_volume_surface)

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        sd = parameter.field.extra_data[f'sd_{self.region.name}']
        coor = parameter.field.get_coor()

        return coor, sg, sd.econn.copy()

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(parameter)

        return (n_fa, 1, 1, 1), parameter.dtype

class SurfaceMomentTerm(Term):
    r"""
    Surface integral of the outer product of the unit outward normal
    :math:`\ul{n}` and the coordinate :math:`\ul{x}` shifted by :math:`\ul{x}_0`

    :Definition:

    .. math::
        \int_{\Gamma} \ul{n} (\ul{x} - \ul{x}_0)

    :Arguments:
        - material  : :math:`\ul{x}_0` (special)
        - parameter : any variable
    """
    name = 'ev_surface_moment'
    arg_types = ('material', 'parameter')
    arg_shapes = {'material' : '.: D', 'parameter' : 'N'}
    integration = 'facet'

    function = staticmethod(terms.di_surface_moment)

    def get_fargs(self, material, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        sd = parameter.field.extra_data[f'sd_{self.region.name}']
        coor = parameter.field.get_coor() \
               - nm.asarray(material, dtype=nm.float64)[None,:]

        return coor, sg, sd.econn.copy()

    def get_eval_shape(self, material, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(parameter)

        return (n_fa, 1, dim, dim), parameter.dtype

class IntegrateMatTerm(Term):
    r"""
    Evaluate material parameter :math:`m` in a volume region.

    Depending on evaluation mode, integrate a material parameter over a
    volume region ('eval'), average it in elements ('el_avg') or
    interpolate it into volume quadrature points ('qp').

    Uses reference mapping of :math:`y` variable.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \int_{\cal{D}} c

    :Arguments:
        - material  : :math:`c` (can have up to two dimensions)
        - parameter : :math:`y`
    """
    name = 'ev_integrate_mat'
    arg_types = ('material', 'parameter')
    arg_shapes = [{'material' : 'N, N', 'parameter' : 'N'}]
    integration = ('cell', 'facet')

    @staticmethod
    def function(out, mat, geo, fmode):
        if fmode == 2:
            out[:] = mat
            status = 0

        else:
            status = geo.integrate(out, mat, fmode)

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(parameter)

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return mat, geo, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)
        n_row, n_col = mat.shape[-2:]

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, n_row, n_col), mat.dtype


class SumNodalValuesTerm(Term):
    r"""
    Sum nodal values.

    :Arguments:
        - parameter : :math:`p` or :math:`\ul{u}`
    """
    name = 'ev_sum_vals'
    arg_types = ('parameter',)
    arg_shapes = {'parameter' : 'N'}

    @staticmethod
    def function(out, vec):
        out[:] = nm.sum(vec, 0)

        return 0

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vec = parameter.get_state_in_region(self.region)

        return vec,

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, n_c), parameter.dtype
