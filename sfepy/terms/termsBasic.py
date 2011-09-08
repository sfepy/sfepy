import numpy as nm

from sfepy.base.base import assert_
from sfepy.linalg import dot_sequences
from sfepy.terms.terms import Term, terms
from sfepy.terms.terms_base import ScalarScalar, ScalarScalarTH

class IntegrateVolumeTerm(Term):
    r"""
    :Description:
    Depending on evaluation mode, integrate a variable over a volume
    region ('eval'), average it in elements or interpolate it into
    volume quadrature points ('qp').

    :Definition:
    .. math::
        \int_\Omega y \mbox{ , } \int_\Omega \ul{y}

    .. math::
        \mbox{vector for } K \from \Ical_h:
        \int_{T_K} y / \int_{T_K} 1 \mbox{ , }
        \int_{T_K} \ul{y} / \int_{T_K} 1

    .. math::
         y|_{qp} \mbox{ , } \ul{y}|_{qp}

    :Arguments:
        parameter : :math:`y` or :math:`\ul{y}`
    """
    name = 'di_volume_integrate'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, val_qp, vg, fmode):
        if fmode == 2:
            out[:] = val_qp
            status = 0

        else:
            status = vg.integrate(out, val_qp, fmode)

        return status

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(parameter)

        val = self.get(parameter, 'val')

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

        return val, vg, fmode

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, n_c, 1), parameter.dtype

class IntegrateSurfaceTerm(Term):
    r"""
    :Description:
    Depending on evaluation mode, integrate a variable over a surface
    region ('eval'), average it in element faces or interpolate it into
    surface quadrature points ('qp'). For vector variables, setting
    `term_mode` to `'flux'` leads to computing corresponding fluxes for
    the three modes instead.

    :Definition:
    .. math::
        \int_\Gamma y \mbox{ , } \int_\Gamma \ul{y}
        \mbox{ , } \int_\Gamma \ul{y} \cdot \ul{n} \mbox{ or }
        \int_\Gamma c y \mbox{ , } \int_\Gamma c \ul{y}
        \mbox{ , } \int_\Gamma c \ul{y} \cdot \ul{n} \mbox{ flux }

    .. math::
        \mbox{vector for } K \from \Ical_h:
        \int_{T_K} y / \int_{T_K} 1 \mbox{ , }
        \int_{T_K} \ul{y} / \int_{T_K} 1 \mbox{ , }
        \int_{T_K} (\ul{y} \cdot \ul{n}) / \int_{T_K} 1

    .. math::
         y|_{qp} \mbox{ , } \ul{y}|_{qp}
        \mbox{ , } (\ul{y} \cdot \ul{n})|_{qp} \mbox{ flux }

    :Arguments:
        material : :math:`c` (optional),
        parameter : :math:`y` or :math:`\ul{y}`,
    """
    name = 'di_surface_integrate'
    arg_types = ('opt_material', 'parameter')
    integration = 'surface'

    @staticmethod
    def function(out, val_qp, sg, fmode):
        if fmode == 2:
            out[:] = val_qp
            status = 0

        elif fmode == 5:
            normal = sg.variable(0)
            out[:] = dot_sequences(val_qp, normal)
            status = 0

        else:
            status = sg.integrate(out, val_qp, fmode)

        return status

    def get_fargs(self, material, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(parameter)

        val_qp = self.get(parameter, 'val')
        if material is not None:
            val_qp *= material

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)
        if term_mode == 'flux':
            n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(parameter)
            if n_c == dim:
                fmode += 3

        return val_qp, sg, fmode

    def get_eval_shape(self, material, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        if term_mode == 'flux':
            n_c = 1

        return (n_fa, n_qp, n_c, 1), parameter.dtype

class IntegrateVolumeOperatorTerm(Term):
    r"""
    :Description:
    Volume integral of a test function weighted by a scalar function
    :math:`c`.

    :Definition:
    .. math::
        \int_\Omega q \mbox{ or } \int_\Omega c q

    :Arguments:
        material : :math:`c` (optional),
        virtual  : :math:`q`
    """
    name = 'dw_volume_integrate'
    arg_types = ('opt_material', 'virtual')

    @staticmethod
    def function(out, material, bf, geo):
        bf_t = nm.tile(bf.transpose((0, 2, 1)), (out.shape[0], 1, 1, 1))
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

class IntegrateSurfaceOperatorTerm(IntegrateVolumeOperatorTerm):
    r"""
    :Description:
    Surface integral of a test function weighted by a scalar function
    :math:`c`.

    :Definition:
    .. math::
        \int_{\Gamma} q \mbox{ or } \int_\Gamma c q

    :Arguments:
        material : :math:`c` (optional),
        virtual  : :math:`q`
    """
    name = 'dw_surface_integrate'
    arg_types = ('opt_material', 'virtual')
    integration = 'surface'

class DotProductVolumeTerm(Term):
    r"""
    :Description:
    Volume :math:`L^2(\Omega)` dot product for both scalar and vector
    fields.

    :Definition:
    .. math::
        \int_\Omega p r \mbox{ , } \int_\Omega \ul{u} \cdot \ul{w}
        \mbox{ or }\int_\Omega c p r \mbox{ , } \int_\Omega c \ul{u} \cdot \ul{w}

    :Arguments:
        material    : :math:`c` (optional),
        parameter_1 : :math:`p` or :math:`\ul{u}`,
        parameter_2 : :math:`r` or :math:`\ul{w}`
    """
    name = 'd_volume_dot'
    arg_types = ('opt_material', 'parameter_1', 'parameter_2')

    @staticmethod
    def function(out, mat, val1, val2, geo):
        if val1.shape[-1] > 1:
            out_qp = nm.sum(val1 * val2, axis=-1)
        else:
            out_qp = val1 * val2

        if mat is not None:
            status = geo.integrate(out, mat * out_qp)
        else:
            status = geo.integrate(out, out_qp)

        return status

    def get_fargs(self, mat, par1, par2,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(par1)

        val1 = self.get(par1, 'val')
        val2 = self.get(par2, 'val')

        return mat, val1, val2, geo

    def get_eval_shape(self, mat, par1, par2,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_cell, n_qp, dim, n_n, n_c = self.get_data_shape(par1)

        return (n_cell, 1, 1, 1), par1.dtype

class DotProductSurfaceTerm(DotProductVolumeTerm):
    r"""
    :Description:
    Surface :math:`L^2(\Gamma)` dot product for both scalar and vector
    fields.

    :Definition:
    .. math::
        \int_\Gamma p r \mbox{ , } \int_\Gamma \ul{u} \cdot \ul{w}
        \mbox{ or } \int_\Gamma c p r \mbox{ , } \int_\Gamma c \ul{u} \cdot \ul{w}

    :Arguments:
        material    : :math:`c` (optional),
        parameter_1 : :math:`p` or :math:`\ul{u}`,
        parameter_2 : :math:`r` or :math:`\ul{w}`
    """
    name = 'd_surface_dot'
    arg_types = ('opt_material', 'parameter_1', 'parameter_2')
    integration = 'surface'

class VolumeTerm(Term):
    r"""
    :Description:
    Volume of a domain. Uses approximation of the parameter variable.

    :Definition:
    .. math::
        \int_\Omega 1

    :Arguments:
        parameter : any variable
    """
    name = 'd_volume'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, geo):
        out[:] = geo.variable(2)

        return 0

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(parameter)

        return geo,

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_cell, n_qp, dim, n_n, n_c = self.get_data_shape(parameter)

        return (n_cell, 1, 1, 1), parameter.dtype

class SurfaceTerm(VolumeTerm):
    r"""
    :Description:
    Surface of a domain. Uses approximation of the parameter variable.

    :Definition:
    .. math::
        \int_\Gamma 1

    :Arguments:
        parameter : any variable
    """
    name = 'd_surface'
    arg_types = ('parameter',)
    integration = 'surface'

class VolumeSurfaceTerm(Term):
    r"""
    :Description:
    Volume of a domain, using a surface integral. Uses approximation of the
    parameter variable.

    :Definition:
    .. math::
        \int_\Gamma \ul{x} \cdot \ul{n}

    :Arguments:
        parameter : any variable
    """
    name = 'd_volume_surface'
    arg_types = ('parameter',)
    integration = 'surface'

    function = staticmethod(terms.d_volume_surface)

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, sg = self.get_approximation(parameter)

        sd = ap.surface_data[self.region.name]
        coor = parameter.field.get_coor()

        return coor, sg.bf, sg, sd.econn.copy()

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(parameter)

        return (n_fa, 1, 1, 1), parameter.dtype

class SurfaceMomentTerm(Term):
    r"""
    :Description:
    Surface integral of the outer product of the unit outward normal
    :math:`\ul{n}` and the coordinate :math:`\ul{x}` shifted by :math:`\ul{x}_0`

    :Definition:
    .. math::
        \int_{\Gamma} \ul{n} (\ul{x} - \ul{x}_0)

    :Arguments:
        parameter : any variable,
        shift     : :math:`\ul{x}_0`
    """
    name = 'di_surface_moment'
    arg_types = ('parameter', 'shift')
    integration = 'surface'

    function = staticmethod(terms.di_surface_moment)

    def get_fargs(self, parameter, shift,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        ap, sg = self.get_approximation(parameter)

        sd = ap.surface_data[self.region.name]
        coor = parameter.field.get_coor() \
               - nm.asarray(shift, dtype=nm.float64)[None,:]

        return coor, sg.bf, sg, sd.econn.copy()

    def get_eval_shape(self, parameter, shift,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_fa, n_qp, dim, n_fn, n_c = self.get_data_shape(parameter)

        return (n_fa, 1, dim, dim), parameter.dtype

class IntegrateMatTerm(Term):
    r"""
    :Description:
    Material parameter :math:`m` integrated over a volume/surface region
    or averaged in its elements/faces. Uses approximation of :math:`y`
    variable.

    :Definition:
    .. math::
        \int_\Omega m

    .. math::
        \mbox{vector for } K \from \Ical_h: \int_{T_K} m / \int_{T_K} 1

    :Arguments:
        material  : :math:`m` (can have up to two dimensions),
        parameter : :math:`y`
    """
    name = 'di_integrate_mat'
    arg_types = ('material', 'parameter')

    @staticmethod
    def function(out, mat, geo, fmode):
        status = geo.integrate(out, mat, fmode)

        return status

    def get_fargs(self, mat, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        geo, _ = self.get_mapping(parameter)

        fmode = {'eval' : 0, 'el_avg' : 1}.get(mode, 1)

        return mat, geo, fmode

    def get_eval_shape(self, mat, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)
        n_row, n_col = mat.shape[-2:]

        return (n_el, 1, n_row, n_col), mat.dtype

class DotProductVolumeTerm(Term):
    r"""
    :Description:
    Volume :math:`L^2(\Omega)` weighted dot product for both scalar
    and vector (not implemented in weak form!) fields. Can be evaluated. Can
    use derivatives.

    :Definition:
    .. math::
        \int_\Omega q p \mbox{ , } \int_\Omega \ul{v} \cdot \ul{u} \mbox{ , }
        \int_\Omega p r \mbox{ , } \int_\Omega \ul{u} \cdot \ul{w}
        \mbox { or }
        \int_\Omega c q p \mbox{ , } \int_\Omega c \ul{v} \cdot \ul{u} \mbox{ , }
        \int_\Omega c p r \mbox{ , } \int_\Omega c \ul{u} \cdot \ul{w}

    :Arguments 1:
        material : optional weight function :math:`c`,
        virtual  : :math:`q` or :math:`\ul{v}`,
        state    : :math:`p` or :math:`\ul{u}`

    :Arguments 2:
        material    : optional weight function :math:`c`,
        parameter_1 : :math:`p` or :math:`\ul{u}`,
        parameter_2 : :math:`r` or :math:`\ul{w}`
    """
    name = 'dw_volume_dot'
    arg_types = (('opt_material', 'virtual', 'state'),
                 ('opt_material', 'parameter_1', 'parameter_2'))
    modes = ('weak', 'eval')

    @staticmethod
    def dw_volume_dot(out, mat, val_qp, vvg, svg, fmode):
        bf_t = vvg.bf.transpose((0, 2, 1))
        if mat is not None:
            if fmode == 0:
                vec = bf_t * mat * val_qp

            else:
                vec = bf_t * mat * svg.bf

        else:
            if fmode == 0:
                vec = bf_t * val_qp

            else:
                vec = bf_t * svg.bf

        status = vvg.integrate(out, vec)

        return status

    @staticmethod
    def d_volume_dot(out, mat, val1_qp, val2_qp, vg):
        if val1_qp.shape[2] > 1:
            vec = nm.sum(val1_qp * val2_qp, axis=-1)

        else:
            vec = val1_qp * val2_qp

        if mat is not None:
            status = vg.integrate(out, mat * vec)
        else:
            status = vg.integrate(out, vec)

        return status

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vvg, _ = self.get_mapping(virtual)
        svg, _ = self.get_mapping(state)

        if mode == 'weak':
            if state.n_components > 1:
                raise NotImplementedError

            if diff_var is None:
                val_qp = self.get(state, 'val')
                fmode = 0

            else:
                val_qp = nm.array([0], ndmin=4, dtype=nm.float64)
                fmode = 1

            return mat, val_qp, vvg, svg, fmode

        elif mode == 'eval':
            val1_qp = self.get(virtual, 'val')
            val2_qp = self.get(state, 'val')

            return mat, val1_qp, val2_qp, vvg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)

        return (n_el, 1, 1, 1), state.dtype

    def set_arg_types(self):
        if self.mode == 'weak':
            self.function = self.dw_volume_dot

        else:
            self.function = self.d_volume_dot

##
# c: 03.04.2008
class DotSProductVolumeOperatorWTHTerm( ScalarScalarTH, Term ):
    r"""
    :Description:
    Fading memory volume :math:`L^2(\Omega)` weighted dot product for
    scalar fields. Can use derivatives.

    :Definition:
    .. math::
        \int_\Omega \left [\int_0^t \Gcal(t-\tau) p(\tau) \difd{\tau} \right] q

    :Arguments:
        ts       : :class:`TimeStepper` instance,
        material : :math:`\Gcal(\tau)`,
        virtual  : :math:`q`,
        state    : :math:`p`
    """
    name = 'dw_volume_dot_w_scalar_th'
    arg_types = ('ts', 'material', 'virtual', 'state')
    use_caches = {'state_in_volume_qp' : [['state', {'state' : (-1,-1)}]]}

    function = staticmethod(terms.dw_volume_wdot_scalar)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mats, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        if (ts.step == 0) and (mode == 0):
            raise StopIteration

        n_el, n_qp = self.data_shape[:2]
        bf = ap.get_base('v', 0, self.integral)

        if mode == 1:
            aux = nm.array( [0], ndmin = 4, dtype = nm.float64 )
            mat = mats[0]
            mat = nm.tile(mat, (n_el, n_qp, 1, 1))
            return (ts.dt, aux, bf, mat, vg), shape, mode

        else:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            def iter_kernel():
                for ii, mat in enumerate( mats ):
                    vec_qp = cache('state', self, ii,
                                   state=state, get_vector=self.get_vector)
                    mat = nm.tile(mat, (n_el, n_qp, 1, 1))
                    yield ii, (ts.dt, vec_qp, bf, mat, vg)
            return iter_kernel, shape, mode

class DotSProductVolumeOperatorWETHTerm( ScalarScalar, Term ):
    r"""
    :Description:
    Fading memory volume :math:`L^2(\Omega)` weighted dot product for
    scalar fields. This term has the same definition as
    dw_volume_dot_w_scalar_th, but assumes an exponential approximation of
    the convolution kernel resulting in much higher efficiency. Can use
    derivatives.

    :Definition:
    .. math::
        \int_\Omega \left [\int_0^t \Gcal(t-\tau) p(\tau) \difd{\tau} \right] q

    :Arguments:
        ts         : :class:`TimeStepper` instance,
        material_0 : :math:`\Gcal(0)`,
        material_1 : :math:`\exp(-\lambda \Delta t)` (decay at :math:`t_1`),
        virtual    : :math:`q`,
        state      : :math:`p`
    """
    name = 'dw_volume_dot_w_scalar_eth'
    arg_types = ('ts', 'material_0', 'material_1', 'virtual', 'state')
    use_caches = {'state_in_volume_qp' : [['state']],
                  'exp_history' : [['material_0', 'material_1', 'state']]}

    function = staticmethod(terms.dw_volume_wdot_scalar)

    def get_fargs( self, diff_var = None, chunk_size = None, **kwargs ):
        ts, mat0, mat1, virtual, state = self.get_args( **kwargs )
        ap, vg = self.get_approximation(virtual)

        self.set_data_shape( ap )
        shape, mode = self.get_shape( diff_var, chunk_size )

        bf = ap.get_base('v', 0, self.integral)
        if diff_var is None:
            cache = self.get_cache( 'state_in_volume_qp', 0 )
            vec_qp = cache('state', self, 0,
                           state=state, get_vector=self.get_vector)

            cache = self.get_cache('exp_history', 0)
            increment = cache('increment', self, 0,
                              decay=mat1, values=vec_qp)
            history = cache('history', self, 0)

            fargs = (ts.dt, history + increment, bf, mat0, vg)
            if ts.step == 0: # Just init the history in step 0.
                raise StopIteration

        else:
            aux = nm.array([0], ndmin=4, dtype=nm.float64)
            fargs = (ts.dt, aux, bf, mat0, vg)

        return fargs, shape, mode

class SumNodalValuesTerm(Term):
    r"""
    :Description:
    Sum nodal values.

    :Arguments:
        parameter : :math:`p` or :math:`\ul{u}`,
    """
    name = 'd_sum_vals'
    arg_types = ('parameter',)

    @staticmethod
    def function(out, vec):
        out[:] = vec

        return 0

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vec = parameter.get_state_in_region(self.region)

        return vec

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        return (n_el, n_c), parameter.dtype
