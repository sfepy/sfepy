import numpy as nm
from sfepy.terms.terms import Term, terms
from sfepy.base.base import Struct
from sfepy import site_config


_msg_missing_data = 'missing family data!'

class HyperElasticFamilyData(Struct):
    """
    Base class for hyperelastic family data.

    The common (family) data are cached in the evaluate cache of state
    variable.
    """
    data_shapes = {
            'mtx_f': ('n_el', 'n_qp', 'dim', 'dim'),
            'det_f': ('n_el', 'n_qp', 1, 1),
            'sym_b': ('n_el', 'n_qp', 'sym', 1),
            'tr_b': ('n_el', 'n_qp', 1, 1),
            'in2_b': ('n_el', 'n_qp', 1, 1),
            'sym_c' : ('n_el', 'n_qp', 'sym', 1),
            'tr_c' : ('n_el', 'n_qp', 1, 1),
            'in2_c' : ('n_el', 'n_qp', 1, 1),
            'sym_inv_c' : ('n_el', 'n_qp', 'sym', 1),
            'green_strain': ('n_el', 'n_qp', 'sym', 1),
            'inv_f': ('n_el', 'n_qp', 'dim', 'dim'),
    }

    def __init__(self, **kwargs):
        Struct.__init__(self, **kwargs)

    def init_data_struct(self, state_shape, name='family_data'):
        from sfepy.mechanics.tensors import dim2sym

        n_el, n_qp, dim, n_en, n_c = state_shape
        sym = dim2sym(dim); sym

        shdict = dict(( (k, v) for k, v in locals().items()\
            if k in ['n_el', 'n_qp', 'dim', 'n_en', 'n_c', 'sym']))

        data = Struct(name=name)
        setattr(data, 'names', self.data_names)
        for key in self.data_names:
            shape = [shdict[sh] if type(sh) is str else sh\
                 for sh in self.data_shapes[key]]
            setattr(data, key, nm.zeros(shape, dtype=nm.float64))

        return data

    def __call__(self, state, region, integral, geometry_type,
                 step=0, derivative=None):
        step_cache = state.evaluate_cache.setdefault(self.cache_name, {})
        cache = step_cache.setdefault(step, {})

        key = (region.name, integral.order, geometry_type)
        data_key = key + (derivative,)

        if data_key in cache:
            data = cache[data_key]

        else:
            vg, _ = state.field.get_mapping(region,
                                            integral, geometry_type[0],
                                            get_saved=True)

            vec = state(step=step, derivative=derivative)
            econn = state.field.get_econn('cell', region)

            st_shape = state.get_data_shape(integral, geometry_type[0],
                                            region.name)
            data = self.init_data_struct(st_shape)

            fargs = tuple([getattr(data, k) for k in self.data_names])
            fargs = fargs + (vec, vg, econn)
            fargs = Term.translate_fargs_mapping(self.family_function,
                                                 list(fargs))

            self.family_function(*fargs)
            cache[data_key] = data

            det_f = data.det_f
            # Minimum over quadrature points.
            jmin = nm.min(det_f, axis=1, keepdims=True)
            jneg = jmin < 0.0
            is_warp = jneg.any()
            if is_warp:
                if site_config.debug_warped_cells():
                    from sfepy.base.base import output
                    from sfepy.discrete.fem import extend_cell_data

                    output('warped (negative volume) elements:')
                    output(nm.unique(nm.nonzero(jneg[:, 0, 0, 0])[0]))

                    mesh = region.domain.mesh
                    jmin = extend_cell_data(jmin,
                                            region.domain, region.name,
                                            val=jmin.max())
                    jneg = extend_cell_data(jneg.astype(nm.float64),
                                            region.domain, region.name,
                                            val=0.0)
                    out = {
                        'jmin' : Struct(name='output_data',
                                        mode='cell', data=jmin),
                        'jneg' : Struct(name='output_data',
                                        mode='cell', data=jneg),
                    }
                    mesh.write('warped_cells.vtk', io='auto', out=out)
                    raise RuntimeError(
                        "inspect 'warped_cells.vtk' for negative volume cells"
                    )

                raise ValueError('warp violation!')

        return data

class HyperElasticBase(Term):
    """
    Base class for all hyperelastic terms in TL/UL formulation.

    `HyperElasticBase.__call__()` computes element contributions given either
    stress (-> residual) or tangent modulus (-> tangent sitffnes matrix),
    i.e. constitutive relation type (CRT) related data. The CRT data are
    computed in subclasses implementing particular CRT (e.g. neo-Hookean
    material), in self.compute_crt_data().

    Modes:

      - 0: total formulation
      - 1: updated formulation

    Notes
    -----
    This is not a proper Term!
    """
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material' : '1, 1', 'virtual' : ('D', 'state'),
                  'state' : 'D'}
    integration = ('cell', 'facet_extra')

    @staticmethod
    def integrate(out, val_qp, vg, fmode):
        if fmode == 2:
            out[:] = val_qp
            status = 0

        else:
            status = vg.integrate(out, val_qp, fmode)

        return status

    @staticmethod
    def function(out, fun, *args):
        return fun(out, *args)

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.stress_cache = None

    def compute_stress(self, mat, family_data, **kwargs):
        out = nm.empty_like(family_data.green_strain)

        get = family_data.get
        fargs = [get(name, msg_if_none=_msg_missing_data)
                 for name in self.family_data_names]

        self.stress_function(out, mat, *fargs, **kwargs)

        return out

    def compute_tan_mod(self, mat, family_data, **kwargs):
        shape = list(family_data.green_strain.shape)
        shape[-1] = shape[-2]
        out = nm.empty(shape, dtype=nm.float64)

        get = family_data.get
        fargs = [get(name, msg_if_none=_msg_missing_data)
                 for name in self.family_data_names]

        self.tan_mod_function(out, mat, *fargs, **kwargs)

        return out

    def get_fargs(self, mat, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg, _ = self.get_mapping(state)

        name = state.name
        fd = self.get_family_data(state, self.region, self.integral,
                                  self.geometry_types[name],
                                  self.arg_steps[name],
                                  self.arg_derivatives[name])

        if mode == 'weak':
            if diff_var is None:
                stress = self.compute_stress(mat, fd, **kwargs)
                self.stress_cache = stress
                tan_mod = nm.array([0], ndmin=4, dtype=nm.float64)

                fmode = 0

            else:
                stress = self.stress_cache
                if stress is None:
                    stress = self.compute_stress(mat, fd, **kwargs)

                tan_mod = self.compute_tan_mod(mat, fd, **kwargs)
                fmode = 1

            return (self.weak_function,
                    stress, tan_mod, fd.mtx_f, fd.det_f, vg, fmode,
                    self.hyperelastic_mode)

        elif mode in ('el_avg', 'qp'):
            if term_mode == 'strain':
                out_qp = fd.green_strain

            elif term_mode == 'stress':
                out_qp = self.compute_stress(mat, fd, **kwargs)

            else:
                raise ValueError('unsupported term mode in %s! (%s)'
                                 % (self.name, term_mode))

            fmode = {'el_avg' : 1, 'qp' : 2}[mode]

            return self.integrate, out_qp, vg, fmode

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.mechanics.tensors import dim2sym

        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = dim2sym(dim)

        if mode != 'qp':
            n_qp = 1

        return (n_el, n_qp, sym, 1), state.dtype

class DeformationGradientTerm(Term):
    r"""
    Deformation gradient :math:`\ull{F}` in quadrature points for
    `term_mode='def_grad'` (default) or the jacobian :math:`J` if
    `term_mode='jacobian'`.

    Supports 'eval', 'el_avg' and 'qp' evaluation modes.

    :Definition:

    .. math::
        \ull{F} = \pdiff{\ul{x}}{\ul{X}}|_{qp}
        = \ull{I} + \pdiff{\ul{u}}{\ul{X}}|_{qp} \;, \\
        \ul{x} = \ul{X} + \ul{u} \;, J = \det{(\ull{F})}

    :Arguments:
        - parameter : :math:`\ul{u}`
    """
    name = 'ev_def_grad'
    arg_types = ('parameter',)
    arg_shapes = {'parameter' : 'D'}
    integration = ('cell', 'facet_extra')

    @staticmethod
    def function(out, vec, vg, econn, term_mode, fmode):
        d = 1 if term_mode == 'jacobian' else vg.dim
        out_qp = nm.empty((out.shape[0], vg.n_qp, d, d), dtype=out.dtype)

        mode = 1 if term_mode == 'jacobian' else 0
        terms.dq_def_grad(out_qp, vec, vg.cmap, econn, mode)

        if fmode == 2:
            out[:] = out_qp
            status = 0

        else:
            status = vg.integrate(out, out_qp, fmode)

        return status

    def get_fargs(self, parameter,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        get_saved = False if term_mode == 'act_conf' else True
        vg, _ = self.get_mapping(parameter, get_saved)

        vec = self.get_vector(parameter)

        fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)
        econn = parameter.field.get_econn(
            self.get_dof_conn_type(self.arg_names[0]), self.region
        )
        return vec, vg, econn, term_mode, fmode

    def get_eval_shape(self, parameter,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(parameter)

        if mode != 'qp':
            n_qp = 1

        if term_mode == 'jacobian':
            return (n_el, n_qp, 1, 1), parameter.dtype

        else: # 'def_grad'
            return (n_el, n_qp, dim, dim), parameter.dtype
