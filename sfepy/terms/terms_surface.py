import numpy as nm

from sfepy.base.base import assert_
from sfepy.terms.terms import Term, terms
from sfepy.linalg import dot_sequences
from sfepy.mechanics.contact_bodies import ContactPlane, ContactSphere
from sfepy.mechanics.tensors import get_full_indices
from sfepy.discrete.common.extmods._geommech import geme_mulAVSB3py

##
# 22.08.2006, c
class LinearTractionTerm(Term):
    r"""
    Linear traction forces, where, depending on dimension of
    'material' argument, :math:`\ull{\sigma} \cdot \ul{n}` is
    :math:`\bar{p} \ull{I} \cdot \ul{n}` for a given scalar pressure,
    :math:`\ul{f}` for a traction vector, and itself for a stress tensor.

    The material parameter can have one of the following shapes: 1 or (1, 1),
    (D, 1), (S, 1) in all modes, or (D, D) in the `eval` mode only.
    The symmetric tensor storage (S, 1) is as follows: in 3D S = 6
    and the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`,
    in 2D S = 3 and the indices ordered as :math:`[11, 22, 12]`.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot \ull{\sigma} \cdot \ul{n},
        \int_{\Gamma} \ul{v} \cdot \ul{n},

    :Arguments:
        - material : :math:`\ull{\sigma}`
        - virtual  : :math:`\ul{v}`
    """
    name = 'dw_surface_ltr'
    arg_types = (('opt_material', 'virtual'),
                 ('opt_material', 'parameter'))
    arg_shapes = [{'opt_material' : 'S, 1', 'virtual' : ('D', None),
                   'parameter' : 'D'},
                  {'opt_material' : 'D, 1'}, {'opt_material' : '1, 1'},
                  {'opt_material' : 'D, D'}, {'opt_material' : None}]
    modes = ('weak', 'eval')
    integration = 'facet'

    @staticmethod
    def d_fun(out, traction, val, sg):
        tdim, tdim2 = traction.shape[2:]
        dim = val.shape[2]
        sym = (dim + 1) * dim // 2

        if tdim == 0:
            aux = dot_sequences(val, sg.normal, 'ATB')

        elif tdim == 1: # Pressure
            aux = dot_sequences(val, traction * sg.normal, 'ATB')

        elif tdim == dim and tdim2 == 1: # Traction vector
            aux = dot_sequences(val, traction, 'ATB')

        elif tdim == dim and tdim2 == dim: # Traction tensor - nonsymmetric
            trn = dot_sequences(traction, sg.normal, 'ATB')
            aux = dot_sequences(val, trn, 'ATB')

        elif tdim == sym: # Traction tensor
            trn, ret = geme_mulAVSB3py(traction, sg.normal)
            aux = dot_sequences(val, trn, 'ATB')

        status = sg.integrate(out, aux)
        return status

    def get_fargs(self, traction, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        if traction is not None:
            n_el, _, _, _, _ = self.get_data_shape(virtual)
            traction = Term.tile_mat(traction, n_el)

        if traction is None:
            traction = nm.zeros((0,0,0,0), dtype=nm.float64)

        if mode == 'weak':
            return traction, sg

        elif mode == 'eval':
            val = self.get(virtual, 'val')
            return traction, val, sg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, traction, virtual,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(virtual)

        return (n_el, 1, 1, 1), virtual.dtype

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = terms.dw_surface_ltr

        else:
            self.function = self.d_fun


class SDLinearTractionTerm(Term):
    r"""
    Sensitivity of the linear traction term.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot (\ull{\sigma}\, \ul{n}),
        \int_{\Gamma} \ul{v} \cdot \ul{n},

    :Arguments:
        - material  : :math:`\ull{\sigma}`
        - parameter : :math:`\ul{v}`
    """

    name = 'ev_sd_surface_ltr'
    arg_types = ('opt_material', 'parameter', 'parameter_mv')
    arg_shapes = [{'opt_material': 'S, 1', 'parameter': 'D',
                   'parameter_mv': 'D'}, {'opt_material': '1, 1'},
                  {'opt_material': 'D, 1'}, {'opt_material': 'D, D'},
                  {'opt_material': None}]
    integration = 'facet'

    @staticmethod
    def d_fun(out, traction, val, grad_mv, div_mv, sg):

        tdim, tdim2 = (None, None) if traction is None else traction.shape[2:]
        dim = val.shape[2]
        sym = (dim + 1) * dim // 2
        n_el, n_qp = div_mv.shape[:2]

        val2 = sg.normal

        if tdim is None:
            trac = nm.tile(nm.eye(dim), (n_el, n_qp, 1, 1))

        elif tdim == 1:
            trac = nm.tile(nm.eye(dim), (n_el, n_qp, 1, 1)) * traction

        elif tdim == dim and tdim2 == 1:  # Traction vector
            trac = nm.tile(nm.eye(dim), (n_el, n_qp, 1, 1))
            val2 = traction

        elif tdim == dim and tdim2 == dim:  # Traction tensor - nonsymmetric
            trac = traction

        elif tdim == sym:  # Traction tensor
            remap = nm.array(get_full_indices(dim)).flatten()
            trac = traction[..., remap, :].reshape((-1, n_qp, dim, dim))

        sa_trac = trac * div_mv
        sa_trac -= nm.einsum('qpik,qpkj->qpij', trac, grad_mv,
                             optimize='greedy')

        aux = dot_sequences(val, sa_trac, 'ATB')
        aux = dot_sequences(aux, val2, 'AB')
        status = sg.integrate(out, aux)
        return status

    def get_fargs(self, traction, par_u, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(par_u)

        val = self.get(par_u, 'val')
        grad_mv = self.get(par_mv, 'grad', integration='facet_extra')
        div_mv = self.get(par_mv, 'div', integration='facet_extra')

        return traction, val, grad_mv, div_mv, sg

    def get_eval_shape(self, traction, par_u, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_u)

        return (n_el, 1, 1, 1), par_u.dtype

    def set_arg_types(self):
        self.function = self.d_fun


class ContactPlaneTerm(Term):
    r"""
    Small deformation elastic contact plane term with penetration penalty.

    The plane is given by an anchor point :math:`\ul{A}` and a normal
    :math:`\ul{n}`. The contact occurs in points that orthogonally project onto
    the plane into a polygon given by orthogonal projections of boundary points
    :math:`\{\ul{B}_i\}`, :math:`i = 1, \dots, N_B` on the plane. In such
    points, a penetration distance :math:`d(\ul{u}) = (\ul{X} + \ul{u} -
    \ul{A}, \ul{n})` is computed, and a force :math:`f(d(\ul{u})) \ul{n}` is
    applied. The force depends on the non-negative parameters :math:`k`
    (stiffness) and :math:`f_0` (force at zero penetration):

    - If :math:`f_0 = 0`:

      .. math::

         f(d) = 0 \mbox{ for } d \leq 0 \;, \\
         f(d) = k d \mbox{ for } d > 0 \;.

    - If :math:`f_0 > 0`:

      .. math::

         f(d) = 0 \mbox{ for } d \leq -\frac{2 r_0}{k} \;, \\
         f(d) = \frac{k^2}{4 r_0} d^2 + k d + r_0
         \mbox{ for } -\frac{2 r_0}{k} < d \leq 0 \;, \\
         f(d) = k d + f_0 \mbox{ for } d > 0 \;.

      In this case the dependence :math:`f(d)` is smooth, and a (small) force
      is applied even for (small) negative penetrations: :math:`-\frac{2
      r_0}{k} < d \leq 0`.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot f(d(\ul{u})) \ul{n}

    :Arguments:
        - material_f : :math:`[k, f_0]`
        - material_n : :math:`\ul{n}` (special)
        - material_a : :math:`\ul{A}` (special)
        - material_b : :math:`\{\ul{B}_i\}`, :math:`i = 1, \dots, N_B` (special)
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_contact_plane'
    arg_types = ('material_f', 'material_n', 'material_a', 'material_b',
                 'virtual', 'state')
    arg_shapes = {'material_f' : '1, 2', 'material_n' : '.: D',
                  'material_a' : '.: D', 'material_b' : '.: N, D',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['3_4', '3_8']
    integration = 'facet'

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.cp = None

    @staticmethod
    def function(out, force, normal, geo, fmode):
        bf = geo.bf[0]
        nbf = bf * normal[None, :, None]
        nbf.shape = (bf.shape[0], bf.shape[2] * normal.shape[0])

        if fmode == 0:
            out_qp = force * nbf[None, :, :, None]

        else:
            nbf2 = nbf[:, :, None] * nbf[:, None, :]
            out_qp = force * nbf2[None, :, :, :]

        status = geo.integrate(out, nm.ascontiguousarray(out_qp))

        return status

    @staticmethod
    def smooth_f(d, k, f0, a, eps, diff):
        ii = nm.where((d > eps) & (d <= 0.0))[0]

        if not diff:
            out = nm.where(d > 0.0, k * d + f0, 0.0)
            if len(ii):
                di = d[ii]
                out[ii] = a[ii] * di**2 + k[ii] * di + f0[ii]

        else:
            out = nm.where(d > 0.0, k, 0.0)
            if len(ii):
                out[ii] = 2 * a[ii] * d[ii] + k[ii]

        return out

    @staticmethod
    def _get_force_pars(force_pars, shape):
        k = force_pars[..., 0]
        f0 = force_pars[..., 1]
        k.shape = f0.shape = shape

        ir = f0 >= 1e-14
        eps = nm.where(ir, - 2.0 * f0 / k, 0.0)
        a = nm.zeros_like(eps)
        a[ir] = k[ir]**2 / (4.0 * f0[ir])

        return k, f0, eps, a

    def get_fargs(self, force_pars, normal, anchor, bounds, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        assert_((force_pars >= 0.0).all(),
                'force parameters must be non-negative!')

        force_pars = Term.tile_mat(force_pars, sg.n_el)

        if self.cp is None:
            self.cp = ContactPlane(anchor, normal, bounds)

        qps = self.get_physical_qps()
        qp_coors = qps.values
        u_qp = self.get(state, 'val').reshape(qp_coors.shape)

        # Deformed QP coordinates.
        coors = u_qp + qp_coors

        force = nm.zeros(coors.shape[0], dtype=nm.float64)

        k, f0, eps, a = self._get_force_pars(force_pars, force.shape)

        # Active points in contact change with displacements!
        ii = self.cp.mask_points(qp_coors)

        if diff_var is None:
            if ii.any():
                d = self.cp.get_distance(coors[ii])
                # Force in the plane normal direction.
                force[ii] = self.smooth_f(d, k[ii], f0[ii], a[ii], eps[ii], 0)

            fmode = 0

        else:
            if ii.any():
                d = self.cp.get_distance(coors[ii])
                # Force in the plane normal direction derivative.
                force[ii] = self.smooth_f(d, k[ii], f0[ii], a[ii], eps[ii], 1)

            fmode = 1

        force.shape = qps.shape[:2] + (1, 1)

        return force, self.cp.normal, sg, fmode

class ContactSphereTerm(ContactPlaneTerm):
    r"""
    Small deformation elastic contact sphere term with penetration penalty.

    The sphere is given by a centre point :math:`\ul{C}` and a radius
    :math:`R`. The contact occurs in points that are closer to :math:`\ul{C}`
    than :math:`R`. In such points, a penetration distance :math:`d(\ul{u}) =
    R - ||\ul{X} + \ul{u} - \ul{C}||` is computed, and a force
    :math:`f(d(\ul{u})) \ul{n}(\ul{u})` is applied, where :math:`\ul{n}(\ul{u})
    = (\ul{X} + \ul{u} - \ul{C}) / ||\ul{X} + \ul{u} - \ul{C}||`. The force
    depends on the non-negative parameters :math:`k` (stiffness) and
    :math:`f_0` (force at zero penetration):

    - If :math:`f_0 = 0`:

      .. math::

         f(d) = 0 \mbox{ for } d \leq 0 \;, \\
         f(d) = k d \mbox{ for } d > 0 \;.

    - If :math:`f_0 > 0`:

      .. math::

         f(d) = 0 \mbox{ for } d \leq -\frac{2 r_0}{k} \;, \\
         f(d) = \frac{k^2}{4 r_0} d^2 + k d + r_0
         \mbox{ for } -\frac{2 r_0}{k} < d \leq 0 \;, \\
         f(d) = k d + f_0 \mbox{ for } d > 0 \;.

      In this case the dependence :math:`f(d)` is smooth, and a (small) force
      is applied even for (small) negative penetrations: :math:`-\frac{2
      r_0}{k} < d \leq 0`.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot f(d(\ul{u})) \ul{n}(\ul{u})

    :Arguments:
        - material_f : :math:`[k, f_0]`
        - material_c : :math:`\ul{C}` (special)
        - material_r : :math:`R` (special)
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_contact_sphere'
    arg_types = ('material_f', 'material_c', 'material_r', 'virtual', 'state')
    arg_shapes = {'material_f' : '1, 2', 'material_c' : '.: D',
                  'material_r' : '.: 1',
                  'virtual' : ('D', 'state'), 'state' : 'D'}
    geometries = ['3_4', '3_8']
    integration = 'facet'

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        self.cs = None

    @staticmethod
    def function(out, force, normals, fd, geo, fmode):
        bf = geo.bf[0]
        nbf = bf * normals
        nbf.shape = (normals.shape[0], normals.shape[1],
                     bf.shape[2] * normals.shape[2])

        if fmode == 0:
            out_qp = force * nbf[..., None]

        else:
            nbf2 = nbf[..., None] * nbf[..., None, :]

            dim = normals.shape[2]
            n_ep = bf.shape[2]
            bb = bf[:, 0]
            vb = nm.zeros((bf.shape[0], dim, dim * n_ep))
            for ii in range(dim):
                vb[:, ii, ii*n_ep:(ii+1)*n_ep] = bb
            ee = nm.eye(3)[None, ...]
            eebf2 = dot_sequences(vb, dot_sequences(ee, vb), 'ATB')

            out_qp = force * nbf2
            if fd is not None:
                out_qp -= fd * (eebf2[None, :, :, :] - nbf2)

        status = geo.integrate(out, nm.ascontiguousarray(out_qp))

        return status

    def get_fargs(self, force_pars, centre, radius, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        assert_((force_pars >= 0.0).all(),
                'force parameters must be non-negative!')

        force_pars = Term.tile_mat(force_pars, sg.n_el)

        if self.cs is None:
            self.cs = ContactSphere(centre, radius)

        qps = self.get_physical_qps()
        qp_coors = qps.values
        u_qp = self.get(state, 'val').reshape(qp_coors.shape)

        # Deformed QP coordinates.
        coors = u_qp + qp_coors

        force = nm.zeros(coors.shape[0], dtype=nm.float64)
        normals = nm.zeros((coors.shape[0], 3), dtype=nm.float64)
        fd = None

        k, f0, eps, a = self._get_force_pars(force_pars, force.shape)

        # Active points in contact change with displacements!
        ii = self.cs.mask_points(coors, nm.abs(eps.min()))

        if diff_var is None:
            if ii.any():
                d, normals[ii] = self.cs.get_distance(coors[ii])
                force[ii] = self.smooth_f(d, k[ii], f0[ii], a[ii], eps[ii], 0)

            fmode = 0

        else:
            if ii.any():
                d, normals[ii] = self.cs.get_distance(coors[ii])
                # Force derivative.
                force[ii] = self.smooth_f(d, k[ii], f0[ii], a[ii], eps[ii], 1)

                # Force / (R - d).
                aux = self.smooth_f(d, k[ii], f0[ii], a[ii], eps[ii], 0)
                fd = nm.zeros_like(force)
                fd[ii] = aux / (self.cs.radius - d)
                fd.shape = qps.shape[:2] + (1, 1)

            fmode = 1

        force.shape = qps.shape[:2] + (1, 1)
        normals.shape = qps.shape[:2] + (3, 1)

        return force, normals, fd, sg, fmode

class SufaceNormalDotTerm(Term):
    r"""
    "Scalar traction" term, (weak form).

    :Definition:

    .. math::
        \int_{\Gamma} q \ul{c} \cdot \ul{n}

    :Arguments:
        - material : :math:`\ul{c}`
        - virtual  : :math:`q`
    """
    name = 'dw_surface_ndot'
    arg_types = (('material', 'virtual'),
                 ('material', 'parameter'))
    arg_shapes = {'material' : 'D, 1', 'virtual' : (1, None), 'parameter' : 1}
    modes = ('weak', 'eval')
    integration = 'facet'

    @staticmethod
    def dw_fun(out, material, bf, sg):
        bf_t = nm.tile(bf.transpose((0, 1, 3, 2)), (out.shape[0], 1, 1, 1))
        bf_t = nm.ascontiguousarray(bf_t)
        aux = dot_sequences(material, sg.normal, 'ATB')
        status = sg.integrate(out, bf_t * aux)
        return status

    @staticmethod
    def d_fun(out, material, val, sg):
        aux = dot_sequences(material, sg.normal, 'ATB')
        status = sg.integrate(out, val * aux)
        return status

    def get_fargs(self, mat, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        if mode == 'weak':
            return mat, sg.bf, sg

        elif mode == 'eval':
            val = self.get(virtual, 'val')
            return mat, val, sg

        else:
            raise ValueError('unsupported evaluation mode in %s! (%s)'
                             % (self.name, mode))

    def get_eval_shape(self, mat, virtual,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(virtual)

        return (n_el, 1, 1, 1), virtual.dtype

    def set_arg_types( self ):
        if self.mode == 'weak':
            self.function = self.dw_fun

        else:
            self.function = self.d_fun

class SDSufaceIntegrateTerm(Term):
    r"""
    Sensitivity of scalar traction.

    :Definition:

    .. math::
        \int_{\Gamma} p \nabla \cdot \ul{\Vcal}

    :Arguments:
        - parameter : :math:`p`
        - parameter_mv : :math:`\ul{\Vcal}`
    """
    name = 'ev_sd_surface_integrate'
    arg_types = ('parameter', 'parameter_mv')
    arg_shapes = {'parameter': 1, 'parameter_mv': 'D'}
    integration = 'facet'

    @staticmethod
    def function(out, val_p, div_v, sg):
        status = sg.integrate(out, val_p * div_v)
        return status

    def get_fargs(self, par, par_v,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(par)

        val_p = self.get(par, 'val')
        div_v = self.get(par_v, 'div', integration='facet_extra')
        return val_p, div_v, sg

    def get_eval_shape(self, par, par_v,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_v)

        return (n_el, 1, 1, 1), par.dtype

class SurfaceJumpTerm(Term):
    r"""
    Interface jump condition.

    :Definition:

    .. math::
        \int_{\Gamma} c\, q (p_1 - p_2)

    :Arguments:
        - material : :math:`c`
        - virtual  : :math:`q`
        - state_1  : :math:`p_1`
        - state_2  : :math:`p_2`
    """
    name = 'dw_jump'
    arg_types = ('opt_material', 'virtual', 'state_1', 'state_2')
    arg_shapes = [{'opt_material' : '1, 1', 'virtual' : (1, None),
                   'state_1' : 1, 'state_2' : 1},
                  {'opt_material' : None}]
    integration = 'facet'

    @staticmethod
    def function(out, jump, mul, bf1, bf2, sg, fmode):
        bf_t = nm.tile(sg.bf.transpose((0, 1, 3, 2)), (out.shape[0], 1, 1, 1))

        if fmode == 0:
            vec = bf_t * jump

        elif fmode == 1:
            vec = bf_t * bf1

        else:
            vec = - bf_t * bf2

        if mul is None:
            status = sg.integrate(out, vec)

        else:
            status = sg.integrate(out, mul * vec)

        return status

    def get_fargs(self, coef, virtual, state1, state2,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)
        sg1, _ = self.get_mapping(state1)
        sg2, _ = self.get_mapping(state2)

        if diff_var is None:
            val1 = self.get(state1, 'val')
            val2 = self.get(state2, 'val')
            jump = val1 - val2
            fmode = 0

        else:
            jump = None

            if diff_var == state1.name:
                fmode = 1

            else:
                fmode = 2

        return jump, coef, sg1.bf, sg2.bf, sg, fmode
