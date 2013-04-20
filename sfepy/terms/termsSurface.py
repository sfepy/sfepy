import numpy as nm

from sfepy.terms.terms import Term, terms
from sfepy.linalg import dot_sequences
from sfepy.mechanics.contact_planes import ContactPlane

##
# 22.08.2006, c
class LinearTractionTerm( Term ):
    r"""
    Linear traction forces (weak form), where, depending on dimension of
    'material' argument, :math:`\ull{\sigma} \cdot \ul{n}` is
    :math:`\bar{p} \ull{I} \cdot \ul{n}` for a given scalar pressure,
    :math:`\ul{f}` for a traction vector, and itself for a stress tensor.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot \ull{\sigma} \cdot \ul{n}

    :Arguments:
        - material : :math:`\ull{\sigma}`
        - virtual  : :math:`\ul{v}`
    """
    name = 'dw_surface_ltr'
    arg_types = ('material', 'virtual')
    arg_shapes = [{'material' : 'S, 1', 'virtual' : ('D', None)},
                  {'material' : 'D, 1'}, {'material' : '1, 1'}]
    integration = 'surface'

    function = staticmethod(terms.dw_surface_ltr)

    def get_fargs(self, traction, virtual,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        return traction, sg

class ContactPlaneTerm(Term):
    r"""
    Small deformation contact plane term with linear penetration penalty.

    The plane is given by an anchor point :math:`\ul{A}` and a normal
    :math:`\ul{n}`. The contact occurs in points that orthogonally project onto
    the plane into a polygon given by orthogonal projections of boundary points
    :math:`\{\ul{B}_i\}`, :math:`i = 1, \dots, N_B` on the plane. In such
    points, a penetration distance :math:`d(\ul{u}) = (\ul{X} + \ul{u} -
    \ul{A}, \ul{n})` is computed, and if it is positive, a force :math:`k
    d(\ul{u})^+ \ul{n}` is applied.

    :Definition:

    .. math::
        \int_{\Gamma} \ul{v} \cdot k d(\ul{u})^+ \ul{n}

    :Arguments:
        - material_k : :math:`k`
        - material_n : :math:`\ul{n}`
        - material_a : :math:`\ul{A}`
        - material_b : :math:`\{\ul{B}_i\}`, :math:`i = 1, \dots, N_B`
        - virtual    : :math:`\ul{v}`
        - state      : :math:`\ul{u}`
    """
    name = 'dw_contact_plane'
    arg_types = ('material_k', 'material_n', 'material_a', 'material_b',
                 'virtual', 'state')
    integration = 'surface'

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

    def get_fargs(self, stiffness, normal, anchor, bounds, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(virtual)

        if self.cp is None:
            self.cp = ContactPlane(anchor, normal, bounds)

        ig = self.char_fun.ig
        qps = self.get_physical_qps()
        qp_coors = qps.values[ig]
        u_qp = self.get(state, 'val').reshape(qp_coors.shape)

        # Deformed QP coordinates.
        coors = u_qp + qp_coors

        force = nm.zeros(coors.shape[0], dtype=nm.float64)

        # Active points in contact change with displacements!
        ii = self.cp.mask_points(coors)

        if diff_var is None:
            if ii.any():
                dist = self.cp.get_distance(coors[ii])
                # Force in the plane normal direction.
                force[ii] = nm.where(dist > 0.0, stiffness * dist, 0.0)

            fmode = 0

        else:
            if ii.any():
                dist = self.cp.get_distance(coors[ii])
                # Force in the plane normal direction derivative.
                force[ii] = nm.where(dist > 0.0, stiffness, 0.0)

            fmode = 1

        force.shape = qps.shape[ig][:2] + (1, 1)

        return force, self.cp.normal, sg, fmode

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
    integration = 'surface'

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

class SDSufaceNormalDotTerm(Term):
    r"""
    Sensitivity of scalar traction.

    :Definition:

    .. math::
        \int_{\Gamma} p \ul{c} \cdot \ul{n} \nabla \cdot \ul{\Vcal}

    :Arguments:
        - material : :math:`\ul{c}`
        - parameter : :math:`p`
        - parameter_mesh_velocity : :math:`\ul{\Vcal}`
    """
    name = 'd_sd_surface_ndot'
    arg_types = ('material', 'parameter', 'parameter_mesh_velocity')
    arg_shapes = {'material' : 'D, 1', 'parameter' : 1,
                  'parameter_mesh_velocity' : 'D'}
    integration = 'surface'

    @staticmethod
    def function(out, material, val_p, div_mv, sg):
        aux = dot_sequences(material, sg.normal, 'ATB')
        aux2 = dot_sequences(aux, div_mv)
        status = sg.integrate(out, val_p * aux2)
        return status

    def get_fargs(self, mat, par, par_mv,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        sg, _ = self.get_mapping(par)

        val_p = self.get(par, 'val')
        div_mv = self.get(par_mv, 'div', integration='surface_extra')
        return mat, val_p, div_mv, sg

    def get_eval_shape(self, mat, par, par_mv,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(par_mv)

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
    integration = 'surface'

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
