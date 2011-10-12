import numpy as nm

from sfepy.base.base import assert_
from sfepy.linalg import dot_sequences
from sfepy.mechanics.tensors import dim2sym, transform_data
import sfepy.mechanics.membranes as membranes
from sfepy.fem.poly_spaces import PolySpace
from sfepy.terms.terms import Term

def eval_membrane_mooney_rivlin(a1, a2, mtx_c, c33, mode):
    """
    Evaluate stress or tangent stiffness of the Mooney-Rivlin membrane.

    [1] Baoguo Wu, Xingwen Du and Huifeng Tan: A three-dimensional FE
    nonlinear analysis of membranes, Computers & Structures 59 (1996),
    no. 4, 601--605.
    """
    a12 = 2.0 * a1[..., 0, 0]
    a22 = 2.0 * a2[..., 0, 0]

    sh = mtx_c.shape
    sym = dim2sym(sh[2])

    c11 = mtx_c[..., 0, 0]
    c12 = mtx_c[..., 0, 1]
    c22 = mtx_c[..., 1, 1]
    pressure = c33 * (a12 + a22 * (c11 + c22))

    if mode == 0:
        out = nm.empty((sh[0], sh[1], sym, 1))

        # S_11, S_22, S_12.
        out[..., 0, 0] = -pressure * c22 * c33 + a12 + a22 * (c22 + c33)
        out[..., 1, 0] = -pressure * c11 * c33 + a12 + a22 * (c11 + c33)
        out[..., 2, 0] = +pressure * c12 * c33 - a22 * c12

    else:
        out = nm.empty((sh[0], sh[1], sym, sym))

        dp11 = a22 * c33 - pressure * c22 * c33
        dp22 = a22 * c33 - pressure * c11 * c33
        dp12 = 2.0 * pressure * c12 * c33

        # D_11, D_22, D_33
        out[..., 0, 0] = - 2.0 * ((a22 - pressure * c22) * c22 * c33**2
                                  + c33 * c22 * dp11)
        out[..., 1, 1] = - 2.0 * ((a22 - pressure * c11) * c11 * c33**2
                                  + c33 * c11 * dp22)
        out[..., 2, 2] = - a22 + pressure * (c33 + 2.0 * c12**2 * c33**2) \
                         + c12 * c33 * dp12

        # D_21, D_31, D_32
        out[..., 1, 0] = 2.0 * ((a22 - pressure * c33
                                 - (a22 - pressure * c11) * c22 * c33**2)
                                - c33 * c11 * dp11)
        out[..., 2, 0] = 2.0 * (-pressure * c12 * c22 * c33**2
                                + c12 * c33 * dp11)
        out[..., 2, 1] = 2.0 * (-pressure * c12 * c11 * c33**2
                                + c12 * c33 * dp22)

        out[..., 0, 1] = out[..., 1, 0]
        out[..., 0, 2] = out[..., 2, 0]
        out[..., 1, 2] = out[..., 2, 1]

        # D_12, D_13, D_23
        ## out[..., 0, 1] = 2.0 * ((a22 - pressure * c33
        ##                          - (a22 - pressure * c22) * c11 * c33**2)
        ##                         - c33 * c22 * dp22)
        ## out[..., 0, 2] = 2.0 * (a22 - pressure * c22) * c12 * c33**2 \
        ##                  - c33 * c22 * dp12
        ## out[..., 1, 2] = 2.0 * (a22 - pressure * c11) * c12 * c33**2 \
        ##                  - c33 * c11 * dp12

    return out

class TLMembraneTerm(Term):
    r"""
    Mooney-Rivlin membrane with plain stress assumption.

    :Arguments:
        virtual  : :math:`\ul{v}`,
        state    : :math:`\ul{u}`
    """
    name = 'dw_tl_membrane'
    arg_types = ('material_a1', 'material_a2', 'material_h0',
                 'virtual', 'state')
    integration = 'surface'

    @staticmethod
    def function(out, fun, *args):
        """
        Notes
        -----
        `fun` is either `weak_function` or `eval_function` according to
        evaluation mode.
        """
        return fun(out, *args)

    @staticmethod
    def weak_function(out, a1, a2, h0, mtx_c, c33, mtx_b, mtx_t, bfg, geo,
                      fmode):
        crt = eval_membrane_mooney_rivlin(a1, a2, mtx_c, c33, fmode)

        n_ep = bfg.shape[3]

        if fmode == 0:
            bts = dot_sequences(mtx_b, crt, 'ATB')

            status = geo.integrate(out, bts * h0)

            # Transform to global coordinate system, one node at
            # a time.
            for iep in range(n_ep):
                ir = slice(iep, None, n_ep)
                fn = out[:, 0, ir, 0]
                fn[:] = dot_sequences(mtx_t, fn, 'AB')


        else:
            btd = dot_sequences(mtx_b, crt, 'ATB')
            btdb = dot_sequences(btd, mtx_b)

            stress = eval_membrane_mooney_rivlin(a1, a2, mtx_c, c33, 0)

            kts =  membranes.get_tangent_stress_matrix(stress, bfg)

            mtx_k = kts + btdb

            status = geo.integrate(out, mtx_k * h0)

            # Transform to global coordinate system, one node at
            # a time.
            dot = dot_sequences
            for iepr in range(n_ep):
                ir = slice(iepr, None, n_ep)
                for iepc in range(n_ep):
                    ic = slice(iepc, None, n_ep)
                    fn = out[:, 0, ir, ic]
                    fn[:] = dot(dot(mtx_t, fn, 'AB'), mtx_t, 'ABT')

        return status

    @staticmethod
    def eval_function(out, a1, a2, h0, mtx_c, c33, mtx_b, mtx_t, geo,
                      term_mode, fmode):

        if term_mode == 'strain':
            out_qp = membranes.get_green_strain_sym3d(mtx_c, c33)

        elif term_mode == 'stress':
            n_el, n_qp, dm, _ = mtx_c.shape
            dim = dm + 1
            sym = dim2sym(dim)
            out_qp = nm.zeros((n_el, n_qp, sym, 1), dtype=mtx_c.dtype)

            stress = eval_membrane_mooney_rivlin(a1, a2, mtx_c, c33, 0)
            out_qp[..., 0:2, 0] = stress[..., 0:2, 0]
            out_qp[..., 3, 0] = stress[..., 2, 0]

        status = geo.integrate(out, out_qp, fmode)
        out[:] = transform_data(out, mtx=mtx_t)

        return status

    def __init__(self, *args, **kwargs):
        Term.__init__(self, *args, **kwargs)

        igs = self.region.igs
        self.mtx_t = {}.fromkeys(igs, None)
        self.membrane_geo = {}.fromkeys(igs, None)
        self.bfg = {}.fromkeys(igs, None)

    def describe_membrane_geometry(self, ig, field, sg, sd):
        # Coordinates of element vertices.
        coors = field.coors[sd.econn[:, :sg.nFP]]

        # Coordinate transformation matrix (transposed!).
        self.mtx_t[ig] = membranes.create_transformation_matrix(coors)

        # Transform coordinates to the local coordinate system.
        coors_loc = dot_sequences((coors - coors[:, 0:1, :]), self.mtx_t[ig])

        # Mapping from transformed element to reference element.
        gel = field.gel.surface_facet
        vm = membranes.create_mapping(coors_loc, gel, 1)

        qp = self.integral.get_qp(gel.name)
        ps = PolySpace.any_from_args(None, gel, field.approx_order)
        self.membrane_geo[ig] = vm.get_mapping(*qp, poly_space=ps)

        # Transformed base function gradient w.r.t. material coordinates
        # in quadrature points.
        self.bfg[ig] = self.membrane_geo[ig].bfg

    def get_fargs(self, a1, a2, h0, virtual, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vv, vu = virtual, state

        ap, sg = self.get_approximation(vv)
        sd = ap.surface_data[self.region.name]

        ig = self.char_fun.ig
        if self.mtx_t[ig] is None:
            self.describe_membrane_geometry(ig, vu.field, sg, sd)

        mtx_t = self.mtx_t[ig]
        bfg = self.bfg[ig]
        geo = self.membrane_geo[ig]

        # Displacements of element nodes.
        vec_u = vu.get_state_in_region(self.region, igs=[self.char_fun.ig])
        el_u = vec_u[sd.leconn]

        # Transform displacements to the local coordinate system.
        # u_new = T^T u
        el_u_loc = dot_sequences(el_u, mtx_t, 'AB')
        ## print el_u_loc

        mtx_c, c33, mtx_b = membranes.describe_deformation(el_u_loc, bfg)

        if mode == 'weak':
            fmode = diff_var is not None

            return (self.weak_function,
                    a1, a2, h0, mtx_c, c33, mtx_b, mtx_t, bfg, geo, fmode)

        else:
            fmode = {'eval' : 0, 'el_avg' : 1, 'qp' : 2}.get(mode, 1)

            assert_(term_mode in ['strain', 'stress'])

            return (self.eval_function,
                    a1, a2, h0, mtx_c, c33, mtx_b, mtx_t, geo, term_mode, fmode)

    def get_eval_shape(self, a1, a2, h0, virtual, state,
                       mode=None, term_mode=None, diff_var=None, **kwargs):
        n_el, n_qp, dim, n_en, n_c = self.get_data_shape(state)
        sym = dim * (dim + 1) / 2

        return (n_el, 1, sym, 1), state.dtype
