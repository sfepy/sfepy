import numpy as nm
from sfepy.terms.terms import Term

from dg_field import get_unraveler, get_raveler


def unravel_sol(state):
    """
    Unpacks solution from flat vector to
    (n_cell, n_el_nod, n_comp)


    :param state:
    """
    u = state.data[0]
    n_cell = state.field.n_cell
    n_el_nod = state.field.poly_space.n_nod
    unravel = get_unraveler(n_el_nod, n_cell)
    return unravel(u)


class AdvVolDGTerm(Term):
    # Can be removed use  DotProductVolumeTerm from sfepy.terms.terms_dot

    name = "dw_dg_volume"
    modes = ("weak",)
    arg_types = ('virtual', 'state')
    # arg_types = ('ts', 'virtual', 'state')
    arg_shapes = {'virtual': 1, #('D', 'state'),
                  'state': 1}
    symbolic = {'expression': '1',
                'map': {'u': 'state'}
                }

    def __init__(self, integral, region, u=None, v=None):
        Term.__init__(self, "adv_vol(v, u)", "v, u", integral, region, u=u, v=v)
        self.u = u
        self.v = v
        self.setup()

    def get_fargs(self, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        if diff_var is not None:
            mtx_mode = True
            u = unravel_sol(state)
        else:
            mtx_mode = False
            u = unravel_sol(state)

        dim = state.field.dim
        n_el_nod = state.field.poly_space.n_nod
        vols = self.region.domain.cmesh.get_volumes(dim)[:, None]
        return u, mtx_mode, n_el_nod, vols

    def function(self, out, u, mtx_mode, n_el_nod, vols):
        if mtx_mode:
            # integral over element with constant test
            # function is just volume of the element
            out[:] = 0.0
            # out[:, 0, 0, 0] = vols
            # out[:, 0, 1, 1] = vols / 3.0
            for i in range(n_el_nod):
                # these are values for legendre basis in 1D!
                out[:, :, i, i] = vols / (2.0 * i + 1.0)
        else:
            out[:] = 0.0
            for i in range(n_el_nod):
                out[:, :, i, 0] = vols / (2.0 * i + 1.0) * u[:, i]
        status = None
        return status


class AdvFluxDGTerm1D(Term):

    def __init__(self, integral, region, u=None, v=None, a=lambda x: 1):
        Term.__init__(self, "adv_lf_flux(a.val, v, u)", "a.val, v, u", integral, region, u=u, v=v, a=a)
        self.u = u
        self.v = v
        self.a = a
        self.setup()

    name = "dw_dg_advect_flux"
    modes = ("weak",)
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material': 'D, 1',
                  'virtual': ('D', 'state'),
                  'state'   : 1}
    symbolic = {'expression' : 'grad(a*u)',
                'map': {'u': 'state', 'a': 'material'}
    }

    def get_fargs(self, a, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        if diff_var is not None:
            # do not eval in matrix mode, we however still need
            # this term to have diff_var in order for it to receive the values
            doeval = False
            return None, None, None, None, doeval, 0, 0
        else:
            doeval = True

            u = unravel_sol(state)
            n_el_nod = state.field.poly_space.n_nod
            n_el_facets = state.field.n_el_facets
            nb_dofs, nb_normals = state.field.get_nbrhd_dofs(state.field.region, state)
            # state variable has dt in it!

            fargs = (u, nb_dofs, nb_normals, a[:, :1, 0, 0], doeval, n_el_nod, n_el_facets)
            return fargs

    def function(self, out, u, nb_u, nb_n, velo, doeval, n_el_nod, n_el_facets):
        if not doeval:
            out[:] = 0.0
            return None

        #  the Lax-Friedrichs flux is

        #       F(a, b) = 1/2(f(a) + f(b)) + max(|f'(w)|) / 2 * (a - b)

        # in our case a and b are values from elements left and right of
        # the respective element boundaries
        # for Legendre basis these are:
        # u_left = U_0 + U_1 + U_2 + ...
        # u_right = U_0 - U_1 + U_2 - U_3 ... = sum_{p=0}^{order} (-1)^p * U_p

        # left flux is calculated in j_-1/2  where U(j-1) and U(j) meet
        # right flux is calculated in j_+1/2 where U(j) and U(j+1) meet

        # fl:
        # fl = velo[:, 0] * (ul[:, 0] + ul[:, 1] +
        #                    (u[:, 0] - u[:, 1])) / 2 + \
        #      nm.abs(velo[:, 0]) * (ul[:, 0] + ul[:, 1] -
        #                            (u[:, 0] - u[:, 1])) / 2
        a = 0
        b = 0
        sign = 1
        for i in range(n_el_nod):
            a += nb_u[:, 0, i]  # integral left
            b += sign * u[:, i]
            sign *= -1

        fl = (velo * a + velo * b) / 2 + \
              nm.abs(velo) * (a - b) / 2

        # fl:
        # fp = velo[:, 0] * (u[:, 0] + u[:, 1] +
        #                            (ur[:, 0] - ur[: , 1])) / 2 + \
        #              nm.abs(velo[:, 0]) * (u[:, 0] + u[:, 1] -
        #                                    (ur[:, 0] - ur[:, 1])) / 2
        a = 0
        b = 0
        sign = 1
        for i in range(n_el_nod):
            a += u[:, i]
            b += sign * nb_u[:, 1 , i]
            sign *= -1

        fp = (velo * a + velo * b) / 2 + \
              nm.abs(velo) * (a - b) / 2

        out[:] = 0.0

        # flux0 = (fl - fp)
        # flux1 = (- fl - fp + intg1)
        # out[:, 0, 0, 0] = -flux0
        # out[:, 0, 1, 0] = -flux1

        # for Legendre basis integral of higher order
        # functions of the basis is zero,
        # hence we calculate integral
        #
        # int_{j-1/2}^{j+1/2} f(u)dx
        #
        # only from the zero order function, over [-1, 1] - hence the 2
        intg1 = velo * u[:, 0] * 2
        intg2 = velo * u[:, 1] * 2 if n_el_nod > 2 else 0
        # i.e. intg1 = a * u0 * reference_el_vol

        flux = list()
        flux[:3] = (fl - fp), (- fl - fp + intg1), (fl - fp + intg2)

        for i in range(n_el_nod):
            out[:, :, i, 0] = -flux[i]


        status = None
        return status


class AdvFluxDGTerm(Term):

    def __init__(self, integral, region, u=None, v=None, a=lambda x: 1):
        Term.__init__(self, "adv_lf_flux(a.val, v, u)", "a.val, v, u", integral, region, u=u, v=v, a=a)
        self.u = u
        self.v = v
        self.a = a
        self.setup()

    name = "dw_dg_advect_flux"
    modes = ("weak",)
    arg_types = ('material', 'virtual', 'state')
    arg_shapes = {'material': 'D, 1',
                  'virtual': 1,  #('D', 'state'),
                  'state': 1}
    symbolic = {'expression' : 'grad(a*u)',
                'map': {'u': 'state', 'a': 'material'}
    }

    def get_fargs(self, a, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        if diff_var is not None:
            # do not eval in matrix mode, we however still need
            # this term to have diff_var in order for it to receive the values
            doeval = False
            dofs = unravel_sol(state)
            return dofs, None, None, None, a[:, 0, :, 0], doeval, 0, 0
        else:
            doeval = True

            field = state.field
            region = field.region

            dofs = unravel_sol(state)
            cell_normals = field.get_cell_normals_per_facet(region)
            facet_base_vals = field.get_facet_base(base_only=True)
            inner_facet_qp_vals, outer_facet_qp_vals, whs = field.get_both_facet_qp_vals(dofs, region)

            fargs = (dofs, inner_facet_qp_vals, outer_facet_qp_vals, facet_base_vals, whs,  cell_normals, a[:, 0, :, 0], doeval)
            return fargs

    def function(self, out, dofs, in_fc_v, out_fc_v, fc_b, whs, fc_n, velo, doeval):
        if not doeval:
            out[:] = 0.0
            return None

        alf = 0

        n_cell = dofs.shape[0]
        n_el_nod = dofs.shape[1]
        n_el_facets = fc_n.shape[-2]
        #  Calculate integrals over facets representing Lax-Friedrichs fluxes
        C = nm.abs(nm.sum(fc_n * velo[:, None, :], axis=-1))[:, None]
        facet_fluxes = nm.zeros((n_cell, n_el_facets, n_el_nod))
        for facet_n in range(n_el_facets):
            for n in range(n_el_nod):
                fc_v_p = in_fc_v[:, facet_n, :] + out_fc_v[:, facet_n, :]
                fc_v_m = in_fc_v[:, facet_n, :] - out_fc_v[:, facet_n, :]
                central = velo[:, None, :] * fc_v_p[:, :, None]/2.
                upwind = (C[:, :, facet_n]*(1 - alf)/2. * fc_n[:, facet_n])[..., None, :] * fc_v_m[:, :, None]
                facet_fluxes[:, facet_n, n] = nm.sum(fc_n[:, facet_n] *
                                                             nm.sum((central + upwind) *
                                                                    (fc_b[None, :, 0, facet_n, 0, n] *
                                                                     whs[:, facet_n, :])[..., None], axis=1),
                                                             axis=1)
        cell_fluxes = nm.sum(facet_fluxes, axis=1)

        if False:
            from my_utils.visualizer import plot_1D_legendre_dofs, reconstruct_legendre_dofs
            import matplotlib.pyplot as plt
            x = self.region.domain.mesh.coors
            plot_1D_legendre_dofs(x, (dofs,))
            ww, xx = reconstruct_legendre_dofs(x, 1, dofs.T[0, ..., None, None])
            plt.plot(xx, ww[:, 0], label="recon")

            # plt.figure("Vals outer")
            # plt.plot(x[:-1], out_fc_v[:, 0], marker=".", label="left outer", color="b",  ls="")
            # plt.plot(x[1:], out_fc_v[:, 1], marker=".",  label="right outer", color="r",  ls="")
            # plt.plot(xx, ww[:, 0], label="recon")
            # plt.legend()
            #
            # # Plot mesh
            # XN = x[-1]
            # X1 = x[0]
            # Xvol = XN - X1
            # X = (x[1:] + x[:-1]) / 2
            # plt.vlines(x[:, 0], ymin=0, ymax=.5, colors="grey")
            # plt.vlines((X1, XN), ymin=0, ymax=.5, colors="k")
            # plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")
            #
            # plt.figure("Vals inner")
            # plt.plot(xx, ww[:, 0], label="recon")
            # plt.plot(x[:-1], in_fc_v[:, 0], marker=".", label="left in", color="b", ls="")
            # plt.plot(x[1:], in_fc_v[:, 1], marker=".", label="right in", color="r", ls="")
            # plt.legend()
            #
            # # Plot mesh
            # XN = x[-1]
            # X1 = x[0]
            # Xvol = XN - X1
            # X = (x[1:] + x[:-1]) / 2
            # plt.vlines(x[:, 0], ymin=0, ymax=.5, colors="grey")
            # plt.vlines((X1, XN), ymin=0, ymax=.5, colors="k")
            # plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")
            # plt.legend()

            for n in range(n_el_nod):
                plt.figure("Flux {}".format(n))
                fig = plt.gcf()
                fig.clear()
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(x[:-1], facet_fluxes[:, 0, n], marker=">", label="flux {} left".format(n), color="b", ls="")
                plt.plot(x[1:], facet_fluxes[:, 1,  n], marker=">", label="flux {} right".format(n), color="r", ls="")


                # Plot mesh
                XN = x[-1]
                X1 = x[0]
                Xvol = XN - X1
                X = (x[1:] + x[:-1]) / 2
                plt.vlines(x[:, 0], ymin=0, ymax=.5, colors="grey")
                plt.vlines((X1, XN), ymin=0, ymax=.5, colors="k")
                plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")
                plt.legend()

            # plt.show()

            for n in range(n_el_nod):
                plt.figure("Flux {}".format(n))
                fig = plt.gcf()
                fig.clear()
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(X, cell_fluxes[:, n], marker="D", label="cell flux {}".format(n), color="r", ls="")


                # Plot mesh
                XN = x[-1]
                X1 = x[0]
                Xvol = XN - X1
                X = (x[1:] + x[:-1]) / 2
                plt.vlines(x[:, 0], ymin=0, ymax=.5, colors="grey")
                plt.vlines((X1, XN), ymin=0, ymax=.5, colors="k")
                plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")
                plt.legend()

        out[:] = 0.0
        for i in range(n_el_nod):
            out[:, :, i, 0] = cell_fluxes[:, i, None]


        status = None
        return status

class ScalarDotMGradScalarDGTerm(Term):
    r"""
    Volume dot product of a scalar gradient dotted with a material vector with
    a scalar.

    :Definition:

    .. math::
        \int_{\Omega} q \ul{y} \cdot \nabla p \mbox{ , }
        \int_{\Omega} p \ul{y} \cdot \nabla q

    :Arguments 1:
        - material : :math:`\ul{y}`
        - virtual  : :math:`q`
        - state    : :math:`p`

    :Arguments 2:
        - material : :math:`\ul{y}`
        - state    : :math:`p`
        - virtual  : :math:`q`
    """
    name = 'dw_s_dot_mgrad_s'
    arg_types = (('material', 'virtual', 'state'),
                 ('material', 'state', 'virtual'))
    arg_shapes = [{'material' : 'D, 1',
                   'virtual/grad_state' : (1, None),
                   'state/grad_state' : 1,
                   'virtual/grad_virtual' : (1, None),
                   'state/grad_virtual' : 1}]
    modes = ('grad_state', 'grad_virtual')

    @staticmethod
    def function(out, out_qp, geo, fmode):
        if fmode == 0:
            status = geo.integrate(out, out_qp)
        else:
            status = None
            out[:] = 0.0
        return status

    def get_fargs(self, mat, var1, var2,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        from sfepy.linalg import dot_sequences

        vg1, _ = self.get_mapping(var1)
        vg2, _ = self.get_mapping(var2)

        if diff_var is None:
            if self.mode == 'grad_state':
                geo = vg1
                bf_t = vg1.bf.transpose((0, 1, 3, 2))
                val_qp = self.get(var2, 'grad')
                out_qp = bf_t * dot_sequences(mat, val_qp, 'ATB')

            else:
                geo = vg2
                val_qp = self.get(var1, 'val')
                out_qp = dot_sequences(vg2.bfg, mat, 'ATB') * val_qp

            fmode = 0

        else:
            if self.mode == 'grad_state':
                geo = vg1
                bf_t = vg1.bf.transpose((0, 1, 3, 2))
                out_qp = bf_t * dot_sequences(mat, vg2.bfg, 'ATB')

            else:
                geo = vg2
                out_qp = dot_sequences(vg2.bfg, mat, 'ATB') * vg1.bf

            fmode = 1

        return out_qp, geo, fmode