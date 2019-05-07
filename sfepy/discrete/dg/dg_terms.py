import numpy as nm


# sfepy imports
from sfepy.terms.terms import Term, terms
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)
from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm

from sfepy.discrete.dg.dg_field import get_unraveler, get_raveler


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

class AdvectDGFluxTerm(Term):
    r"""
    Lax-Friedrichs flux term for advection of scalar quantity :math:`p` with the advection velocity
    :math:`\ul{a}` given as a material parameter (a known function of space and time).

    Part of the discretization of
        math::



    Is residual only! Use with state variable history [-1].

    :Definition:

    .. math::

        \int_{\partial{T_K}} \vec{n} \cdot f^{*} (p_{in}, p_{out})\cdot q

        where
            f^{*}(p_{in}, p_{out}) =  \vec{a}  \frac{p_{in} + p_{out}}{2}  + (1 - \alpha) \vec{n}
            C \frac{ p_{in} - p_{out}}{2},

        $\alpha \in [0, 1]$, $\alpha = 0$ for upwind scheme, $\alpha = 1$ for central scheme,  and

            C = \max_{u \in [?, ?]} \abs{n_x \pdiff{a_1}{u} + n_y \pdiff{a_2}{u}}

        the $p_{in}$ resp. $p_{out}$ is solution on the boundary of the element provided
        by element itself resp. its neighbor and a is advection velocity.

    :Arguments 1:
        - material : :math:`\ul{a}`
        - virtual  : :math:`q`
        - state    : :math:`p`

    :Arguments 3:
        - material    : :math:`\ul{a}`
        - virtual     : :math:`q`
        - state       : :math:`p`
        - opt_material : :math: `\alpha`
    """

    alf = 0
    name = "dw_dg_advect_laxfrie_flux"
    modes = ("weak",)
    arg_types = ('opt_material', 'material_advelo', 'virtual', 'state')
    arg_shapes = [{'opt_material': '.: 1',
                  'material_advelo': 'D, 1',
                  'virtual': (1, 'state'),
                  'state': 1
                  },
                {'opt_material': None}]
    integration = 'volume'
    symbolic = {'expression' : 'div(a*u)*w',
                'map': {'u': 'state', 'a': 'material', 'v' : 'virtual'}
    }

    def get_fargs(self, alpha, advelo, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        if alpha is not None:
            self.alf = alpha # extract alpha value regardless of shape

        if diff_var is not None:
            output("Diff var is not None in residual only term {} ! Skipping.".format(self.name))
            return None, None, None, None, None, advelo[:, 0, :, 0], 0
        else:
            field = state.field
            region = field.region

            if not "DG" in field.family_name:
                raise ValueError("Used DG term with non DG field {} of family {}".format(field.name, field.family_name))

            shapes = (field.n_cell, field.n_el_nod)
            cell_normals = field.get_cell_normals_per_facet(region)
            facet_base_vals = field.get_facet_base(base_only=True)
            inner_facet_qp_vals, outer_facet_qp_vals, weights = field.get_both_facet_qp_vals(state, region)

            self.dofs = unravel_sol(state)

            fargs = (shapes, inner_facet_qp_vals, outer_facet_qp_vals,
                     facet_base_vals, weights, cell_normals, advelo[:, 0, :, 0])
            return fargs

    # noinspection PyUnreachableCode
    def function(self, out, shapes, in_fc_v, out_fc_v, fc_b, whs, fc_n, advelo):
        """

        :param out:
        :param shapes: (n_cell, n_el_nod)
        :param in_fc_v: inner values for facets per cell, shape = (n_cell, n_el_faces, 1)
        :param out_fc_v: outer values for facets per cell, shape = (n_cell, n_el_faces, 1)
        :param fc_b: values of basis in facets qps, shape = (1, n_el_facet, n_qp)
        :param whs: weights of the qps on facets, shape = (n_cell, n_el_facet, n_qp
        :param fc_n: facet normals, shape = (n_cell, n_el_facets)
        :param advelo: advection velocity, shape = (n_cell, 1)
        :return:
        """
        if shapes is None:
            out[:] = 0.0
            return None

        n_cell = shapes[0]
        n_el_nod = shapes[1]
        n_el_facets = fc_n.shape[-2]
        #  Calculate integrals over facets representing Lax-Friedrichs fluxes
        C = nm.abs(nm.sum(fc_n * advelo[:, None, :], axis=-1))[:, None]
        facet_fluxes = nm.zeros((n_cell, n_el_facets, n_el_nod))
        for facet_n in range(n_el_facets):
            for mode_n in range(n_el_nod):
                fc_v_p = in_fc_v[:, facet_n, :] + out_fc_v[:, facet_n, :]
                fc_v_m = in_fc_v[:, facet_n, :] - out_fc_v[:, facet_n, :]
                central = advelo[:, None, :] * fc_v_p[:, :, None] / 2.
                upwind = ((1 - self.alf)/2. * C[:, :, facet_n] * fc_n[:, facet_n])[..., None, :] * fc_v_m[:, :, None]

                facet_fluxes[:, facet_n, mode_n] = nm.sum(fc_n[:, facet_n] *
                                                     nm.sum((central + upwind) *
                                                            (fc_b[None, :, 0, facet_n, 0, mode_n] *
                                                             whs[:, facet_n, :])[..., None], axis=1),
                                                     axis=1)

        cell_fluxes = nm.sum(facet_fluxes, axis=1)

        # 1D plots
        if False:
            dofs = self.dofs
            # TODO remove 1D plots
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

            for mode_n in range(n_el_nod):
                plt.figure("Flux {}".format(mode_n))
                fig = plt.gcf()
                fig.clear()
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(x[:-1], facet_fluxes[:, 0, mode_n], marker=">", label="flux {} left".format(mode_n), color="b", ls="")
                plt.plot(x[1:], facet_fluxes[:, 1,  mode_n], marker=">", label="flux {} right".format(mode_n), color="r", ls="")


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

            for mode_n in range(n_el_nod):
                plt.figure("Flux {}".format(mode_n))
                fig = plt.gcf()
                fig.clear()
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(xx, ww[:, 0], label="recon")
                plt.plot(X, cell_fluxes[:, mode_n], marker="D", label="cell flux {}".format(mode_n), color="r", ls="")


                # Plot mesh
                XN = x[-1]
                X1 = x[0]
                Xvol = XN - X1
                X = (x[1:] + x[:-1]) / 2
                plt.vlines(x[:, 0], ymin=0, ymax=.5, colors="grey")
                plt.vlines((X1, XN), ymin=0, ymax=.5, colors="k")
                plt.vlines(X, ymin=0, ymax=.3, colors="grey", linestyles="--")
                plt.legend()
                plt.close('all')


        # 2D plots
        if False:
            # TODO remove 2D plots
            import matplotlib.pyplot as plt
            import sfepy.postprocess.plot_cmesh as pc

            cmesh = self.region.domain.mesh.cmesh

            def plot_facet_normals(ax, cmesh, normals, scale, col='m'):
                dim = cmesh.dim
                ax = pc._get_axes(ax, dim)

                edim = dim - 1
                coors = cmesh.get_centroids(edim)
                coors = pc._to2d(coors)

                cmesh.setup_connectivity(dim, edim)
                c2f = cmesh.get_conn(dim, edim)
                for ic, o1 in enumerate(c2f.offsets[:-1]):
                    o2 = c2f.offsets[ic + 1]
                    for ifal, ifa in enumerate(c2f.indices[o1:o2]):
                        # print(ic, ifal, ifa)
                        cc = nm.array([coors[ifa], coors[ifa] + scale * normals[ic, ifal]])
                        # print(cc)
                        ax.plot(*cc.T, color=col)
                        # ax.plot([cc[1, 0]], [cc[1, 1]], color=col, marker=">")


                return ax


            ax = pc.plot_cmesh(
                None, cmesh,
                wireframe_opts={'color': 'k', 'linewidth': 2},
                entities_opts=[
                    {'color': 'k', 'label_global': 12, 'label_local': 8, 'size': 20},  # vertices
                    {'color': 'b', 'label_global': 6, 'label_local': 8, 'size': 10},  # faces
                    {'color': 'r', 'label_global': 12, 'size': 1},  # cells
                ])
            # for i in range(n_el_nod):
            #     ax = plot_facet_normals(ax, cmesh, facet_fluxes[:, :, i, None] * fc_n, 1, col='r')

            ax = plot_facet_normals(ax, cmesh, fc_n, .01, col='m')
            plt.show()

        out[:] = 0.0
        for i in range(n_el_nod):
            out[:, :, i, 0] = cell_fluxes[:, i, None]

        status = None
        return status


class DiffusionDGFluxTerm(Term):
    """
    Combines Lax-Friedrichs fluxes to get diffusion
    """


    alf = 0
    name = "dw_dg_difusion_laxfrie_flux"
    modes = ("weak",)
    arg_types = ('opt_material', 'material_diff_tensor', 'virtual', 'state')
    arg_shapes = [{'opt_material': '.: 1',
                   'material_diff_tensor': 'D, D',
                   'virtual': (1, 'state'),
                   'state': 1
                   },
                  {'opt_material': None}]
    integration = 'volume'
    symbolic = {'expression': 'div(D*grad(u))',
                'map': {'u': 'state', 'a': 'material', 'v': 'virtual'}
                }

    def get_fargs(self, alpha, diff_tensor, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        # TODO this should be kinda the same, use some abstraction - in NonlinearHyperDGFluxTerm too!
        if alpha is not None:
            self.alf = alpha # extract alpha value regardless of shape
        if diff_var is not None:
            # TODO will this term be residual only?
            output("Diff var is not None in residual only term {} ! Skipping.".format(self.name))
            return None, None, None, None, None, diff_tensor[:, 0, :, :], 0
        else:
            field = state.field
            region = field.region

            if not "DG" in field.family_name:
                raise ValueError("Used DG term with non DG field {} of family {}".format(field.name, field.family_name))

            dofs = unravel_sol(state)
            cell_normals = field.get_cell_normals_per_facet(region)
            facet_base_vals = field.get_facet_base(base_only=True)
            facet_base_dvals = field.get_facet_base(base_only=True, derivative=True)

            inner_facet_qp_vals, outer_facet_qp_vals, weights = field.get_both_facet_qp_vals(dofs, region)
            inner_facet_qp_dvals, outer_facet_qp_dvals, weights = field.get_both_facet_qp_vals(dofs, region, derivative=True)

            fargs = (dofs, inner_facet_qp_dvals, outer_facet_qp_dvals,
                     # TODO get correct shape of diffusion tensor
                     facet_base_dvals, weights, cell_normals, diff_tensor[:, 0, :, :])
            return fargs

        def function(self, out, dofs, in_fc_dv, out_fc_dv, fc_db, whs, fc_n):
            ...


class NonlinearHyperDGFluxTerm(Term):

    alf = 0
    name = "dw_dg_nonlinear_laxfrie_flux"
    modes = ("weak",)
    arg_types = ('opt_material', 'material_fun', 'material_fun_d', 'virtual', 'state')
    arg_shapes = [{'opt_material': '.: 1',
                   'material_fun': '.: 1',
                   'material_fun_d': '.: 1',
                   'virtual': (1, 'state'),
                   'state': 1
                   },
                  {'opt_material': None}]
    integration = 'volume'
    symbolic = {'expression' : 'div(f(u))*w',
                'map': {'u': 'state', 'v' : 'virtual', 'f': 'function'}
    }
    # def __init__(self, integral, region, f=None, df=None, **kwargs):
    #     AdvectDGFluxTerm.__init__(self, self.name + "(v, u)", "v, u[-1]", integral, region, **kwargs)
    #     # TODO how to pass f and df as Term arguments i.e. from conf examples?
    #     if f is not None:
    #         self.fun = f
    #     else:
    #         raise ValueError("Function f not provided to {}!".format(self.name))
    #     if df is not None:
    #         self.dfun = df
    #     else:
    #         raise ValueError("Derivative of function {} no provided to {}".format(self.fun, self.name))

    def get_fargs(self, alpha, fun, dfun, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        # FIXME - just  to get it working, remove code duplicity!

        if alpha is not None:
            # FIXME this is only hotfix to get scalar!
            self.alf = nm.max(alpha)  # extract alpha value regardless of shape

        self.fun = fun
        self.dfun = dfun

        if diff_var is not None:
            output("Diff var is not None in residual only term {} ! Skipping.".format(self.name))
            return None, None, None, None, None, 0
        else:
            field = state.field
            region = field.region

            if not "DG" in field.family_name:
                raise ValueError(
                    "Used DG term with non DG field {} of family {}!".format(field.name, field.family_name))

            dofs = unravel_sol(state)
            cell_normals = field.get_cell_normals_per_facet(region)
            facet_base_vals = field.get_facet_base(base_only=True)
            inner_facet_qp_vals, outer_facet_qp_vals, weights = field.get_both_facet_qp_vals(dofs, region)

            fargs = (dofs, inner_facet_qp_vals, outer_facet_qp_vals,
                     facet_base_vals, weights, cell_normals)
            return fargs

    # noinspection PyUnreachableCode
    def function(self, out, dofs, in_fc_v, out_fc_v, fc_b, whs, fc_n):
        """

        :param out:
        :param dofs: NOT necessary but was good for debugging, shape = (n_cell, n_el_nod)
        :param in_fc_v: inner values for facets per cell, shape = (n_cell, n_el_faces, 1)
        :param out_fc_v: outer values for facets per cell, shape = (n_cell, n_el_faces, 1)
        :param fc_b: values of basis in facets qps, shape = (1, n_el_facet, n_qp)
        :param whs: weights of the qps on facets, shape = (n_cell, n_el_facet, n_qp
        :param fc_n: facet normals, shape = (n_cell, n_el_facets)
        :param advelo: advection velocity, shape = (n_cell, 1)
        :return:
        """
        # FIXME - just  to get it working, remove code duplicity!
        if dofs is None:
            out[:] = 0.0
            return None

        f = self.fun
        df = self.dfun

        n_cell = dofs.shape[0]
        n_el_nod = dofs.shape[1]
        n_el_facets = fc_n.shape[-2]
        #  Calculate integrals over facets representing Lax-Friedrichs fluxes
        facet_fluxes = nm.zeros((n_cell, n_el_facets, n_el_nod))

        # get maximal wave speeds at facets
        df_in = df(in_fc_v)
        df_out = df(out_fc_v)
        fc_n__dot__df_in = nm.sum(fc_n[:,:,None, :] * df_in, axis=-1 )
        fc_n__dot__df_out = nm.sum(fc_n[:, :, None, :] * df_out, axis=-1)
        dfdn = nm.stack((fc_n__dot__df_in, fc_n__dot__df_out), axis=-1)
        C = nm.amax(nm.abs(dfdn), axis=(-2, -1))[..., None, None]

        for facet_n in range(n_el_facets):
            for mode_n in range(n_el_nod):
                fc_v_m = in_fc_v[:, facet_n, :] - out_fc_v[:, facet_n, :]

                central = (f(in_fc_v[:, facet_n, :]) + f(out_fc_v[:, facet_n, :])) / 2.
                upwind = (1 - self.alf) / 2. * C[:, facet_n, :, :] * fc_n[:, facet_n][..., None, :] * fc_v_m[:, :, None]

                facet_fluxes[:, facet_n, mode_n] = nm.sum(fc_n[:, facet_n] *
                                                          nm.sum((central + upwind) *
                                                                 (fc_b[None, :, 0, facet_n, 0, mode_n] *
                                                                  whs[:, facet_n, :])[..., None], axis=1), # sum over qps to get integral
                                                          axis=1) # sum over dimension to get scalar product

        cell_fluxes = nm.sum(facet_fluxes, axis=1)

        out[:] = 0.0
        for i in range(n_el_nod):
            out[:, :, i, 0] = cell_fluxes[:, i, None]

        status = None
        return status


from sfepy.linalg import dot_sequences

class NonlinScalarDotGradTerm(Term):
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
    name = 'dw_ns_dot_grad_s'
    arg_types = (( 'material_fun', 'material_fun_d','virtual', 'state'),
                 ( 'material_fun', 'material_fun_d', 'state', 'virtual'))
    arg_shapes = [{ 'material_fun': '.: 1',
                   'material_fun_d': '.: 1',
                    'virtual/grad_state' : (1, None),
                   'state/grad_state' : 1,
                   'virtual/grad_virtual' : (1, None),
                   'state/grad_virtual' : 1}]
    modes = ('grad_state', 'grad_virtual')

    # def __init__(self, integral, region, f=None, df=None, **kwargs):
    #     ScalarDotMGradScalarTerm.__init__(self, self.name + "(v, u)", "u[-1], v", integral, region, **kwargs)
    #     # TODO how to pass f and df as Term arguments i.e. from conf examples?
    #     if f is not None:
    #         self.fun = f
    #     else:
    #         raise ValueError("Function f not provided to {}!".format(self.name))
    #     if df is not None:
    #         self.dfun = df
    #     else:
    #         raise ValueError("Derivative of function {} no provided to {}".format(self.fun, self.name))


    @staticmethod
    def function(out, out_qp, geo, fmode):
        status = geo.integrate(out, out_qp)
        return status

    def get_fargs(self, fun, fund, var1, var2,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg1, _ = self.get_mapping(var1)
        vg2, _ = self.get_mapping(var2)

        if diff_var is None:
            if self.mode == 'grad_state':
                # TODO chceck correct shapes for integration
                geo = vg1
                bf_t = vg1.bf.transpose((0, 1, 3, 2))
                val_qp = dfun(self.get(var2, 'val')[..., 0])
                val_grad_qp = self.get(var2, 'grad')
                val = dot_sequences(val_qp, val_grad_qp, 'ATB')
                out_qp = dot_sequences(bf_t , val_grad_qp, 'ATB')

            else:
                # TODO chceck correct shapes for integration
                geo = vg2
                val_qp = fun(self.get(var1, 'val'))[..., 0, :].swapaxes(-2, -1)
                out_qp = dot_sequences(vg2.bfg, val_qp, 'ATB')

            fmode = 0

        else:
            # TODO what in matrix mode?
            if self.mode == 'grad_state':
                geo = vg1
                bf_t = vg1.bf.transpose((0, 1, 3, 2))
                out_qp = dot_sequences(bf_t , vg2.bfg, 'ATB')

            else:
                geo = vg2
                out_qp = dot_sequences(vg2.bfg, vg1.bf, 'ATB')

            fmode = 1

        return out_qp, geo, fmode






