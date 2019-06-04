import numpy as nm

# sfepy imports
from sfepy.terms.terms import Term, terms
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)
from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm

from sfepy.discrete.dg.dg_field import get_unraveler, get_raveler, DGField


class DGTerm(Term):
    def call_function(self, out, fargs):
        try:
            out, status = self.function(out, *fargs)

        except (RuntimeError, ValueError):
            terms.errclear()
            raise

        if status:
            terms.errclear()
            raise ValueError('term evaluation failed! (%s)' % self.name)

        return out, status

    def eval_real(self, shape, fargs, mode='eval', term_mode=None,
                  diff_var=None, **kwargs):
        out = nm.empty(shape, dtype=nm.float64)

        if mode == 'weak':
            out, status = self.call_function(out, fargs)

        else:
            status = self.call_function(out, fargs)

        return out, status

    def _get_nbrhd_dof_indexes(self, active_cells, active_nrbhs, field):
        """
        Get indexes of DOFs of neighbouring cells, maybe move to DGField?
        :param active_cells:
        :param active_nrbhs:
        :param field:
        :return:
        """
        inner_iels = field.bubble_dofs
        inner_iels = nm.stack((nm.repeat(inner_iels, field.n_el_nod), nm.tile(inner_iels, field.n_el_nod).flatten()),
                              axis=-1)
        outer_iels = nm.stack((nm.repeat(field.bubble_dofs[active_cells], field.n_el_nod),
                               nm.tile(field.bubble_dofs[active_nrbhs], field.n_el_nod).flatten()), axis=-1)
        iels = nm.vstack((inner_iels, outer_iels))
        return iels

class AdvectDGFluxTerm(DGTerm):
    r"""
    Lax-Friedrichs flux term for advection of scalar quantity :math:`p` with the advection velocity
    :math:`\ul{a}` given as a material parameter (a known function of space and time).

    :Definition:

    .. math::

        \int_{\partial{T_K}} \vec{n} \cdot f^{*} (p_{in}, p_{out})q

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

    alpha = 0
    name = "dw_dg_advect_laxfrie_flux"
    modes = ("weak",)
    arg_types = ('opt_material', 'material_advelo', 'virtual', 'state')
    arg_shapes = [{'opt_material'   : '.: 1',
                   'material_advelo': 'D, 1',
                   'virtual'        : (1, 'state'),
                   'state'          : 1
                   },
                  {'opt_material': None}]
    integration = 'volume'
    symbolic = {'expression': 'div(a*u)*w',
                'map'       : {'u': 'state', 'a': 'material', 'v': 'virtual'}
                }

    def get_fargs(self, alpha, advelo, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        if alpha is not None:
            self.alpha = alpha  # extract alpha value regardless of shape

        field = state.field
        region = field.region

        if "DG" not in field.family_name:
            raise ValueError("Used DG term with non DG field {} of family {}".format(field.name, field.family_name))

        fargs = (state, test, diff_var, field, region, advelo[:, 0, :, 0])
        return fargs

    # noinspection PyUnreachableCode
    def function(self, out, state, test, diff_var, field : DGField, region, advelo):

        if diff_var is not None:

            fc_n = field.get_cell_normals_per_facet(region)
            C = nm.abs(nm.einsum("ifk,ik->if", fc_n, advelo))

            nbrhd_idx = field.get_facet_neighbor_idx(region, state.eq_map)
            active_cells, active_facets = nm.where(nbrhd_idx[:, :, 0] >= 0)
            active_nrbhs = nbrhd_idx[active_cells, active_facets, 0]
            active_nrbhs_facets = nbrhd_idx[active_cells, active_facets, 1]

            in_fc_b, out_fc_b, whs = field.get_both_facet_base_vals(state, region)

            # compute values
            inner_diff = nm.einsum("nfk, nfk->nf",
                                   fc_n,
                                   advelo[:, None, :] + nm.einsum("nfk, nf->nfk", (1 - self.alpha) * fc_n, C)) / 2.
            # FIXME broadcast advelo to facets - maybe somehow get values of advelo at them?

            outer_diff = nm.einsum("nfk, nfk->nf",
                                   fc_n,
                                   advelo[:, None, :] - nm.einsum("nfk, nf->nfk", (1 - self.alpha) * fc_n, C)) / 2.

            inner_vals = nm.einsum("nf, ndfq, nbfq, nfq -> ndb", inner_diff, in_fc_b, in_fc_b, whs)
            outer_vals = nm.einsum("i, idq, ibq, iq -> idb",
                                             outer_diff[active_cells, active_facets],
                                             in_fc_b[active_cells, :, active_facets],
                                             out_fc_b[active_cells, :, active_facets],
                                             whs[active_cells, active_facets])

            vals = nm.vstack((inner_vals, outer_vals))
            vals = vals.flatten()

            # compute postions within matrix
            iels = self._get_nbrhd_dof_indexes(active_cells, active_nrbhs, field)

            out = (vals, iels[:, 0], iels[:, 1], state, state)
        else:
            fc_n = field.get_cell_normals_per_facet(region)
            facet_base_vals = field.get_facet_base(base_only=True)
            in_fc_v, out_fc_v, weights = field.get_both_facet_state_vals(state, region)
            # get sane facet base shape
            fc_b = facet_base_vals[:, 0, :, 0, :].T  # (n_el_nod, n_el_facet, n_qp)

            # get maximal wave speeds at facets
            C = nm.abs(nm.einsum("ifk,ik->if", fc_n, advelo))

            fc_v_avg = (in_fc_v + out_fc_v)/2.
            fc_v_jmp = in_fc_v - out_fc_v

            central = nm.einsum("ik,ifq->ifkq", advelo, fc_v_avg)
            upwind = (1 - self.alpha) / 2. * nm.einsum("if,ifk,ifq->ifkq", C, fc_n, fc_v_jmp)

            cell_fluxes = nm.einsum("ifk,ifkq,dfq,ifq->id", fc_n, central + upwind, fc_b, weights)

            out[:] = 0.0
            n_el_nod = field.n_el_nod
            for i in range(n_el_nod):
                out[:, :, i, 0] = cell_fluxes[:, i, None]

        status = None
        return out, status


class DiffusionDGFluxTerm(DGTerm):
    r"""
    Basic diffusion term for scalar quantity.
    Offers two modes: avg_state and avg_virtual

    """
    name = "dw_dg_diffusion_flux"
    arg_types = (('material_diff_tensor', 'state', 'virtual'),
                 ('material_diff_tensor', 'virtual', 'state')
                 )
    arg_shapes = [{'material_diff_tensor': '1, 1',
                   'virtual/avg_state': (1, None),
                   'state/avg_state' : 1,
                   'virtual/avg_virtual': (1, None),
                   'state/avg_virtual' : 1,
                   }]
    integration = 'volume'
    modes = ('avg_state', 'avg_virtual')

    def get_fargs(self, diff_tensor, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        field = state.field
        region = field.region

        if "DG" not in field.family_name:
            raise ValueError("Used DG term with non DG field {} of family {}".format(field.name, field.family_name))
        if self.mode == "avg_state":
            # put state where it is expected by function
            state = test
        fargs = (state, diff_var, field, region, diff_tensor[:, 0, :, :])
        return fargs

    def function(self, out, state, diff_var, field, region, D):
        if diff_var is not None:
            fc_n = field.get_cell_normals_per_facet(region)

            nbrhd_idx = field.get_facet_neighbor_idx(region, state.eq_map)
            active_cells, active_facets = nm.where(nbrhd_idx[:, :, 0] >= 0)
            active_nrbhs = nbrhd_idx[active_cells, active_facets, 0]

            inner_facet_base, outer_facet_base, whs = field.get_both_facet_base_vals(state, region,
                                                                                     derivative=False
                                                                                     )

            inner_facet_base_d, outer_facet_base_d, _ = field.get_both_facet_base_vals(state, region,
                                                                                       derivative=True)
            if self.mode == 'avg_state':
                inner_vals = nm.einsum("nkl, nfk, ndfq, nbfkq, nfq->ndb", D, fc_n,
                                            inner_facet_base,  # test
                                            inner_facet_base_d/2,  # state
                                            whs)

                outer_vals = nm.einsum("ikl, ik, idq, ibkq, iq->idb",
                                  D[active_cells],
                                  fc_n[active_cells, active_facets],
                                  inner_facet_base[active_cells, :, active_facets],  # test
                                  outer_facet_base_d[active_cells, :, active_facets]/2,  # state
                                  whs[active_cells, active_facets])
            elif self.mode == 'avg_virtual':
                inner_vals = nm.einsum("nkl, nfk, ndfkq, nbfq, nfq->ndb", D, fc_n,
                                             inner_facet_base_d/2,  # test
                                             inner_facet_base,  # state
                                             whs)

                outer_vals = nm.einsum("ikl, ik, idkq, ibq, iq->idb",
                                             D[active_cells],
                                             fc_n[active_cells, active_facets],
                                             inner_facet_base_d[active_cells, :, active_facets]/2,  # test
                                             - outer_facet_base[active_cells, :, active_facets],  # state
                                             whs[active_cells, active_facets])

            iels = self._get_nbrhd_dof_indexes(active_cells, active_nrbhs, field)

            vals = nm.vstack((inner_vals, outer_vals))
            vals = vals.flatten()

            out = (vals, iels[:, 0], iels[:, 1], state, state)
        else:
            fc_n = field.get_cell_normals_per_facet(region)
            inner_facet_base, outer_facet_base, _ = field.get_both_facet_base_vals(state, region,
                                                                                   derivative=False
                                                                                   )
            inner_facet_state_d, outer_facet_state_d, _ = field.get_both_facet_state_vals(state, region,
                                                                                          derivative=True
                                                                                          )
            inner_facet_base_d, outer_facet_base_d, _ = field.get_both_facet_base_vals(state, region,
                                                                                       derivative=True)
            inner_facet_state, outer_facet_state, weights = field.get_both_facet_state_vals(state, region,
                                                                                            derivative=False
                                                                                            )
            if self.mode == 'avg_state':
                avgDdState = (nm.einsum("ikl,ifkq->ifkq", D, inner_facet_state_d) +
                              nm.einsum("ikl,ifkq->ifkq", D, outer_facet_state_d)) / 2.
                jmpBase = inner_facet_base  # outer_facet_base is in DG zero - hence the jump is inner value
                vals = nm.einsum("ifkq , ifk, idfq, ifq -> id", avgDdState, fc_n, jmpBase, weights)

            elif self.mode == 'avg_virtual':
                avgDdbase = (nm.einsum("ikl,idfkq->idfkq", D, inner_facet_base_d)) / 2.
                            # in DG test function is non zero only inside element - hence we average with zero
                jmpState = inner_facet_state - outer_facet_state
                vals = nm.einsum("idfkq, ifk, ifq , ifq -> id", avgDdbase, fc_n, jmpState, weights)

            cell_fluxes = vals

            out[:] = 0.0
            n_el_nod = field.n_el_nod
            for i in range(n_el_nod):
                out[:, :, i, 0] = cell_fluxes[:, i, None]

        status = None
        return out, status


class DiffusionInteriorPenaltyTerm(DGTerm):
    name = "dw_dg_interior_penal"
    modes = ("weak",)
    arg_types = ('material_Cw', 'virtual', 'state')
    arg_shapes = [{'material_Cw': '.: 1',
                   'virtual'    : (1, 'state'),
                   'state'      : 1
                   }]

    def get_fargs(self, Cw, test, state, mode=None, term_mode=None, diff_var=None, **kwargs):

        field = state.field
        region = field.region

        if "DG" not in field.family_name:
            raise ValueError("Used DG term with non DG field {} of family {}".format(field.name, field.family_name))

        fargs = (state, test, diff_var, field, region, Cw)
        return fargs

    def function(self, out, state, test, diff_var, field, region, Cw):

        if diff_var is not None:

            nbrhd_idx = field.get_facet_neighbor_idx(region, state.eq_map)
            active_cells, active_facets = nm.where(nbrhd_idx[:, :, 0] >= 0)
            active_nrbhs = nbrhd_idx[active_cells, active_facets, 0]

            inner_facet_base, outer_facet_base, whs = field.get_both_facet_base_vals(state, region,
                                                                                   derivative=False)

            facet_vols = nm.sum(whs, axis=-1)
            inv_facet_vols = 1./nm.sum(whs, axis=-1)

            sigma = Cw / facet_vols

            inner = nm.einsum("nf, ndfq, nbfq, nfq -> ndb",
                               sigma,
                               inner_facet_base,  # test
                               inner_facet_base,  # state
                               whs)

            outer = nm.einsum("i, idq, ibq, iq -> idb",
                               sigma[active_cells, active_facets],
                               inner_facet_base[active_cells, :, active_facets],  # test
                               - outer_facet_base[active_cells, :, active_facets],  # state
                               whs[active_cells, active_facets])

            vals = nm.vstack((inner, outer))
            vals = vals.flatten()

            iels = self._get_nbrhd_dof_indexes(active_cells, active_nrbhs, field)
            out = (vals, iels[:, 0], iels[:, 1], state, state)

            from scipy.sparse import coo_matrix
            extra = coo_matrix((vals, (iels[:, 0], iels[:, 1])),
                               shape=2 * (field.n_el_nod * field.n_cell,))
            fextra = extra.toarray()
            Mu = nm.dot(fextra, state.data[0]).reshape((field.n_cell, field.n_el_nod))
        else:
            inner_facet_state, outer_facet_state, whs = field.get_both_facet_state_vals(state, region,
                                                                                            derivative=False
                                                                                            )
            inner_facet_base, outer_facet_base, _ = field.get_both_facet_base_vals(state, region,
                                                                                   derivative=False)
            facet_vols = nm.sum(whs, axis=-1)

            jmp_state = inner_facet_state - outer_facet_state
            jmp_base = inner_facet_base  # - outer_facet_base
            sigma = Cw / facet_vols

            n_el_nod = nm.shape(inner_facet_base)[1]
            cell_penalty = nm.einsum("nf,nfq,ndfq,nfq->nd", sigma, jmp_state, jmp_base, whs)

            out[:] = 0.0
            for i in range(n_el_nod):
                out[:, :, i, 0] = cell_penalty[:, i, None]

        status = None
        return out, status


class NonlinearHyperDGFluxTerm(DGTerm):
    alf = 0
    name = "dw_dg_nonlinear_laxfrie_flux"
    modes = ("weak",)
    arg_types = ('opt_material', 'material_fun', 'material_fun_d', 'virtual', 'state')
    arg_shapes = [{'opt_material'  : '.: 1',
                   'material_fun'  : '.: 1',
                   'material_fun_d': '.: 1',
                   'virtual'       : (1, 'state'),
                   'state'         : 1
                   },
                  {'opt_material': None}]
    integration = 'volume'
    symbolic = {'expression': 'div(f(u))*w',
                'map'       : {'u': 'state', 'v': 'virtual', 'f': 'function'}
                }

    def get_fargs(self, alpha, fun, dfun, test, state,
                  mode=None, term_mode=None, diff_var=None, **kwargs):

        if alpha is not None:
            self.alf = nm.max(alpha)  # extract alpha value regardless of shape

        self.fun = fun
        self.dfun = dfun

        if diff_var is not None:
            output("Diff var is not None in nonlinear residual only term {} ! Skipping.".format(self.name))
            return None, None, None, 0, 0
        else:
            field = state.field
            region = field.region

            if "DG" not in field.family_name:
                raise ValueError(
                        "Used DG term with non DG field {} of family {}!".format(field.name, field.family_name))

            fargs = (state, field, region, fun, dfun)
            return fargs

    # noinspection PyUnreachableCode
    def function(self, out, state, field, region, f, df):
        if state is None:
            out[:] = 0.0
            return None

        fc_n = field.get_cell_normals_per_facet(region)
        facet_base_vals = field.get_facet_base(base_only=True)
        in_fc_v, out_fc_v, weights = field.get_both_facet_state_vals(state, region)

        fc_b = facet_base_vals[:, 0, :, 0, :].T  # (n_el_nod, n_el_facet, n_qp)

        n_cell = field.n_cell
        n_el_nod = field.n_el_nod
        n_el_facets = field.n_el_facets

        # get maximal wave speeds at facets
        df_in = df(in_fc_v)
        df_out = df(out_fc_v)
        fc_n__dot__df_in = nm.einsum("ifk,ifqk->ifq", fc_n, df_in)
        fc_n__dot__df_out = nm.einsum("ifk,ifqk->ifq", fc_n, df_out)
        dfdn = nm.stack((fc_n__dot__df_in, fc_n__dot__df_out), axis=-1)
        C = nm.amax(nm.abs(dfdn), axis=(-2, -1))

        fc_f_avg = (f(in_fc_v) + f(out_fc_v)) / 2.
        fc_v_jmp = in_fc_v - out_fc_v

        central = fc_f_avg
        upwind = (1 - self.alf) / 2. * nm.einsum("if,ifk,ifq->ifqk", C, fc_n, fc_v_jmp)

        cell_fluxes = nm.einsum("ifk,ifqk,dfq,ifq->id", fc_n, central + upwind, fc_b, weights)

        out[:] = 0.0
        for i in range(n_el_nod):
            out[:, :, i, 0] = cell_fluxes[:, i, None]

        status = None
        return out, status


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
    arg_types = (('material_fun', 'material_fun_d', 'virtual', 'state'),
                 ('material_fun', 'material_fun_d', 'state', 'virtual'))
    arg_shapes = [{'material_fun'        : '.: 1',
                   'material_fun_d'      : '.: 1',
                   'virtual/grad_state'  : (1, None),
                   'state/grad_state'    : 1,
                   'virtual/grad_virtual': (1, None),
                   'state/grad_virtual'  : 1}]
    modes = ('grad_state', 'grad_virtual')

    @staticmethod
    def function(out, out_qp, geo, fmode):
        status = geo.integrate(out, out_qp)
        return status

    def get_fargs(self, fun, dfun, var1, var2,
                  mode=None, term_mode=None, diff_var=None, **kwargs):
        vg1, _ = self.get_mapping(var1)
        vg2, _ = self.get_mapping(var2)

        if diff_var is None:
            if self.mode == 'grad_state':
                # TODO check correct shapes for integration
                geo = vg1
                bf_t = vg1.bf.transpose((0, 1, 3, 2))
                val_qp = dfun(self.get(var2, 'val')[..., 0])
                val_grad_qp = self.get(var2, 'grad')
                val = dot_sequences(val_qp, val_grad_qp, 'ATB')
                out_qp = dot_sequences(bf_t, val_grad_qp, 'ATB')

            else:
                geo = vg2
                val_qp = fun(self.get(var1, 'val'))[..., 0, :].swapaxes(-2, -1)
                out_qp = dot_sequences(vg2.bfg, val_qp, 'ATB')

            fmode = 0

        else:
            # TODO what in matrix mode?
            if self.mode == 'grad_state':
                geo = vg1
                bf_t = vg1.bf.transpose((0, 1, 3, 2))
                out_qp = dot_sequences(bf_t, vg2.bfg, 'ATB')

            else:
                geo = vg2
                out_qp = dot_sequences(vg2.bfg, vg1.bf, 'ATB')

            fmode = 1

        return out_qp, geo, fmode
