# -*- coding: utf-8 -*-
"""
Test all terms in dg_terms. Performs numerical test on simple meshes with
hand calculated values.
"""
import numpy as nm
import numpy.testing as nmts
import scipy.sparse as sp

from base.base import Struct
from discrete import Variables, Materials
from discrete.common.dof_info import EquationMap
from sfepy.discrete import (DGFieldVariable, Material, Integral,
                            Function, Equation, Equations, Problem)

from sfepy.discrete.dg.fields import DGField
from sfepy.terms.terms_dg import DGTerm, \
    AdvectionDGFluxTerm, NonlinearHyperbolicDGFluxTerm, NonlinearScalarDotGradTerm, \
    DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm

from sfepy.discrete.dg.tests.test_dg_field import prepare_field_1D, \
    prepare_field_2D


class DGTermTestScope:
    """
    Class for  easy creation of all the data needed for testing terms.
    """

    def burg_fun(self, u):
        vu = self.burg_velo * u[..., None] ** 2
        return vu

    def burg_fun_d(self, u):
        v1 = 2 * self.burg_velo * u[..., None]
        return v1

    def __init__(self, dim, approx_order, **kwargs):
        """
        Creates Struct object with all the data necessary to test terms

        :param dim: dimension
        :param approx_order: approximation order
        :param kwargs: velo, diffusion or penalty for prepare_materials
        :return: term test scope
        """

        if dim == 1:
            (field, regions), mesh = prepare_field_1D(approx_order)
        elif dim == 2:
            (field, regions), mesh = prepare_field_2D(approx_order)

        self.field = field
        self.regions = regions
        self.mesh = mesh

        self.n_cell = field.n_cell
        self.n_nod = field.n_nod
        self.n_el_nod = field.n_el_nod

        self.u, self.v = self.prepare_variables(field)
        self.u.data = [(nm.zeros(self.n_nod))]
        self.variables = Variables([ self.u, self.v])

        self.integral = Integral('i', order=approx_order * 2)
        self.a, self.D, self.Cw = self.prepare_materials(field, **kwargs)

        if dim == 1:
            velo = nm.array(1.0)
        elif dim == 2:
            velo = nm.array([1.0, 0])

        self.burg_velo = velo.T / nm.linalg.norm(velo)

        self.nonlin = Material('nonlin',
                               values={'.fun': self.burg_fun,
                                       '.dfun': self.burg_fun_d})

        self.out = nm.zeros((self.n_cell, 1, self.n_el_nod, 1))

    def prepare_variables(self, field):
        """
        Prepares state and test variables, adds empty
        eq_map to state variable

        :param field:
        :return: state, test
        """
        n_nod = field.n_nod

        u = DGFieldVariable('u', 'unknown', field, history=1)
        v = DGFieldVariable('v', 'test', field, primary_var_name='u')
        var_di = Struct(
            details=Struct(dpn=1, n_nod=n_nod,
                           name="field_var_dof_details"),
            indx=slice(0, n_nod, None), n_dof=n_nod, name='u_dof_info',
            var_name="u")

        u.eq_map = EquationMap("eq_map", ["u.0"], var_di)
        u.eq_map._init_empty(field)

        return u, v

    def prepare_materials(self, field, velo=1.0, diffusion=0.1, penalty=100):
        """
        Crates material objects with data attribute, containing properly shaped
        data to pass to terms

        :param field: DGField
        :param velo: optional values for velocity a
        :param diffusion: optional value for diffusion tensor D
        :param penalty: optional value for diffusion penalty Cw
        :return: a, D, Cw
        """

        a = Material('a', val=[velo])
        a.data = nm.ones((field.n_cell, 1)) * velo

        D = Material('D', val=[diffusion])
        D.data = nm.ones((field.n_cell, 1, 1)) * diffusion

        Cw = Material("Cw", values={".val": penalty})
        Cw.data = penalty

        return a, D, Cw


class TestAdvectDGFluxTerm:

    def test_function_explicit_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = AdvectionDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u[-1], v",
                                   ts.integral, ts.regions["omega"],
                                   u=ts.u, v=ts.v, a=ts.a)

        # ts.u.data[0][::ts.n_el_nod] = 1

        result = nm.zeros(ts.out.shape)

        out, _ = term.function(ts.out,
                               ts.u,
                               None,  # diff_var
                               ts.field,
                               ts.regions["omega"],
                               ts.a.data
                               )

        nmts.assert_almost_equal(out, result)

    def test_function_implicit_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = AdvectionDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u, v",
                                   ts.integral, ts.regions["omega"],
                                   u=ts.u, v=ts.v, a=ts.a)

        # ts.u.data[0][::ts.n_el_nod] = 1
        result = nm.zeros(((ts.n_cell * ts.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            ts.out,  # out, note that for implicit mode the out
                     # argument is ignored
            ts.u,  # state
            "u",  # diff_var
            ts.field,
            ts.regions["omega"],
            ts.a.data,  # advelo
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((ts.n_cell * ts.n_el_nod),) * 2).toarray()

        nmts.assert_almost_equal(out, result)


class TestNonlinearHyperDGFluxTerm:

    def test_function_explicit_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = NonlinearHyperbolicDGFluxTerm("adv_stiff(f, df u, v)",
                                        "nonlin.f, nonlin.df, u[-1], v",
                                             ts.integral, ts.regions["omega"],
                                             u=ts.u, v=ts.v, nonlin=ts.nonlin)

        # ts.u.data[0][::ts.n_el_nod] = 1
        result = nm.zeros(ts.out.shape)

        out, _ = term.function(ts.out,
                               ts.u,
                               ts.field,
                               ts.regions["omega"],
                               ts.burg_fun,
                               ts.burg_fun_d
                               )

        nmts.assert_almost_equal(out, result)


class TestDiffusionDGFluxTerm:

    def test_function_explicit_right_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)",
                                   "D.val, v,  u[-1]",
                                   ts.integral, ts.regions["omega"],
                                   u=ts.u, v=ts.v, D=ts.D)
        term.mode = "avg_state"

        result = nm.zeros(ts.out.shape)

        out, _ = term.function(ts.out,  # out
                               ts.u,  # state
                               None,  # diff_var, explicit
                               ts.field,
                               ts.regions["omega"],
                               ts.D.data,  # advelo
                               )

        nmts.assert_almost_equal(out, result)

    def test_function_explicit_left_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, u, v)",
                                   "D.val,  u[-1], v",
                                   ts.integral, ts.regions["omega"],
                                   u=ts.u, v=ts.v, D=ts.D)
        term.mode = "avg_virtual"

        result = nm.zeros(ts.out.shape)

        out, _ = term.function(ts.out,  # out
                               ts.u,  # state
                               None,  # diff_var, explicit
                               ts.field,
                               ts.regions["omega"],
                               ts.D.data,  # advelo
                               )

        nmts.assert_almost_equal(out, result)

    def test_function_implicit_right_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)",
                                   "D.val, v,  u",
                                   ts.integral, ts.regions["omega"],
                                   u=ts.u, v=ts.v, D=ts.D)
        term.mode = "avg_state"

        result = nm.zeros(((ts.n_cell * ts.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            ts.out,  # out
            ts.u,  # state
            "u",  # diff_var, explicit
            ts.field,
            ts.regions["omega"],
            ts.D.data,  # advelo
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((ts.n_cell * ts.n_el_nod),) * 2).toarray()

        nmts.assert_almost_equal(out, result)

    def test_function_implicit_left_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, u, v)",
                                   "D.val,  u[-1], v",
                                   ts.integral, ts.regions["omega"],
                                   u=ts.u, v=ts.v, D=ts.D)
        term.mode = "avg_virtual"

        result = nm.zeros(((ts.n_cell * ts.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            ts.out,  # out
            ts.u,  # state
            "u",  # diff_var, explicit
            ts.field,
            ts.regions["omega"],
            ts.D.data,  # advelo
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((ts.n_cell * ts.n_el_nod),) * 2).toarray()

        nmts.assert_almost_equal(out, result)


class TestDiffusionInteriorPenaltyTerm:

    def test_function_explicit_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = DiffusionInteriorPenaltyTerm("adv_stiff(Cw.val, u, v)",
                                            "Cw.val, u[-1], v",
                                            ts.integral, ts.regions["omega"],
                                            u=ts.u, v=ts.v, Cw=ts.Cw)

        # ts.u.data[0][::ts.n_el_nod] = 1

        result = nm.zeros(ts.out.shape)

        out, _ = term.function(ts.out,
                               ts.u,
                               None,  # diff_var
                               ts.field,
                               ts.regions["omega"],
                               ts.Cw.data
                               )

        nmts.assert_almost_equal(out, result)

    def test_function_implicit_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = DiffusionInteriorPenaltyTerm("adv_stiff(a.val, u, v)",
                                            "Cw.val, u, v",
                                            ts.integral, ts.regions["omega"],
                                            u=ts.u, v=ts.v, a=ts.Cw)

        # ts.u.data[0][::ts.n_el_nod] = 1

        result = nm.zeros(((ts.n_cell * ts.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            ts.out,  # out, note that for implicit mode the out
                     # argument is ignored
            ts.u,  # state
            "u",  # diff_var
            ts.field,
            ts.regions["omega"],
            ts.a.data,  # advelo
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((ts.n_cell * ts.n_el_nod),) * 2).toarray()

        nmts.assert_almost_equal(out, result)


class TestNonlinScalarDotGradTerm:

    def test_function_explicit_1D(self):
        ts = DGTermTestScope(dim=1, approx_order=3)

        term = NonlinearScalarDotGradTerm("adv_stiff(f, df u, v)",
                                       "nonlin.f, nonlin.df, u[-1], v",
                                          ts.integral, ts.regions["omega"],
                                          u=ts.u, v=ts.v, nonlin=ts.nonlin)

        # ts.u.data[0][::ts.n_el_nod] = 1
        result = nm.zeros(ts.out.shape)

        fargs = term.get_fargs(
            ts.burg_fun,
            ts.burg_fun_d,
            ts.u,
            ts.v
        )
        out, _ = term.function(*fargs)

        nmts.assert_almost_equal(out, result)
