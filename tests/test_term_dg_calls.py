# -*- coding: utf-8 -*-
"""
Test all terms in terms_dg. Performs numerical test on simple mesh.
"""
import functools
import inspect

import numpy as nm
import numpy.testing as nmts
import scipy.sparse as sp

from base.base import Struct
from base.testing import TestCommon
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


class Test(TestCommon):

    def capture_assertion_decorator(self, method):

        @functools.wraps(method)
        def captured_assertion_method(_):
            try:
                method()
            except AssertionError:
                return False
            return True

        return captured_assertion_method.__get__(self, self.__class__)

    @staticmethod
    def from_conf(conf, options):
        """
        Filters out terms test classes and gathers their test methods in
        resulting object.
        """
        term_test_classes = [(key, var) for key, var in dict(globals()).items()
                   if (key.startswith("Test") and key.endswith("Term"))]

        all_test = Test()
        for cname, term_test_cls in term_test_classes:
            term_test = term_test_cls()
            methods = inspect.getmembers(term_test, inspect.ismethod)
            all_test.update({f"{mname}_{cname[4:]}":
                             all_test.capture_assertion_decorator(meth)
                         for mname, meth in methods})
        return all_test

class DGTermTestEnvornment:
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
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = AdvectionDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u[-1], v",
                                   te.integral, te.regions["omega"],
                                   u=te.u, v=te.v, a=te.a)

        # te.u.data[0][::te.n_el_nod] = 1

        result = nm.zeros(te.out.shape)

        out, _ = term.function(te.out,
                               te.u,
                               None,  # diff_var
                               te.field,
                               te.regions["omega"],
                               te.a.data
                               )

        nmts.assert_almost_equal(out, result)

        return True

    def test_function_implicit_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = AdvectionDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u, v",
                                   te.integral, te.regions["omega"],
                                   u=te.u, v=te.v, a=te.a)

        # te.u.data[0][::ts.n_el_nod] = 1
        expected = nm.zeros(((te.n_cell * te.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            te.out,  # out, note that for implicit mode the out
                     # argument is ignored
            te.u,  # state
            "u",  # diff_var
            te.field,
            te.regions["omega"],
            te.a.data,  # advelo
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((te.n_cell * te.n_el_nod),) * 2).toarray()

        assert expected.shape == out.shape

        return True



class TestNonlinearHyperDGFluxTerm:

    def test_function_explicit_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = NonlinearHyperbolicDGFluxTerm("adv_stiff(f, df u, v)",
                                        "nonlin.f, nonlin.df, u[-1], v",
                                             te.integral, te.regions["omega"],
                                             u=te.u, v=te.v, nonlin=te.nonlin)

        # te.u.data[0][::ts.n_el_nod] = 1
        result = nm.zeros(te.out.shape)

        out, _ = term.function(te.out,
                               te.u,
                               te.field,
                               te.regions["omega"],
                               te.burg_fun,
                               te.burg_fun_d
                               )

        nmts.assert_almost_equal(out, result)

        return True



class TestDiffusionDGFluxTerm:

    def test_function_explicit_right_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)",
                                   "D.val, v,  u[-1]",
                                   te.integral, te.regions["omega"],
                                   u=te.u, v=te.v, D=te.D)
        term.mode = "avg_state"

        result = nm.zeros(te.out.shape)

        out, _ = term.function(te.out,  # out
                               te.u,  # state
                               None,  # diff_var, explicit
                               te.field,
                               te.regions["omega"],
                               te.D.data,  # advelo
                               )

        nmts.assert_almost_equal(out, result)

        return True


    def test_function_explicit_left_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, u, v)",
                                   "D.val,  u[-1], v",
                                   te.integral, te.regions["omega"],
                                   u=te.u, v=te.v, D=te.D)
        term.mode = "avg_virtual"

        result = nm.zeros(te.out.shape)

        out, _ = term.function(te.out,  # out
                               te.u,  # state
                               None,  # diff_var, explicit
                               te.field,
                               te.regions["omega"],
                               te.D.data,  # advelo
                               )

        nmts.assert_almost_equal(out, result)

        return True


    def test_function_implicit_right_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)",
                                   "D.val, v,  u",
                                   te.integral, te.regions["omega"],
                                   u=te.u, v=te.v, D=te.D)
        term.mode = "avg_state"

        expected = nm.zeros(((te.n_cell * te.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            te.out,  # out
            te.u,  # state
            "u",  # diff_var, explicit
            te.field,
            te.regions["omega"],
            te.D.data,  # advelo
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((te.n_cell * te.n_el_nod),) * 2).toarray()

        assert expected.shape == out.shape

        return True


    def test_function_implicit_left_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, u, v)",
                                   "D.val,  u[-1], v",
                                   te.integral, te.regions["omega"],
                                   u=te.u, v=te.v, D=te.D)
        term.mode = "avg_virtual"

        expected = nm.zeros(((te.n_cell * te.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            te.out,  # out
            te.u,  # state
            "u",  # diff_var, explicit
            te.field,
            te.regions["omega"],
            te.D.data,  # advelo
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((te.n_cell * te.n_el_nod),) * 2).toarray()

        assert expected.shape == out.shape

        return True



class TestDiffusionInteriorPenaltyTerm:

    def test_function_explicit_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = DiffusionInteriorPenaltyTerm("adv_stiff(Cw.val, u, v)",
                                            "Cw.val, u[-1], v",
                                            te.integral, te.regions["omega"],
                                            u=te.u, v=te.v, Cw=te.Cw)

        # te.u.data[0][::ts.n_el_nod] = 1

        result = nm.zeros(te.out.shape)

        out, _ = term.function(te.out,
                               te.u,
                               None,  # diff_var
                               te.field,
                               te.regions["omega"],
                               te.Cw.data,
                               te.D.data
                               )

        nmts.assert_almost_equal(out, result)

        return True


    def test_function_implicit_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = DiffusionInteriorPenaltyTerm("adv_stiff(D.val, a.val, u, v)",
                                            "Cw.val, u, v",
                                            te.integral, te.regions["omega"],
                                            u=te.u, v=te.v, a=te.Cw)

        # te.u.data[0][::ts.n_el_nod] = 1


        expected = nm.zeros(((te.n_cell * te.n_el_nod),) * 2)

        (out, iel1, iel2, _, _), _ = term.function(
            te.out,  # out, note that for implicit mode the out
                     # argument is ignored
            te.u,  # state
            "u",  # diff_var
            te.field,
            te.regions["omega"],
            te.Cw.data,
            te.D.data,
        )

        out = sp.csr_matrix((out, (iel1, iel2)),
                            shape=((te.n_cell * te.n_el_nod),) * 2).toarray()
        assert expected.shape == out.shape

        return True



class TestNonlinScalarDotGradTerm:

    def test_function_explicit_1D(self):
        te = DGTermTestEnvornment(dim=1, approx_order=3)

        term = NonlinearScalarDotGradTerm("adv_stiff(f, df u, v)",
                                       "nonlin.f, nonlin.df, u[-1], v",
                                          te.integral, te.regions["omega"],
                                          u=te.u, v=te.v, nonlin=te.nonlin)
        term.setup()

        # te.u.data[0][::ts.n_el_nod] = 1
        expected = nm.zeros(te.out.shape)

        out = nm.zeros(te.out.shape)

        fargs = term.get_fargs(
            te.burg_fun,
            te.burg_fun_d,
            te.u,
            te.v
        )
        fargs = (out,) + fargs
        out = term.function(*fargs)

        nmts.assert_almost_equal(out, expected)

        return True


if __name__ == '__main__':
    t = Test()
    t.test_dg_term_calls()