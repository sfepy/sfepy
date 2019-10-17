import numpy as nm
import numpy.testing as nmts

from base.base import Struct
from discrete.common.dof_info import EquationMap
from sfepy.discrete import (DGFieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.variables import Variables
from sfepy.discrete.state import State
from sfepy.mesh.mesh_generators import gen_block_mesh

from sfepy.discrete.dg.dg_basis import get_n_el_nod


from sfepy.discrete.dg.dg_field import DGField
from sfepy.discrete.dg.dg_terms import DGTerm, \
    AdvectDGFluxTerm, NonlinearHyperDGFluxTerm, NonlinScalarDotGradTerm, \
    DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm

from sfepy.discrete.dg.tests.test_dg_field import prepare_field_1D, \
    prepare_field_2D


def prepare_variables(field):
    """
    Prepares state and test variable, adds empty
    eq_map to state variable

    :param field:
    :return: state, test
    """
    n_nod = field.n_nod

    u = DGFieldVariable('u', 'unknown', field, history=1)
    v = DGFieldVariable('v', 'test', field, primary_var_name='u')
    var_di = Struct(
        details=Struct(dpn=1, n_nod=n_nod, name="field_var_dof_details"),
        indx=slice(0, n_nod, None), n_dof=n_nod, name='u_dof_info',
        var_name="u")

    u.eq_map = EquationMap("eq_map", ["u.0"], var_di)
    u.eq_map._init_empty(field)

    return u, v


class TestAdvectDGFluxTerm:

    def test_function_explicit_1D(self):
        approx_order = 3
        (field, regions), mesh = prepare_field_1D(approx_order)
        u, v = prepare_variables(field)
        integral = Integral('i', order=approx_order * 2)

        n_el_nod = field.n_el_nod
        n_cell = field.n_cell
        n_nod = n_cell * n_el_nod

        a = Material('a', val=[1.0])
        a.data = nm.ones((n_cell, 1)) * a.function()["val"][0]

        term = AdvectDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u[-1], v",
                                integral, regions["omega"],
                                  u=u, v=v, a=a)

        u.data = [(nm.ones(n_nod))]
        out = nm.zeros((n_cell, 1, n_el_nod, 1))

        out, _ = term.function(out,
                      u,
                      None,  # diff_var
                      field,
                      regions["omega"],
                      a.data
                      )
        pass



    def test_function_implicit_1D(self):
        approx_order = 3
        (field, regions), mesh = prepare_field_1D(approx_order)
        u, v = prepare_variables(field)
        integral = Integral('i', order=approx_order * 2)

        n_el_nod = field.n_el_nod
        n_cell = field.n_cell
        n_nod = n_cell * n_el_nod

        a = Material('a', val=[1.0])
        a.data = nm.ones((n_cell, 1)) * a.function()["val"][0]

        term = AdvectDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u, v",
                                integral, regions["omega"],
                                u=u, v=v, a=a)

        u.data = [(nm.ones(n_nod))]
        out = nm.zeros((n_cell, 1, n_el_nod, 1))


        out, _ = term.function(out,  # out
                      u,  # state
                      "u",  # diff_var
                      field,
                      regions["omega"],
                      a.data,  # advelo
                      )

        pass


class TestDiffusionDGFluxTerm:

    def test_function_explicit_right_1D(self):
        approx_order = 3
        (field, regions), mesh = prepare_field_1D(approx_order)
        u, v = prepare_variables(field)

        integral = Integral('i', order=approx_order * 2)

        n_el_nod = field.n_el_nod
        n_cell = field.n_cell
        n_nod = n_cell * n_el_nod

        u.data = [(nm.ones(n_nod))]
        out = nm.zeros((n_cell, 1, n_el_nod, 1))


        diffusion_tensor = .02
        D = Material('D', val=[diffusion_tensor])
        D.data = nm.ones((n_cell, 1, 1)) * diffusion_tensor

        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)",
                                   "D.val, v,  u[-1]",
                                   integral, regions["omega"], u=u, v=v, D=D)
        term.mode = "avg_state"

        out, _ = term.function(out,  # out
                      u,  # state
                      None,  # diff_var, explicit
                      field,
                      regions["omega"],
                      D.data,  # advelo
                      )
        pass