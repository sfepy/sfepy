import numpy as nm
import numpy.testing as nmts

from base.base import Struct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.mesh.mesh_generators import gen_block_mesh

from sfepy.discrete.dg.dg_basis import get_n_el_nod


from sfepy.discrete.dg.dg_field import DGField
from sfepy.discrete.dg.dg_terms import DGTerm, \
    AdvectDGFluxTerm, NonlinearHyperDGFluxTerm, NonlinScalarDotGradTerm, \
    DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm

from sfepy.discrete.dg.tests.test_dg_field import prepare_field_1D, prepare_field_2D


class TestAdvectDGFluxTerm:

    def test_function_explicit_1D(self):
        approx_order = 3
        (field, regions), mesh = prepare_field_1D(approx_order)
        u = FieldVariable('u', 'unknown', field, history=1)
        v = FieldVariable('v', 'test', field, primary_var_name='u')
        integral = Integral('i', order=approx_order * 2)
        a = Material('a', val=[1.0])

        term = AdvectDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u[-1], v",
                                integral, regions["omega"],
                                  u=u, v=v, a=a)
        # TODO create state
        term.function([],  # out
                      [],  # state
                      None, # diff_var
                      field,
                      regions["omega"],
                      [], # advelo
                      )

    def test_function_implicit_1D(self):
        approx_order = 3
        (field, regions), mesh = prepare_field_1D(approx_order)
        u = FieldVariable('u', 'unknown', field, history=1)
        v = FieldVariable('v', 'test', field, primary_var_name='u')
        integral = Integral('i', order=approx_order * 2)
        a = Material('a', val=[1.0])

        term = AdvectDGFluxTerm("adv_stiff(a.val, u, v)",
                                "a.val, u, v",
                                integral, regions["omega"],
                                u=u, v=v, a=a)
        term.function([],  # out
                      [],  # state
                      "u",  # diff_var
                      field,
                      regions["omega"],
                      [],  # advelo
                      )


class TestDiffusionDGFluxTerm:

    def test_function_explicit_right_1D(self):
        approx_order = 3
        (field, regions), mesh = prepare_field_1D(approx_order)
        u = FieldVariable('u', 'unknown', field, history=1)
        v = FieldVariable('v', 'test', field, primary_var_name='u')
        integral = Integral('i', order=approx_order * 2)

        diffusion_tensor = .02
        D = Material('D', val=[diffusion_tensor])
        term = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)", "D.val, v,  u[-1]",
                                integral, regions["omega"], u=u, v=v, D=D)
        term.function([],  # out
                      [],  # state
                      None, # diff_var, explicit
                      field,
                      regions["omega"],
                      D, # advelo
                      )