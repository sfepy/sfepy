r"""
DG FEM convergence test for 2D problems
"""
import sys
import time
import numpy as nm
import sympy as sm
import pandas as pd

from os.path import join as pjoin

from matplotlib import pyplot as plt

from sfepy.base.ioutils import ensure_path
from sfepy.base.conf import ProblemConf
from sfepy.discrete import (Integral, Integrals, Material, Problem)
from sfepy.discrete.common.mappings import get_jacobian
from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder

from discrete.fem.utils import refine_mesh

from my_utils.visualizer import reconstruct_legendre_dofs

from examples.dg.example_dg_common import get_1Dmesh_hook, output,\
    get_gen_block_mesh_hook, calculate_num_order, build_attrs_string, plot_conv_results

mesh_center = (0.5, 0.5)
mesh_size = (1.0, 1.0)


# import problem definitions
# from examples.dg.diffusion.example_dg_diffusion2D_Hartmann import define
# from examples.dg.burgess.example_dg_kucera1 import define, mesh_center, mesh_size
# from examples.dg.diffusion.example_dg_laplace1 import define
from examples.dg.advection.example_dg_advection2D import define


def refine_square_tens(filename, refine):
    shape0 = (3, 3)  # always start with (3,3) for refine 0 it reduces to (2, 2), cahnge refine to shift refinement
    ashape = nm.array(shape0)
    shape = (ashape - 1) ** refine + 1

    gen_mesh = get_gen_block_mesh_hook(mesh_size, shape, mesh_center)
    return gen_mesh


def main():

    mesh = "mesh/mesh_tens_2D_01_20.vtk"

    refines = [0, 1, 2, 3, 4]
    orders = [1, 2, 3, 4, 5]
    # orders = [1]

    mod = sys.modules[__name__]

    results = []
    for ir, refine in enumerate(refines):

        gen_mesh = refine_mesh(mesh, refine)
        for io, order in enumerate(orders):

            conf = ProblemConf.from_dict(define(gen_mesh, order,
                                                Cw=10, diffusion_coef=0.002,
                                                # CFL=0.1,
                                                dt=1,
                                                ), mod)

            output("----------------------------------------------------")
            output(conf.example_name + ": " + time.asctime())
            output('n_nod:', refine, 'order:', order)

            try:
                conf.options.save_times = 0
            except AttributeError:
                pass
            pb = Problem.from_conf(conf)
            try:
                conf.options.pre_process_hook(pb)
            except AttributeError:
                pass
            tt = time.clock()
            sol = pb.solve()
            elapsed = time.clock() - tt

            n_cells = pb.domain.shape.n_el
            vols = pb.domain.cmesh.get_volumes(2)
            h = nm.mean(vols)
            if "2_3" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(4*vols))
            elif "2_4" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(2*vols))

            output('shape:', n_cells, 'order:', order)

            base_output_folder = pjoin("output", conf.example_name)
            output_folder = pjoin(base_output_folder, "h" + str(n_cells))
            output_folder = pjoin(output_folder, "o" + str(order))

            output_format = pjoin(output_folder, "sol-h{:02d}o{:02d}.*.{}".format(n_cells, order, "msh"))
            output("Output set to {}, clearing.".format(output_format))
            output("----------------------------------------------------\n\n")

            clear_folder(output_format, confirm=False)
            ensure_path(output_format)

            pb.save_state(output_format.replace("*", "0"), state=sol)

            idiff = Integral('idiff', 20)

            num_qp = pb.evaluate(
                'ev_volume_integrate.idiff.Omega(u)',
                integrals=Integrals([idiff]), mode='qp', copy_materials=False, verbose=False
            )

            aux = Material('aux', function=conf.sol_fun)
            ana_qp = pb.evaluate(
                'ev_volume_integrate_mat.idiff.Omega(aux.u, u)',
                aux=aux, integrals=Integrals([idiff]), mode='qp',
                copy_materials=False, verbose=False
            )

            field = pb.fields['f']
            det = get_jacobian(field, idiff)

            diff_l2 = nm.sqrt((((num_qp - ana_qp)**2) * det).sum())
            ana_l2 = nm.sqrt(((ana_qp ** 2) * det).sum())
            error = diff_l2 / ana_l2

            n_dof = field.n_nod

            result = (h, n_cells, nm.mean(vols), order, n_dof, ana_l2, diff_l2, error, elapsed,
                      pb.ts_conf.cour, pb.ts_conf.dt)

            results.append(result)

    results = nm.array(results)
    output(results)

    err_df = pd.DataFrame(results,
                          columns=["h", "n_cells", "mean_vol", "order", "n_dof", "ana_l2", "diff_l2", "err_rel",
                                   "elapsed", "cour", "dt"])
    err_df = calculate_num_order(err_df)
    err_df.to_csv(pjoin(base_output_folder, conf.example_name + "results" + build_attrs_string(conf) + ".csv"))

    plot_conv_results(base_output_folder, conf, err_df, save=True)

    plt.show()


if __name__ == '__main__':
    main()
