r"""
DG FEM convergence test for 1D diffusion only problem
"""
import numpy as nm
import sympy as sm
import pandas as pd

from matplotlib import pyplot as plt

from examples.dg.example_dg_common import get_1Dmesh_hook, \
    get_gen_block_mesh_hook, calculate_num_order, build_attrs_string, plot_conv_results
from os.path import join as pjoin

mesh_center = (0.5, 0.5)
mesh_size = (1.0, 1.0)

from discrete.fem.utils import refine_mesh

from my_utils.visualizer import reconstruct_legendre_dofs

# from examples.dg.diffusion.example_dg_diffusion2D_Hartmann import define
from examples.dg.burgess.example_dg_kucera1 import define, mesh_center, mesh_size
# from examples.dg.diffusion.example_dg_laplace1 import define
# from examples.dg.advection.example_dg_advection2D import define


def refine_square_tens(filename, refine):
    shape0 = (3, 3)  # always start with (3,3) for refine 0 it reduces to (2, 2), cahnge refine to shift refinement
    ashape = nm.array(shape0)
    shape = (ashape - 1) ** refine + 1

    gen_mesh = get_gen_block_mesh_hook(mesh_size, shape, mesh_center)
    return gen_mesh


def main():
    import sys
    import time
    from sfepy.base.ioutils import ensure_path
    from sfepy.base.base import output
    from sfepy.base.conf import ProblemConf
    from sfepy.discrete import (Integral, Integrals, Material, Problem)
    from sfepy.discrete.common.mappings import get_jacobian
    from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder

    mesh = "mesh/mesh_simp_2D_11_8.vtk"

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
                                                dt=1e-5,
                                                ), mod)
            try:
                conf.options.save_times = 0
            except AttributeError:
                pass
            pb = Problem.from_conf(conf)
            try:
                conf.options.pre_process_hook(pb)
            except AttributeError:
                pass

            n_cells = pb.domain.shape.n_el
            vols = pb.domain.cmesh.get_volumes(2)
            if "2_3" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(4*vols))
            elif "2_4" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(2*vols))

            output('shape:', n_cells, 'order:', order)

            tt = time.clock()
            sol = pb.solve()
            elapsed = time.clock() - tt

            base_output_folder = pjoin("output", conf.example_name)
            output_folder = pjoin(base_output_folder, "h" + str(n_cells))
            output_folder = pjoin(output_folder, "o" + str(order))

            output_format = pjoin(output_folder, "sol-h{:02d}o{:02d}.*.{}".format(n_cells, order, "msh"))
            output("Output set to {}, clearing ...".format(output_format))

            clear_folder(output_format, confirm=False)
            ensure_path(output_format)

            pb.save_state(output_format.replace("*", "0"), state=sol)

            idiff = Integral('idiff', 20)

            num_qp = pb.evaluate(
                'ev_volume_integrate.idiff.Omega(u)',
                integrals=Integrals([idiff]), mode='qp', copy_materials=False,
            )

            aux = Material('aux', function=conf.sol_fun)
            ana_qp = pb.evaluate(
                'ev_volume_integrate_mat.idiff.Omega(aux.u, u)',
                aux=aux, integrals=Integrals([idiff]), mode='qp',
                copy_materials=False,
            )

            field = pb.fields['f']
            det = get_jacobian(field, idiff)

            diff_l2 = nm.sqrt((((num_qp - ana_qp)**2) * det).sum())
            ana_l2 = nm.sqrt((((ana_qp)**2) * det).sum())
            error = diff_l2 / ana_l2

            n_dof = field.n_nod

            result = (h, n_cells, nm.mean(vols), order, n_dof, ana_l2, diff_l2, error, elapsed)
            results.append(result)

    results = nm.array(results)

    err_df = pd.DataFrame(results,
                          columns=["h", "n_cells", "mean_vol", "order", "n_dof", "ana_l2", "diff_l2", "err_rel", "elapsed"])
    err_df = calculate_num_order(err_df)
    err_df.to_csv(
        pjoin(base_output_folder, conf.example_name + "results-cw{}_d{}.csv".format(conf.Cw, conf.diffusion_coef)))

    plot_conv_results(base_output_folder, conf, err_df)


    plt.show()


if __name__ == '__main__':
    main()
