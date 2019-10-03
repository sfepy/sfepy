r"""
DG FEM convergence test for 1D diffusion only problem
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
from matplotlib import pyplot as plt

from my_utils.visualizer import reconstruct_legendre_dofs


from examples.dg.example_dg_common import get_1Dmesh_hook, calculate_num_order, plot_conv_results, build_attrs_string, output
from examples.dg.conv_tests.fem_conv_test import SimpleExpression, define_functions, eval_expr

mstart = -1
mend = 1

# import problem definitions
# from examples.dg.burgess.example_dg_burgess1D_Hesthaven import define
# from examples.dg.diffusion.example_dg_diffusion1D import define
from examples.dg.advection.example_dg_advection1D import define, mstart, mend

def main():


    n_nod = 5
    n_refine = 8
    orders = [1, 2, 3, 4, 5]
    orders = [1, 2]

    mod = sys.modules[__name__]

    sol_fig, axs = plt.subplots(len(orders), n_refine, figsize=(18, 10))

    results = []
    for ir, refine in enumerate(range(n_refine)):
        gen_mesh = get_1Dmesh_hook(mstart, mend, n_nod)
        for io, order in enumerate(orders):

            conf = ProblemConf.from_dict(
                define(gen_mesh, order,
                       dt=1,
                       # CFL=0.002
                )
                , mod, verbose=False)

            output("----------------------------------------------------")
            output(conf.example_name + ": " + time.asctime())
            output('n_nod:', n_nod, 'order:', order)

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
            pb.sol = pb.solve()
            elapsed = time.clock() - tt

            n_cells = pb.domain.shape.n_el
            vols = pb.domain.cmesh.get_volumes(1)
            h = nm.mean(vols)
            if "2_3" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(4*vols))
            elif "2_4" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(2*vols))

            output('shape:', n_cells, 'order:', order)

            base_output_folder = pjoin("output", conf.example_name)
            output_folder = pjoin(base_output_folder, "h" + str(n_cells))
            output_folder = pjoin(output_folder, "o" + str(order))

            output_format = pjoin(output_folder, "sol-h{:02d}o{:02d}.*.{}".format(n_nod, order, "vtk"))
            output("Output set to {}, clearing.".format(output_format))
            output("----------------------------------------------------\n\n")

            clear_folder(output_format, confirm=False)
            ensure_path(output_format)

            pb.save_state(output_format.replace("*", "0"), state=pb.sol)

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

            plot_1D_snr(conf, pb, ana_qp, num_qp, io, order, orders, refine, n_nod, sol_fig, axs)

        n_nod = n_nod + n_nod - 1

    sol_fig.savefig(pjoin(base_output_folder, "err-sol-i20" + build_attrs_string(conf) + ".png"), dpi=100)

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
