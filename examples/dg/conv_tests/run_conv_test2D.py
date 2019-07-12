r"""
DG FEM convergence test for 1D diffusion only problem
"""
import numpy as nm
import sympy as sm


from examples.dg.example_dg_common import get_1Dmesh_hook, get_gen_block_mesh_hook
from examples.dg.conv_tests.fem_conv_test import SimpleExpression, define_functions, eval_expr
from os.path import join as pjoin

from my_utils.visualizer import reconstruct_legendre_dofs

# from examples.dg.diffusion.example_dg_diffusion2D_Hartmann import define
# from examples.dg.diffusion.example_dg_diffusion2D_Hartmann import define
# from examples.dg.burgess.example_dg_kucera1 import define
from examples.dg.diffusion.example_dg_laplace1 import define

def main():
    import sys
    import time
    from sfepy.base.ioutils import ensure_path
    from sfepy.base.base import output
    from sfepy.base.conf import ProblemConf
    from sfepy.discrete import (Integral, Integrals, Material, Problem)
    from sfepy.discrete.common.mappings import get_jacobian
    from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder
    from matplotlib import pyplot as plt

    shape0 = (3, 3)
    n_refine = 6
    orders = [1, 2, 3, 4, 5]
    # orders = [1]

    mod = sys.modules[__name__]

    # sol_fig, axs = plt.subplots(len(orders), n_refine, figsize=(18, 10))

    results = []
    for ir, refine in enumerate(range(n_refine)):
        ashape = nm.array(shape0)
        shape = (ashape - 1) ** refine + 1

        gen_mesh = get_gen_block_mesh_hook((1., 1.), shape, (.5, .5))
        for io, order in enumerate(orders):
            output('shape:', shape, 'order:', order)

            conf = ProblemConf.from_dict(define(gen_mesh, order, Cw=1000, diffusion_coef=1), mod)
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

            n_cell = shape[0] - 1

            base_output_folder = pjoin("output", conf.example_name)
            output_folder = pjoin(base_output_folder, "h" + str(n_cell))
            output_folder = pjoin(output_folder, "o" + str(order))

            output_format = pjoin(output_folder, "sol-h{:02d}o{:02d}.*.{}".format(n_cell, order, "msh"))
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
            # sol_fig.suptitle(
            #     "Numerical and exact solutions, Cw: {}, diffusion: {}".format(conf.Cw, conf.diffusion_coef))
            #
            # qps = pb.fields["f"].mapping.get_physical_qps(idiff.get_qp("1_2")[0])
            # fqps = qps.flatten()
            # coors = pb.domain.mesh.coors
            # u = pb.fields["f"].unravel_sol(sol.vec)
            # uu, xx = reconstruct_legendre_dofs(coors, None, u.swapaxes(0, 1)[:, :, None])
            #
            # ax = axs[order-1][refine]
            #
            # xs = nm.linspace(nm.min(0), nm.max(1), 500)[:, None]
            # ax.set_title("o: {}, h: {}".format(order, n_nod - 1))
            # ax.plot(xs, conf.analytic_sol(xs), label="fun-ex", color="grey")
            # ax.plot(xx[:, 0], uu[:, 0, 0], alpha=.5)
            # ax.plot(fqps, ana_qp.flatten(), "--", color="grey")
            # ax.plot(fqps, num_qp.flatten())
            # ax2 = ax.twinx()
            # ax2.plot(fqps, nm.abs(num_qp.flatten() - ana_qp.flatten()), color="red")
            # if io < len(orders) - 1:
            #     ax.set_xticks([])
            # # if ir > 0:
            # #     ax.set_yticks([])

            result = (n_cell, order, n_dof, ana_l2, diff_l2, error, elapsed)
            results.append(result)


    # sol_fig.savefig("err-sol-i20cw{}_d{}.tif".format(conf.Cw, conf.diffusion_coef), dpi=100)
    results = nm.array(results)
    output(results)

    conv_fig = plt.figure()
    conv_fig.suptitle("Convergences by order, Cw: {}, diffusion: {}".format(conf.Cw, conf.diffusion_coef))
    for o in orders:
        curr_results = results[results[:, 1] == o]
        co = plt.loglog(curr_results[:, 0]**2, curr_results[:, 4], 'o', label=str(o))[0].get_color()
        plt.loglog(curr_results[:, 0]**2, curr_results[:, 4], color=co)
        plt.grid()
        plt.xlabel("cells")
        plt.ylabel("L^2 error")
    plt.legend()
    conv_fig.savefig(pjoin(base_output_folder, conf.example_name + "-cells-cw{}_d{}.tif".format(conf.Cw, conf.diffusion_coef)), dpi=200)


    import pandas as pd
    err_df = pd.DataFrame(results,
                          columns=["h_coef", "order", "n_dof", "ana_l2", "diff_l2", "err_rel", "elapsed"])
    err_df.to_csv(pjoin(base_output_folder, conf.example_name + "results-cw{}_d{}.csv".format(conf.Cw, conf.diffusion_coef)))


    plt.show()

if __name__ == '__main__':
    main()
