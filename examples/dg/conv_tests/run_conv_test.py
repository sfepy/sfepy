r"""
DG FEM convergence tests for 1D and 2D problems
"""
import sys
import os
from os.path import join as pjoin
import time
import numpy as nm
import sympy as sm
import pandas as pd
import importlib
import argparse
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib import pyplot as plt



# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.base.ioutils import ensure_path
from sfepy.base.conf import ProblemConf
from sfepy.discrete import (Integral, Integrals, Material, Problem)
from sfepy.discrete.common.mappings import get_jacobian
from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder
from sfepy.discrete.fem.utils import refine_mesh as refine_mesh

# DG imports
from sfepy.discrete.dg.my_utils.visualizer import reconstruct_legendre_dofs

from examples.dg.example_dg_common import get_1Dmesh_hook, calculate_num_order, \
    plot_conv_results, build_attrs_string, output


parser = argparse.ArgumentParser(description='DG FEM convergence tests for 1D and 2D problems',
                                 epilog='(c) 2019  Man-machine Interaction at NTC UWB, ' +
                                        '\nauthor: Tomas Zitka, email: zitkat@ntc.zcu.cz')
parser.add_argument("problem_file",
                    help="File containing specification for convergence test, must include "
                         "define method and dim variable",
                    metavar='path')

parser.add_argument("-m", "--mesh", help="File with starting mesh to use",
                    default=None, metavar='path', action='store', dest='mesh_file',)

parser.add_argument("-ps", "--plot-solutions", help="Plot 1D solutions and errors",
                    default=False, action='store_true', dest='do1Dplot',)


def main(argv):
    if argv is None:
        argv = sys.argv[1:]
    args = parser.parse_args(argv)

    problem_module_name = "examples.dg." + args.problem_file.replace(".py", "")\
        .replace("\\", ".").replace("/", ".")
    mesh = os.path.abspath(args.mesh_file)

    problem_module = importlib.import_module(problem_module_name)

    refines = [1, 2, 3, 4]
    orders = [1, 2, 3, 4]

    if args.do1Dplot and problem_module.dim == 1:
        sol_fig, axs = plt.subplots(len(orders), len(refines), figsize=(18, 10))

    mod = sys.modules[__name__]
    results = []
    for ir, refine in enumerate(refines):
        gen_mesh = refine_mesh(mesh, refine)
        for io, order in enumerate(orders):
            conf = ProblemConf.from_dict(
                problem_module.define(
                    gen_mesh, order,
                    # Cw=10,
                    # diffusion_coef=0.002,
                    # CFL=0.1,
                    # dt=1,
                ), mod, verbose=False)

            output("----------------Running--------------------------")
            output(conf.example_name + ": " + time.asctime())
            output('refine:', refine, 'order:', order)

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
            vols = pb.domain.cmesh.get_volumes(1)
            h = nm.mean(vols)
            if "2_3" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(4 * vols))
            elif "2_4" in pb.domain.geom_els:
                h = nm.mean(nm.sqrt(2 * vols))

            tt = time.clock()
            pb.sol = pb.solve()
            elapsed = time.clock() - tt

            base_output_folder = pjoin("conv_tests_output", conf.example_name)
            output_folder = pjoin(base_output_folder, "h" + str(n_cells))
            output_folder = pjoin(output_folder, "o" + str(order))

            output_format = pjoin(output_folder, "sol-h{:02d}o{:02d}.*.{}"
                                  .format(n_cells, order,
                                          "vtk" if problem_module.dim == 1 else "msh"))
            output("Output set to {}, clearing.".format(output_format))
            output("------------------Finished------------------\n\n")

            clear_folder(output_format, confirm=False)
            ensure_path(output_format)

            pb.save_state(output_format.replace("*", "0"), state=pb.sol)

            ana_l2, ana_qp, diff_l2, error, num_qp = compute_erros(conf, pb)

            n_dof = pb.fields["f"].n_nod

            result = (h, n_cells, nm.mean(vols), order, n_dof,
                      ana_l2, diff_l2, error, elapsed,
                      getattr(pb.ts_conf, "cour", nm.NAN),
                      getattr(pb.ts_conf, "dt", nm.NAN))

            results.append(result)

            if args.do1Dplot and problem_module.dim == 1:
                plot_1D_snr(conf, pb, ana_qp, num_qp,
                            io, order, orders, ir,
                            sol_fig, axs)
                sol_fig.savefig(pjoin(base_output_folder,
                                      "err-sol-i20" + build_attrs_string(
                                          conf) + ".png"), dpi=100)

    results = nm.array(results)
    output(results)

    err_df = pd.DataFrame(results,
                          columns=["h", "n_cells", "mean_vol", "order", "n_dof",
                                   "ana_l2", "diff_l2", "err_rel",
                                   "elapsed", "cour", "dt"])
    err_df = calculate_num_order(err_df)
    err_df.to_csv(pjoin(base_output_folder,
                        conf.example_name + "results" + build_attrs_string(conf) + ".csv"))

    plot_conv_results(base_output_folder, conf, err_df, save=True)

    plt.show()


def compute_erros(conf, pb):
    """
    Compute errors from analytical solution in conf.sol_fun and numerical
    solution saved in pb
    :param conf:
    :param pb:
    :return:
    """
    idiff = Integral('idiff', 20)
    num_qp = pb.evaluate(
        'ev_volume_integrate.idiff.Omega(u)',
        integrals=Integrals([idiff]), mode='qp',
        copy_materials=False, verbose=False
    )
    aux = Material('aux', function=conf.sol_fun)
    ana_qp = pb.evaluate(
        'ev_volume_integrate_mat.idiff.Omega(aux.u, u)',
        aux=aux, integrals=Integrals([idiff]), mode='qp',
        copy_materials=False, verbose=False
    )
    field = pb.fields['f']
    det = get_jacobian(field, idiff)
    diff_l2 = nm.sqrt((((num_qp - ana_qp) ** 2) * det).sum())
    ana_l2 = nm.sqrt(((ana_qp ** 2) * det).sum())
    error = diff_l2 / ana_l2
    return ana_l2, ana_qp, diff_l2, error, num_qp


def plot_1D_snr(conf, pb, ana_qp, num_qp, io, order, orders, ir, sol_fig, axs):
    """
    Plot 1D solutions and errors

    :param conf:
    :param io:
    :param ana_qp:
    :param num_qp:
    :param order:
    :param orders:
    :param pb:

    :param sol_fig:
    :param axs:
    :return:
    """
    idiff = Integral('idiff', 20)
    sol_fig.suptitle(
        "Numerical and exact solutions" +
        build_attrs_string(conf, remove_dots=False, sep=", "))
    n_cells = pb.domain.shape.n_el

    qps = pb.fields["f"].mapping.get_physical_qps(idiff.get_qp("1_2")[0])
    fqps = qps.flatten()
    coors = pb.domain.mesh.coors
    u = pb.fields["f"].unravel_sol(pb.sol.vec)
    uu, xx = reconstruct_legendre_dofs(coors, None, u.swapaxes(0, 1)[:, :, None])
    ax = axs[io][ir]
    xs = nm.linspace(conf.mstart, conf.mend, 500)[:, None]
    ax.set_title("o: {}, h: {}".format(order, n_cells))
    ax.plot(xs, conf.analytic_sol(xs, t=nm.array(1)), label="fun-ex", color="grey")
    ax.plot(xx[:, 0], uu[:, 0, 0], alpha=.5)
    ax.plot(fqps, ana_qp.flatten(), "--", color="grey")
    ax.plot(fqps, num_qp.flatten())
    ax2 = ax.twinx()
    ax2.plot(fqps, nm.abs(num_qp.flatten() - ana_qp.flatten()), color="red")
    if io < len(orders) - 1:
        ax.set_xticks([])
    # if ir > 0:
    #     ax.set_yticks([])


if __name__ == '__main__':
    main(sys.argv[1:])
