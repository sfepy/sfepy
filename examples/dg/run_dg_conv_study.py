#!/usr/bin/env python
r"""
DG FEM convergence tests for 1D and 2D problems
"""
import sys
import os
from os.path import join as pjoin
import time
import numpy as nm
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
from sfepy.discrete import Integral, Problem
from run_dg_utils import clear_folder
from sfepy.discrete.fem.utils import refine_mesh as refine_mesh

# DG imports
from sfepy.discrete.dg.my_utils.visualizer import reconstruct_legendre_dofs

from examples.dg.run_dg_utils import  calculate_num_order, outputs_folder,\
    plot_conv_results, build_attrs_string, output, compute_erros, configure_output, \
    add_dg_arguments

from examples.dg.example_dg_common import get_1Dmesh_hook, diffusion_schemes_explicit


def create_argument_parser():
    parser = argparse.ArgumentParser(
        description='DG FEM convergence tests for 1D and 2D parametrized problems',
        epilog='(c) 2019  Man-machine Interaction at NTC UWB, ' +
            '\nauthor: Tomas Zitka, email: zitkat@ntc.zcu.cz')

    parser.add_argument("problem_file",
                        help="File containing specification for convergence test, "+
                             "must include define method and dim variable",
                        metavar='path')

    parser.add_argument("-o", "--output",
                        help="Output directory, can contain {}, " +
                             "exampe name is than plugged in.",
                        dest="output_dir",
                        metavar='path', type=str, default=None)

    parser.add_argument("-m", "--mesh", help="File with starting mesh to use",
                        default=None, metavar='path', action='store',
                        dest='mesh_file',)

    parser.add_argument("-dp", "--display-plots", help="Show interactive plots.",
                        default=False, action='store_true', dest='doplot',)

    parser.add_argument("-v", "--verbose", help="To be verbose or",
                        default=False, action='store_true', dest='verbose',)

    # parser.add_argument("--noscreenlog", help="Do not print log to screen",
    #                     default=False, action='store_true', dest='no_output_screen',)
    #
    # parser.add_argument('--logfile', type=str,
    #                     action='store', dest='output_log_name',
    #                     default="last_run.txt", help="Path to log output file")

    parser.add_argument('--orders', metavar="[1, 2, 3 ...]" ,default=None,
                        help='List of orders to try', type=str)

    parser.add_argument('--refines', metavar="[1, 2, 3 ...]", default=None,
                        help='List of refinement degrees', type=str)

    add_dg_arguments(parser)

    return parser


def parse_str2tuple_default(l, default=(1, 2, 3, 4)):
    if l is None:
        return default
    return tuple(int(x) for x in l.strip("()[],").split(","))


# soops-run -o .\soops_tests\ "mesh='mesh/mesh_tensr_2D_01_2.vtk', problem_file='advection/example_dg_advection2D', output_dir='soops_tests/{}/%s', --cfl=[0.1, .5], --limit=[@defined, @undefined],
# --orders='[1,2,3,4]'" .\conv_tests\run_dg_conv_study.py


# soops-run -o .\soops_tests\ "mesh='mesh/mesh_tensr_2D_01_2.vtk',
# problem_file='advection/example_dg_advection2D', output_dir='soops_tests/adv2D/%s',
# --adflux=[0.0, 0.5 ,1.0] ,--cfl=[0.01, 0.1, .5], --limit=[@defined, @undefined],
# --orders='[1,2,3,4]', --refines='[1,2,3,4]'" .\conv_tests\run_dg_conv_study.py

def get_run_info():
    """ Run info for soops """
    run_cmd = """
    python conv_tests/run_dg_conv_study.py {problem_file}
    --mesh={mesh} --output {output_dir}"""
    run_cmd = ' '.join(run_cmd.split())

    opt_args = {
        # advection parameters
        "--advelo": " --advelo={--advelo}",
        "--adflux": " --adflux={--adflux}",
        "--limit" : " --limit",

        # diffusion parameters
        "--cw"  : " --cw={--cw}",
        "--diffcoef"   : " --diffcoef={--diffcoef}",
        "--diffscheme" : " --diffscheme={--diffscheme}",

        # time parameteres
        "--cfl": " --cfl={--cfl}",
        "--dt" : " --dt={--dt}",

        # refining options
        "--orders"  : " --orders={--orders}",
        "--refines" : " --refines={--refines}",
        
        "--verbose" : " --verbose"
    }

    output_dir_key = "output_dir"
    is_finished_basename = ""

    return run_cmd, opt_args, output_dir_key, is_finished_basename


def main(argv):
    if argv is None:
        argv = sys.argv[1:]

    parser = create_argument_parser()
    args = parser.parse_args(argv)

    # configure_output({'output_screen': not args.no_output_screen,
    #                   'output_log_name': "last_run.txt"})

    problem_module_name = "examples.dg." + args.problem_file.replace(".py", "")\
        .replace("\\", ".").replace("/", ".")
    mesh = os.path.abspath(args.mesh_file)

    problem_module = importlib.import_module(problem_module_name)

    refines = parse_str2tuple_default(args.refines, (1, 2, 3, 4)[:2])
    orders = parse_str2tuple_default(args.orders, (1, 2, 3, 4)[:2])

    if problem_module.dim == 1:
        sol_fig, axs = plt.subplots(len(orders), len(refines), figsize=(18, 10))

    mod = sys.modules[__name__]
    results = []
    for ir, refine in enumerate(refines):
        gen_mesh = refine_mesh(mesh, refine)
        for io, order in enumerate(orders):
            conf = ProblemConf.from_dict(
                problem_module.define(
                    filename_mesh=gen_mesh,
                    approx_order=order,

                    flux=args.adflux,
                    limit=args.limit,

                    Cw=args.cw,
                    diffusion_coef=args.diffcoef,
                    diff_scheme_name=args.diffscheme,

                    CFL=args.cfl,
                    dt=args.dt,
                ), mod, verbose=args.verbose)

            output("----------------Running--------------------------")
            output("{}: {}".format(conf.example_name, time.asctime()))
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

            if args.output_dir is None:
                base_output_folder = pjoin(outputs_folder, "conv_tests_out",
                                           conf.example_name)
            elif "{}" in args.output_dir:
                base_output_folder = args.output_dir.format(conf.example_name)
            else:
                base_output_folder = pjoin(args.output_dir)

            output_folder = pjoin(base_output_folder, "h" + str(n_cells))
            output_folder = pjoin(output_folder, "o" + str(order))

            output_format = pjoin(output_folder, "sol-h{:02d}o{:02d}.*.{}"
                                  .format(n_cells, order,
                                          "vtk" if problem_module.dim == 1 else "msh"))
            output("Output set to {}, clearing.".format(output_format))


            clear_folder(output_format, confirm=False)
            ensure_path(output_format)

            pb.save_state(output_format.replace("*", "0"), state=pb.sol)
            output(
                "{}: {}".format(conf.example_name, time.asctime()))
            output("------------------Finished------------------\n\n")

            ana_l2, ana_qp, diff_l2, rel_l2, num_qp = compute_erros(conf.sol_fun, pb)

            n_dof = pb.fields["f"].n_nod

            result = (h, n_cells, nm.mean(vols), order, n_dof,
                      ana_l2, diff_l2, rel_l2, elapsed,
                      getattr(pb.ts_conf, "cour", nm.NAN),
                      getattr(pb.ts_conf, "dt", nm.NAN))

            results.append(result)

            if  problem_module.dim == 1:
                plot_1D_snr(conf, pb, ana_qp, num_qp,
                            io, order, orders, ir,
                            sol_fig, axs)
                sol_fig.savefig(pjoin(base_output_folder,
                                      "err-sol-i20" + build_attrs_string(conf)
                                      + ".png"), dpi=100)

    results = nm.array(results)


    err_df = pd.DataFrame(results,
                          columns=["h", "n_cells", "mean_vol", "order", "n_dof",
                                   "ana_l2", "diff_l2", "err_rel",
                                   "elapsed", "cour", "dt"])
    err_df = calculate_num_order(err_df)
    err_df.to_csv(pjoin(base_output_folder,
                        conf.example_name + "results" + build_attrs_string(conf)
                        + ".csv"))
    output(err_df)

    plot_conv_results(base_output_folder, conf, err_df, save=True)

    if args.doplot:
        plt.show()


def plot_1D_snr(conf, pb, ana_qp, num_qp, io, order, orders, ir, sol_fig, axs):
    """
    Plot 1D solutions and errors

    :param conf:
    :param io: index of order
    :param ir: inder of refirement
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
    ax.plot(xs, conf.analytic_sol(xs, t=nm.array(1)), label="exact", color="grey")
    ax.plot(xx[:, 0], uu[:, 0, 0], alpha=.5, label="num-lin")
    ax.plot(fqps, ana_qp.flatten(), "--", color="grey", label="exact-qp")
    ax.plot(fqps, num_qp.flatten(), label="num-qp")
    ax2 = ax.twinx()
    ax2.plot(fqps, nm.abs(num_qp.flatten() - ana_qp.flatten()), color="red",
             label="error-qp")

    if io < len(orders) - 1:
        ax.set_xticks([])
    else:
        if ir == 0:  # put only one legend
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            sol_fig.legend(lines + lines2, labels + labels2,
                           loc='center right')

if __name__ == '__main__':
    main(sys.argv[1:])
