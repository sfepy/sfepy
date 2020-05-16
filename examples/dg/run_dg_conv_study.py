#!/usr/bin/env python
r"""
DG FEM convergence tests for 1D and 2D problems
"""
import sys
import os
from os.path import join as pjoin
from pathlib import Path
import time
import numpy as nm
import pandas as pd
import importlib
import argparse
import matplotlib

matplotlib.use("Qt5Agg")
from matplotlib import pyplot as plt


# sfepy imports
sys.path.append('.')
from sfepy.discrete.fem import Mesh
from sfepy.base.ioutils import ensure_path
from sfepy.base.conf import ProblemConf
from sfepy.discrete import Integral, Problem
from examples.dg.run_dg_utils import clear_folder, param_names
from sfepy.discrete.fem.utils import refine_mesh as refine_mesh
from sfepy.discrete.fem.meshio import GmshIO

# DG imports
from sfepy.discrete.dg.my_utils.visualizer import reconstruct_legendre_dofs

from examples.dg.run_dg_utils import outputs_folder,\
    plot_conv_results, build_attrs_string, output, compute_erros, configure_output, \
    add_dg_arguments


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

    parser.add_argument("-nr", "--dot-not-recalculate",
                        help="Force racalcualting, note that loading previous" +
                             " results destroys detailed information about run",
                        default=False, action='store_true', dest='no_recalc', )

    parser.add_argument("-v", "--verbose", help="To be verbose or",
                        default=False, action='store_true', dest='verbose',)

    parser.add_argument("--noscreenlog", help="Do not print log to screen",
                        default=False, action='store_true', dest='no_output_screen',)
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


def get_run_info():
    """ Run info for soops """
    run_cmd = """
    python run_dg_conv_study.py {problem_file}
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
    is_finished_basename = "results.csv"

    return run_cmd, opt_args, output_dir_key, is_finished_basename


def main(argv):
    if argv is None:
        argv = sys.argv[1:]

    parser = create_argument_parser()
    args = parser.parse_args(argv)

    prefix = "examples.dg."

    problem_module_name = args.problem_file.replace(".py", "").strip("\\.") \
        .replace("\\", ".").replace("/", ".")
    if not problem_module_name.startswith(prefix):
        problem_module_name = prefix + problem_module_name

    problem_module = importlib.import_module(problem_module_name)

    mesh = str(Path(args.mesh_file))

    refines = parse_str2tuple_default(args.refines, (1, 2, 3, 4))
    orders = parse_str2tuple_default(args.orders, (0, 1, 2, 3, 4))

    if problem_module.dim == 1:
        sol_fig, axs = plt.subplots(len(orders), len(refines), figsize=(18, 10))

    mod = sys.modules[problem_module_name]
    results = []
    for ir, refine in enumerate(refines):
        gen_mesh = refine_mesh(mesh, refine)
        for io, order in enumerate(orders):
            conf = ProblemConf.from_dict(
                problem_module.define(
                    filename_mesh=gen_mesh,
                    approx_order=order,

                    adflux=args.adflux,
                    limit=args.limit,

                    cw=args.cw,
                    diffcoef=args.diffcoef,
                    diffscheme=args.diffscheme,

                    cfl=args.cfl,
                    dt=args.dt,
                ), mod, verbose=args.verbose)
            conf.options.absolute_mesh_path = True

            if args.output_dir is None:
                base_output_folder = Path(outputs_folder) / "conv_tests_out" /\
                                           conf.example_name
            elif "{}" in args.output_dir:
                base_output_folder = Path(args.output_dir.format(conf.example_name))
            else:
                base_output_folder = Path(args.output_dir)

            output_folder = base_output_folder / ("h" + str(refine))
            output_folder = output_folder / ("o" + str(order))

            configure_output({'output_screen': not args.no_output_screen,
                              'output_log_name': str(output_folder / "last_run.txt")})

            output("----------------Running--------------------------")
            output("{}: {}".format(conf.example_name, time.asctime()))
            output('refine:', refine, 'order:', order)

            h, n_cells, pb, vols = create_problem(conf)

            output_format = pjoin(str(output_folder), "sol-h{:02d}o{:02d}.*.{}"
                                  .format(n_cells, order,
                                          "vtk" if problem_module.dim == 1 else "msh"))
            output("Output set to {}, clearing.".format(output_format))

            clear_folder(output_format, confirm=False, doit=True)
            ensure_path(output_format)

            pb, elapsed = run_calc(pb, output_format)

            output("{}: {}".format(conf.example_name, time.asctime()))
            output("------------------Finished------------------\n\n")

            ana_l2, ana_qp, diff_l2, rel_l2, num_qp = compute_erros(conf.sol_fun, pb)

            n_dof = pb.fields["f"].n_nod

            result = (h, n_cells, nm.mean(vols), order, n_dof,
                      ana_l2, diff_l2, rel_l2, elapsed,
                      getattr(pb.ts_conf, "cour", nm.NAN),
                      getattr(pb.ts_conf, "dt", nm.NAN),
                      getattr(pb.solver.status.nls_status, "err", nm.NAN),
                      getattr(pb.solver.status.nls_status, "n_iter", nm.NAN)
                      )

            results.append(result)

            if problem_module.dim == 1:
                plot_1D_snr(conf, pb, ana_qp, num_qp,
                            io, order, orders, ir,
                            sol_fig, axs)
                sol_fig.savefig(base_output_folder /
                                ("err-sol-i20" + build_attrs_string(conf) + ".png"),
                                dpi=100)

    err_df = create_error_df(conf, pb, results)

    err_df.to_csv(base_output_folder / "results.csv")

    err_df.to_csv(base_output_folder.parent /
                  (base_output_folder.name + "_results.csv"))

    output(err_df)

    conv_fig = plot_conv_results(base_output_folder, conf, err_df, save=True)
    conv_fig.savefig(base_output_folder / "results.png", dpi=200)

    if args.doplot:
        plt.show()


def create_error_df(conf, pb, results):
    results = nm.array(results)
    err_df = pd.DataFrame(results,
                          columns=["h", "n_cells", "mean_vol", "order", "n_dof",
                                   "ana_l2", "diff_l2", "err_rel",
                                   "elapsed", "cour", "actual_dt",
                                   "nls_error", "nls_iter"])
    err_df = calculate_num_order(err_df)
    for name in param_names:
        err_df[name] = conf.__dict__[name]
    err_df["gel"] = pb.domain.mesh.descs[0]
    return err_df


def create_problem(conf):
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
    return h, n_cells, pb, vols


def run_calc(pb, output_format):
    tt = time.process_time()
    pb.sol = pb.solve()
    elapsed = time.process_time() - tt
    pb.save_state(output_format.replace("*", "0"), state=pb.sol)
    return pb, elapsed


def calculate_num_order(err_df):
    """
    Uses diff_l2 and n_rows columns of the dataframe to calculate num_order,
    splits dataframe on order clumn
    :param err_df: dataframe, columns: ["n_rows", "order", diff_l2]
    :return:
    """
    res_df = pd.DataFrame()
    for order in err_df["order"].unique():
        order_err_df = err_df[err_df["order"] == order].sort_values("n_cells")

        num_orders = [nm.NAN]

        last_err = order_err_df.iloc[0]["diff_l2"]
        last_h = order_err_df.iloc[0]["h"]
        #         print(order_err_df.iloc[1:, :])
        for i, row in order_err_df.iloc[1:, :].iterrows():
            num_order = nm.log(row["diff_l2"] / last_err) \
                        / nm.log(row["h"] / last_h)
            #             print(row["err_l2"] / last_err)
            #             print(row["n_rows] / last_h)
            #             print("-------------------")

            last_err = row["diff_l2"]
            last_h = row["h"]
            num_orders.append(num_order)

        order_err_df["num_order"] = num_orders
        res_df = res_df.append(order_err_df)
    return res_df


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
