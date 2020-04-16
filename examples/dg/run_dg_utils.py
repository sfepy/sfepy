"""
Script with utility functions for running DG examples and convergence studies
"""
from glob import glob
import os

import pandas as pd

from matplotlib import pyplot as plt
import numpy as nm
from os.path import join as pjoin


from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.discrete import Integral, Material, Integrals
from sfepy.discrete.common.mappings import get_jacobian
from sfepy.base.base import (get_default, output, configure_output, assert_,
                             Struct, basestr, IndexedStruct)

outputs_folder = "outputs"

configure_output({'output_screen': True,
                  'output_log_name': pjoin(outputs_folder, "last_run.txt")})


def plot_conv_results(base_output_folder, conf, err_df,
                      plot_title_attrs=None, save=False):
    """

    :param base_output_folder:
    :param conf: conf structure defined in declarative mode, or object
        containing: dt, CFL or diffusion_coef attributese, example_name
    :plot_title_attrs: attributes to list in the figure sup title
    :param err_df: pandas dataframe containing at least:

    "order", "num_order", "n_cells", "diff_l2" columns

    :return:
    """
    if plot_title_attrs is None:
        plot_title_attrs = ["Cw", "diffusion_coef", "dt", "CFL"]
    fig_sup_title = "Convergence by order"
    file_name = conf.example_name + "-cells"
    fig_sup_title += build_attrs_string(conf, plot_title_attrs, sep=", ",
                                        remove_dots=False)
    file_name += build_attrs_string(conf, plot_title_attrs)

    conv_fig, ax = plt.subplots(1, 1)
    conv_fig.suptitle(fig_sup_title)

    orders = sorted(err_df["order"].unique())
    for o in orders:
        ax.set_xscale('log', basex=2)
        ax.set_yscale('log', basey=10)
        curr_df = err_df[err_df["order"] == o]
        co = ax.plot(curr_df["n_cells"], curr_df["diff_l2"], 'o',
                       label=str(int(o)))[0].get_color()

        ax.plot(curr_df["n_cells"], curr_df["diff_l2"], color=co, label="")
        for i, r in curr_df.iloc[1:, :].iterrows():
            ax.text(r["n_cells"], r["diff_l2"], "{:.2}".format(r["num_order"]))

        ax.grid(True)
        ax.set_xlabel("cells")
        ax.set_ylabel("L^2 error")
    ax.legend(title="Order")

    if save:
        conv_fig.savefig(
            pjoin(base_output_folder,
                  file_name + ".png"),
            dpi=200)

    return conv_fig

def compute_erros(analytic_fun, pb):
    """
    Compute errors from analytical solution in conf.sol_fun and numerical
    solution saved in pb
    :param analytic_fun: analytic solution
    :param pb: problem to with numerical solution
    :return: analytic_fun L2 norm,
              vaules of analytic_fun in qps
              L2 norm of difference between analytic and numerical solution
              relative difference
              values of numerical solution in qps
    """
    idiff = Integral('idiff', 20)
    num_qp = pb.evaluate(
        'ev_volume_integrate.idiff.Omega(u)',
        integrals=Integrals([idiff]), mode='qp',
        copy_materials=False, verbose=False
    )
    aux = Material('aux', function=analytic_fun)
    ana_qp = pb.evaluate(
        'ev_volume_integrate_mat.idiff.Omega(aux.u, u)',
        aux=aux, integrals=Integrals([idiff]), mode='qp',
        copy_materials=False, verbose=False
    )
    field = pb.fields['f']
    det = get_jacobian(field, idiff)
    diff_l2 = nm.sqrt((((num_qp - ana_qp) ** 2) * det).sum())
    ana_l2 = nm.sqrt(((ana_qp ** 2) * det).sum())
    rel_l2 = diff_l2 / ana_l2
    return ana_l2, ana_qp, diff_l2, rel_l2, num_qp


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


def build_attrs_string(conf, attrs=("Cw", "diffusion_coef", "dt", "CFL"),
                       sep="_", ret_form=False, remove_dots=True):
    """
    Builds string from intersection of attributes
    list attrs and conf attributes, removes "." !
    :param conf:
    :param attrs:
    :param sep:
    :return:
    """
    form = ""
    attr_vals = []
    for attr in attrs:
        attr_val = getattr(conf, attr, None)
        if attr_val is not None:
            form += sep + attr + "{}"
            attr_vals += [attr_val]
    if not remove_dots:
        res = form.format(*attr_vals)
    else:
        res = form.format(*attr_vals).replace(".", "")

    if ret_form:
        return res, form, attr_vals
    else:
        return res


def clear_folder(clear_format, confirm=False):
    """
    Deletes files matching the format
    :param clear_format:
    :param confirm:
    :return:
    """
    files = glob(clear_format)
    doit = True
    if confirm:
        for file in files:
            output("Will delete file {}".format(file))
        doit = input("--------------\nDelete files [Y/n]? ").strip() == "Y"

    if doit:
        for file in files:
            os.remove(file)