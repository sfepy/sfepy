import numpy as nm
import pandas as pd
from matplotlib import pyplot as plt

from os.path import join as pjoin


from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO

from sfepy.base.base import (get_default, output, configure_output, assert_,
                             Struct, basestr, IndexedStruct)

from sfepy.terms import register_term
from sfepy.solvers import register_solver

# import various ICs
from sfepy.discrete.dg.my_utils.inits_consts import ghump, gsmooth, \
    left_par_q, left_cos, superic, three_step_u, sawtooth_q, const_q, quadr_cub


from sfepy.discrete.dg.dg_tssolver import TVDRK3StepSolver, RK4StepSolver, \
    EulerStepSolver

from sfepy.discrete.dg.dg_limiters import IdentityLimiter, MomentLimiter1D

from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm, \
    NonlinScalarDotGradTerm, NonlinearHyperDGFluxTerm
from sfepy.discrete.dg.dg_terms import DiffusionDGFluxTerm, \
    DiffusionInteriorPenaltyTerm


configure_output({'output_screen': True, 'output_log_name': "last_run.txt"})

register_term(AdvectDGFluxTerm)
register_term(NonlinScalarDotGradTerm)
register_term(NonlinearHyperDGFluxTerm)
register_term(DiffusionDGFluxTerm)
register_term(DiffusionInteriorPenaltyTerm)

register_solver(TVDRK3StepSolver)
register_solver(RK4StepSolver)
register_solver(EulerStepSolver)

functions = {}

diffusion_schemes_implicit = {"symmetric" :
                                  "- dw_dg_diffusion_flux.i.Omega(D.val, u, v)"
                                + "- dw_dg_diffusion_flux.i.Omega(D.val, v, u)",
                     "non-symmetric":
                         "- dw_dg_diffusion_flux.i.Omega(D.val, u, v)"
                       + "+ dw_dg_diffusion_flux.i.Omega(D.val, v, u)",
                     "incomplete": " dw_dg_diffusion_flux.i.Omega(D.val, u, v)"}

diffusion_schemes_explicit = {"symmetric" :
                              "- dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)"
                            + "- dw_dg_diffusion_flux.i.Omega(D.val, v, u[-1])",
                         "non-symmetric":
                             "- dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)"
                           + "+ dw_dg_diffusion_flux.i.Omega(D.val, v, u[-1])",
                     "incomplete":
                            " dw_dg_diffusion_flux.i.Omega(D.val, u[-1], v)"}


def local_register_function(fun):
    try:
        functions.update({fun.__name__: (fun,)})

    except AttributeError:  # Already a sfepy Function.
        fun = fun.function
        functions.update({fun.__name__: (fun,)})

    return fun


def get_cfl_setup(CFL=None, dt=None):
    """
    Provide either CFL or dt to create preprocess hook that sets up


    Params
    ------
    CFL : float, optional
    dt: float, optional


    Returns
    -------
    setup_cfl_condition(sfepy.discrete.problem)

    """

    if CFL is None and dt is None:
        raise ValueError("Specifiy either CFL or dt in CFL setup")

    def setup_cfl_condition(problem):
        """
        Sets up CFL condition for problem ts_conf in problem
        :param problem: discrete.problem.Problem
        :return:
        """
        ts_conf = problem.ts_conf
        mesh = problem.domain.mesh
        dim = mesh.dim
        first_field = list(problem.fields.values())[0]
        first_field_name = list(problem.fields.keys())[0]
        approx_order = first_field.approx_order
        mats = problem.create_materials('a')
        try:
            # TODO make this more general,
            #  maybe require velocity material name in parameters
            velo = problem.conf_materials['material_a__0'].values["val"]
            max_velo = nm.max(nm.linalg.norm(velo))
        except KeyError:
            max_velo = 1
        dx = nm.min(problem.domain.mesh.cmesh.get_volumes(dim))
        order_corr = 1. / (2 * approx_order + 1)

        if dt is None:
            _dt = dx / max_velo * CFL * order_corr
            if not(nm.isfinite(_dt)):
                _dt = 1
        else:
            _dt = dt
        # time_steps_N = int((tf - t0) / dt) * 2
        tn = int(nm.ceil((ts_conf.t1 - ts_conf.t0) / _dt))
        dtdx = _dt / dx

        ts_conf.dt = _dt
        ts_conf.n_step = tn
        ts_conf.cour = max_velo * dtdx

        output("Preprocess hook - setup_cfl_condition:...")
        output("Approximation order of field {}({}) is {}"
               .format(first_field_name, first_field.family_name, approx_order))
        output("Space divided into {0} cells, {1} steps, step size is {2}"
               .format(mesh.n_el, len(mesh.coors), dx))
        output("Time divided into {0} nodes, {1} steps, step size is {2}"
               .format(tn - 1, tn, _dt))
        if dt is None:
            output("CFL coefficient was {0} and order correction 1/{1} = {2}"
                   .format(CFL,  (2 * approx_order + 1), order_corr))
        else:
            output("CFL coefficient {0} was ignored, dt specified directly"
                   .format(CFL))
        output("Courant number c = max(norm(a)) * dt/dx = {0}"
               .format(ts_conf.cour))
        output("Time stepping solver is {}".format(ts_conf.kind))
        output("... done.")

    return setup_cfl_condition


def get_1Dmesh_hook(XS, XE, n_nod):
    def mesh_hook(mesh, mode):
        """
        Generate the 1D mesh.
        """
        if mode == 'read':

            coors = nm.linspace(XS, XE, n_nod).reshape((n_nod, 1))
            conn = nm.arange(n_nod, dtype=nm.int32) \
                        .repeat(2)[1:-1].reshape((-1, 2))
            mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
            descs = ['1_2']

            mesh = Mesh.from_data('laplace_1d', coors, None,
                                  [conn], [mat_ids], descs)
            return mesh

        elif mode == 'write':
            pass

    return mesh_hook


def get_gen_block_mesh_hook(dims, shape, centre, mat_id=0, name='block',
                   coors=None, verbose=True):
    def mesh_hook(mesh, mode):
        """
        Generate the 1D mesh.
        """
        if mode == 'read':

            mesh = gen_block_mesh(dims, shape, centre, mat_id=mat_id, name=name,
                   coors=coors, verbose=verbose)
            return mesh

        elif mode == 'write':
            pass

    return mesh_hook


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
                        / nm.log( row["h"] / last_h)
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
        curr_df = err_df[err_df["order"] == o]
        co = ax.loglog(curr_df["n_cells"], curr_df["diff_l2"], 'o',
                       label=str(int(o)))[0].get_color()
        ax.loglog(curr_df["n_cells"], curr_df["diff_l2"], color=co, label="")
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


