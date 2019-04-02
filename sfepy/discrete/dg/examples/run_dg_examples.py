import numpy as nm

from os.path import join as pjoin


from sfepy.applications.pde_solver_app import PDESolverApp
from sfepy.base.conf import ProblemConf
from sfepy.base.ioutils import ensure_path
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)


from sfepy.discrete.dg.my_utils.read_plot_1Ddata import load_and_plot_fun
from sfepy.discrete.dg.my_utils.read_plot_1Ddata import clear_output_folder

def get_cfl_setup(CFL):

    def setup_cfl_condition(problem):
        """
        Sets up CFL condition for problem ts_conf in problem
        :param problem: discrete.problem.Problem
        :return:
        """
        ts_conf = problem.ts_conf
        mesh = problem.domain.mesh
        dim = mesh.dim
        approx_order = problem.fields['density'].approx_order
        mats = problem.create_materials('a')
        velo = problem.conf_materials['material_a__0'].values["val"]
        max_velo = nm.max(nm.linalg.norm(velo))
        dx = nm.min(problem.domain.mesh.cmesh.get_volumes(dim))
        order_corr = 1. / (2 * approx_order + 1)
        dt = dx / max_velo * CFL * order_corr
        # time_steps_N = int((tf - t0) / dt) * 2
        tn = int(nm.ceil((ts_conf.t1 - ts_conf.t0) / dt))
        dtdx = dt / dx

        ts_conf += Struct(dt=dt, n_step=tn)
        output("Preprocessing hook setup_cfl_condition:")
        output("Approximation order of field {} order is {}".format(problem.fields['density'].family_name, approx_order))
        output("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
        output("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
        output("CFL coefficient was {0} and order correction 1/{1} = {2}".format(CFL,  (2 * approx_order + 1), order_corr))
        output("Courant number c = max(norm(a)) * dt/dx = {0}".format(max_velo * dtdx))
        output("------------------------------------------")
        output("Time stepping solver is {}".format(ts_conf.kind))



    return setup_cfl_condition

if __name__ == '__main__':
    # TODO improve to accept options, get problem conf file from args
    pc = ProblemConf.from_file('example_dg_advection2D_tens.py')

    output_folder = "output"
    output_name_trunk_folder = pjoin(output_folder, pc.example_name, str(pc.approx_order)+"/")
    output_name_trunk_name = pc.example_name + str(pc.approx_order)
    output_name_trunk = pjoin(output_name_trunk_folder, output_name_trunk_name)
    ensure_path(output_name_trunk_folder)
    clear_output_folder(output_name_trunk_folder)

    sa = PDESolverApp(pc, Struct(output_filename_trunk=output_name_trunk,
                                 save_ebc=False,
                                 save_ebc_nodes=False,
                                 save_region=False,
                                 save_regions=False,
                                 save_regions_as_groups=False,
                                 save_field_meshes=False,
                                 solve_not=False), "sfepy")
    sa()


    if pc.dim == 1:
        load_times = min(pc.options.save_times, sa.problem.ts.n_step)

        load_and_plot_fun(output_name_trunk_folder, output_name_trunk_name,
                          pc.t0, pc.t1, load_times, pc.approx_order, pc.get_ic)