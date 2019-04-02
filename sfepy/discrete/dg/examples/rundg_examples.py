from os.path import join as pjoin


from sfepy.applications.pde_solver_app import PDESolverApp
from sfepy.base.conf import ProblemConf
from sfepy.base.base import Struct
from sfepy.base.ioutils import ensure_path


from sfepy.discrete.dg.my_utils.read_plot_1Ddata import load_and_plot_fun
from sfepy.discrete.dg.my_utils.read_plot_1Ddata import clear_output_folder

def setup_cfl_condition():
    # TODO
    pass

if __name__ == '__main__':
    # TODO improve to accept options, get problem conf file from args
    pc = ProblemConf.from_file("dict_dg_adv1.py")

    output_folder = "output"
    output_name_trunk_folder = pjoin(output_folder, pc.example_name, str(pc.approx_order)+"/")
    output_name_trunk_name = pc.example_name + str(pc.approx_order)
    output_name_trunk = pjoin(output_name_trunk_folder, output_name_trunk_name)
    ensure_path(output_name_trunk_folder)
    clear_output_folder(output_name_trunk_folder)

    # tODO print everything using sfepy output, maybe even log?
    # print("Solving equation \n\n\t\t u_t - div(au)) = 0\n")
    # print("With IC: {}".format(ic_fun.name))
    # # print("and EBCs: {}".format(pb.ebcs.names))
    # # print("and EPBCS: {}".format(pb.epbcs.names))
    # print("-------------------------------------")
    # print("Approximation order is {}".format(pc.approx_order))
    # print("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
    # print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
    # print("CFL coefficient was {0} and order correction {1}".format(CFL, 1 / (2 * approx_order + 1)))
    # print("Courant number c = max(abs(u)) * dt/dx = {0}".format(max_velo * dtdx))
    # print("------------------------------------------")
    # print("Time stepping solver is {}".format(tss.name))
    # print("Limiter used: {}".format(limiter.name))
    # print("======================================")

    sa = PDESolverApp(pc, Struct(output_filename_trunk=output_name_trunk,
                                 save_ebc=False,
                                 save_ebc_nodes=False,
                                 save_region=False,
                                 save_regions=False,
                                 save_regions_as_groups=False,
                                 save_field_meshes=False,
                                 solve_not=False) , "sfepy")
    sa()


    if pc.dim == 1:
        load_times = min(pc.options.save_times, pc.solvers["solvers_tss__0"].n_step)

        load_and_plot_fun(output_name_trunk_folder, output_name_trunk_name,
                          pc.t0, pc.t1, load_times, pc.approx_order, pc.get_ic)