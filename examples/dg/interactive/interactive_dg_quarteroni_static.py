import numpy as nm
from numpy.linalg import norm
import matplotlib.pyplot as plt
from os.path import join as pjoin

# sfepy imports
from sfepy.discrete.fem.periodic import match_x_line, match_y_line
from sfepy.discrete.functions import Functionize
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct
from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions, PeriodicBC,  DGPeriodicBC, DGEssentialBC
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_solvers import SimpleTimeSteppingSolver
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.meshio import VTKMeshIO
from sfepy.base.ioutils import ensure_path


from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm
from terms.terms_diffusion import LaplaceTerm
from terms.terms_volume import LinearVolumeForceTerm



# local imports
from sfepy.discrete.variables import DGFieldVariable

from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm, DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm
from sfepy.discrete.dg.dg_tssolver import EulerStepSolver, TVDRK3StepSolver
from sfepy.discrete.dg.dg_field import DGField
from sfepy.discrete.dg.dg_limiters import IdentityLimiter


from sfepy.discrete.dg.my_utils.inits_consts \
    import left_par_q, gsmooth, const_u, ghump, superic

from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder

# ┌---------------┐
# |   Parameters  |
# |vvvvvvvvvvvvvvv|
approx_order = 3
CFL = .8
t0 = 0
t1 = 1
# |^^^^^^^^^^^^^^^|
# └---------------┘

# Setup  names
domain_name = "domain_2D"
problem_name = "iquarteroni2_static"
output_folder = pjoin("output", problem_name, str(approx_order))
output_format = "msh"
mesh_output_folder = "output/mesh"
save_timestn = 100
clear_folder(pjoin(output_folder, "*." + output_format))


# ┌----------┐
# | Get mesh |
# └----------┘
mesh = gen_block_mesh((1., 1.), (50, 50), (0.5, 0.5))

mesh_name = "tens_2D_mesh"
mesh = Mesh.from_file("mesh/" + mesh_name + ".vtk")

angle = - nm.pi / 5
rotm = nm.array([[nm.cos(angle), -nm.sin(angle)],
                 [nm.sin(angle), nm.cos(angle)]])
velo = -nm.sum(rotm.T * nm.array([1., 0.]), axis=-1)[:, None]
velo = nm.array([[1., 1.]]).T
max_velo = nm.max(nm.linalg.norm(velo))


# ┌---------------------------┐
# | Create problem components |
# └---------------------------┘
integral = Integral('i', order=approx_order * 2)
domain = FEDomain(domain_name, mesh)
omega = domain.create_region('Omega', 'all')

left = domain.create_region('left',
                            'vertices in x == 0',
                            'edge')

right = domain.create_region('right',
                             'vertices in x == 1',
                             'edge')

top = domain.create_region('top',
                           'vertices in y == 1',
                           'edge')

bottom = domain.create_region('bottom',
                              'vertices in y == 0',
                              'edge')

field = DGField('dgfu', nm.float64, 'scalar', omega,
                approx_order=approx_order)

u = DGFieldVariable('u', 'unknown', field, history=1)
v = DGFieldVariable('v', 'test', field, primary_var_name='u')


# ┌----------------------------┐
# |   Define advection terms   |
# └----------------------------┘
MassT = DotProductVolumeTerm("adv_vol(v, u)", "v, u", integral, omega, u=u, v=v)

a = Material('a', val=[velo])
StiffT = ScalarDotMGradScalarTerm("adv_stiff(a.val, u, v)", "a.val, u[-1], v", integral, omega,
                                  u=u, v=v, a=a)

alpha = Material('alpha', val=[.0])
FluxT = AdvectDGFluxTerm("adv_lf_flux(a.val, v, u)", "a.val, v,  u[-1]",
                         integral, omega, u=u, v=v, a=a, alpha=alpha)


# ┌----------------------------┐
# |   Define Diffusion terms   |
# └----------------------------┘
diffusion_tensor = 0.02  # nm.array([[.002, 0],
#           [0, .002]]).T
D = Material('D', val=[diffusion_tensor])
Cw = Material("Cw", values={".val": 10})

DivGrad = LaplaceTerm("diff_lap(D.val, v, u)", "D.val, v, u[-1]",
                      integral, omega, u=u, v=v, D=D)

DiffFluxT = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)", "D.val, v,  u[-1]",
                                integral, omega, u=u, v=v, D=D)
DiffPen = DiffusionInteriorPenaltyTerm("diff_pen(Cw.val, v, u)", "Cw.val, v, u[-1]",
                                       integral, omega, u=u, v=v, Cw=Cw)


# ┌----------------------------┐
# |    Define source term      |
# └----------------------------┘
@Functionize
def source_fun(ts, coors, mode="qp", **kwargs):
    # t = ts.dt * ts.step
    x_1 = coors[..., 0]
    x_2 = coors[..., 1]
    sin = nm.sin
    cos = nm.cos
    exp = nm.exp
    sqrt = nm.sqrt
    if mode == "qp":
        eps = diffusion_tensor
        res = (-1024 * eps * (8 * (4 * (2 * x_1 - 1) ** 2 + 
                                   4 * (2 * x_2 - 1) ** 2 - 1) * (2 * x_1 - 1) ** 2 /
                              (((4 * (2 * x_1 - 1) ** 2 + 4 * (2 * x_2 - 1) ** 2 - 1) ** 2 / eps + 256) ** 2 * eps ** (3 / 2))
                              + 8 * (4 * (2 * x_1 - 1) ** 2 + 
                                     4 * (2 * x_2 - 1) ** 2 - 1) * (2 * x_2 - 1) ** 2 /
                              (((4 * (2 * x_1 - 1) ** 2 + 4 * (2 * x_2 - 1) ** 2 - 1) ** 2 / eps + 256) ** 2 * eps ** (3 / 2))
                              - 1 / (((4 * (2 * x_1 - 1) ** 2 + 4 * (2 * x_2 - 1) ** 2 - 1) ** 2 / eps + 256) * sqrt(eps)))
               - 256 * (2 * x_1 - 1) / (((4 * (2 * x_1 - 1) ** 2 + 4 * (2 * x_2 - 1) ** 2 - 1) ** 2 / eps + 256) * sqrt(eps)) 
               - 256 * (2 * x_2 - 1) / (((4 * (2 * x_1 - 1) ** 2 + 4 * (2 * x_2 - 1) ** 2 - 1) ** 2 / eps + 256) * sqrt(eps))
               )
        return {"val": res[..., None, None]}


g = Material('g', function=source_fun)
SourceTerm = LinearVolumeForceTerm("source_g(g.val, v)", "g.val, v", integral, omega, v=v, g=g)

# ┌----------------------------┐
# |===== Create equation ======|
# └----------------------------┘
eq = Equation('balance',
              (DivGrad - StiffT)
              - DiffFluxT
              + (diffusion_tensor * DiffPen + FluxT)
              +  SourceTerm
              )
eqs = Equations([eq])

# ┌----------------------------┐
# | Create boundary conditions |
# └----------------------------┘
# right_fix_u = EssentialBC('right_fix_u', right, {'u.all' : 0.0})


# ┌--------------------------┐
# | Create initial condition |
# └--------------------------┘
@Functionize
def ic_fun(x, ic=None):
    return gsmooth(x[..., 0:1] - .3) * gsmooth(x[..., 1:] - .3)


ics = InitialCondition('ic', omega, {'u.0': ic_fun})

# ┌----------------┐
# | Create problem |
# └----------------┘
pb = Problem(problem_name, equations=eqs, conf=Struct(options={"save_times": save_timestn}, ics={},
                                                      ebcs={}, epbcs={}, lcbcs={}, materials={},
                                                      ),
             active_only=False)
pb.setup_output(output_dir=output_folder, output_format=output_format)
pb.set_ics(Conditions([ics]))
pb.set_bcs(ebcs=Conditions([
    # dirichlet_bc_u
]))

# ┌----------------┐
# | Choose limiter |
# └----------------┘
limiter = IdentityLimiter


max_velo = nm.max(nm.linalg.norm(velo))
dx = nm.min(mesh.cmesh.get_volumes(2))
dt = dx / max_velo * CFL / (2 * approx_order + 1)
tn = int(nm.ceil((t1 - t0) / dt))
dtdx = dt / dx

# ┌---------------┐
# | Create solver |
# └---------------┘
ls = ScipyDirect({})
nls_status = IndexedStruct()
nls = Newton({'is_linear': True}, lin_solver=ls, status=nls_status)

tss_conf = {'t0'     : t0,
            't1'     : t1,
            'n_step' : tn,
            "limiter": limiter}

tss = EulerStepSolver(tss_conf,
                      nls=nls, context=pb, verbose=True)

# ┌-------┐
# | Solve |
# └-------┘
print("Solving equation \n\n\t\t u_t - div(au)) = 0\n")
print("With IC: {}".format(ic_fun.name))
# print("and EBCs: {}".format(pb.ebcs.names))
# print("and EPBCS: {}".format(pb.epbcs.names))
print("-------------------------------------")
print("Approximation order is {}".format(approx_order))
print("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
print("CFL coefficient was {0} and order correction {1}".format(CFL, 1 / (2 * approx_order + 1)))
print("Courant number c = max(abs(u)) * dt/dx = {0}".format(max_velo * dtdx))
print("------------------------------------------")
print("Time stepping solver is {}".format(tss.name))
# print("Limiter used: {}".format(limiter.name))
print("======================================")

pb.set_solver(tss)
pb.solve()
