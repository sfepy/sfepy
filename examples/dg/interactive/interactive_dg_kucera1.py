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
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions, PeriodicBC, DGPeriodicBC, DGEssentialBC
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.mesh.mesh_generators import gen_block_mesh
from sfepy.discrete.fem.meshio import VTKMeshIO
from sfepy.base.ioutils import ensure_path

# sfepy terms
from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm
from sfepy.terms.terms_volume import LinearVolumeForceTerm
from terms.terms_diffusion import LaplaceTerm

from sfepy.discrete.variables import DGFieldVariable


# local imports
from sfepy.discrete.dg.dg_terms import AdvectDGFluxTerm, DiffusionDGFluxTerm, DiffusionInteriorPenaltyTerm, \
                                        NonlinScalarDotGradTerm, NonlinearHyperDGFluxTerm
from sfepy.discrete.dg.dg_tssolver import EulerStepSolver, TVDRK3StepSolver
from sfepy.discrete.dg.dg_field import DGField
from sfepy.discrete.dg.dg_limiters import IdentityLimiter, MomentLimiter1D


from sfepy.discrete.dg.my_utils.inits_consts \
    import left_par_q, gsmooth, const_u, ghump, superic

from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder

# ┌---------------┐
# |   Parameters  |
# |vvvvvvvvvvvvvvv|
approx_order = 2
diffusion_coef = 0.002
Cw = 15
CFL = .8
t0 = 0
t1 = 1
# |^^^^^^^^^^^^^^^|
# └---------------┘

# Setup  names
domain_name = "domain_2D"
problem_name = "ikucera1_tri"
output_folder = pjoin("output", problem_name, str(approx_order))
output_format = "msh"
mesh_output_folder = "output/mesh"
save_timestn = 100
clear_folder(pjoin(output_folder, "*." + output_format))

# ┌----------┐
# | Get mesh |
# └----------┘
mesh = gen_block_mesh((2., 2.), (21, 21), (0, 0))

mesh_name = "square_tri2"
mesh = Mesh.from_file("mesh/" + mesh_name + ".mesh")

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
                            'vertices in x == -1',
                            'edge')

right = domain.create_region('right',
                             'vertices in x == 1',
                             'edge')

top = domain.create_region('top',
                           'vertices in y == 1',
                           'edge')

bottom = domain.create_region('bottom',
                              'vertices in y == -1',
                              'edge')

field = DGField('dgfu', nm.float64, 'scalar', omega,
                approx_order=approx_order)

u = DGFieldVariable('u', 'unknown', field, history=1)
v = DGFieldVariable('v', 'test', field, primary_var_name='u')

MassT = DotProductVolumeTerm("adv_vol(v, u)", "v, u", integral, omega, u=u, v=v)

# ┌----------------------------┐
# |   Define advection terms   |
# └----------------------------┘
# a = Material('a', val=[velo])
# StiffT = ScalarDotMGradScalarTerm("adv_stiff(a.val, u, v)", "a.val, u[-1], v", integral, omega,
#                                   u=u, v=v, a=a)
#
# alpha = Material('alpha', val=[.0])
# FluxT = AdvectDGFluxTerm("adv_lf_flux(a.val, v, u)", "a.val, v,  u[-1]",
#                          integral, omega, u=u, v=v, a=a, alpha=alpha)


# ┌----------------------------┐
# |    Define Burgess terms    |
# └----------------------------┘
burg_velo = velo.T / nm.linalg.norm(velo)


def burg_fun(u):
    vu = 1/2*burg_velo * u[..., None] ** 2
    return vu


def burg_fun_d(u):
    v1 = burg_velo * u[..., None]
    return v1


nonlin = Material('nonlin', values={'.fun': burg_fun, '.dfun': burg_fun_d})

BurgStiffT = NonlinScalarDotGradTerm("burgess_stiff(f, df, u, v)", "nonlin.fun , nonlin.dfun, u[-1], v",
                                     integral, omega, u=u, v=v, nonlin=nonlin)

BurgFluxT = NonlinearHyperDGFluxTerm("burgess_lf_flux(f, df, u, v)", "nonlin.fun , nonlin.dfun, v, u[-1]",
                                     integral, omega, u=u, v=v, nonlin=nonlin)


# ┌----------------------------┐
# |   Define Diffusion terms   |
# └----------------------------┘
D = Material('D', val=[diffusion_coef])
Cwmat = Material("Cw", values={".val": Cw})

DivGrad = LaplaceTerm("diff_lap(D.val, v, u)", "D.val, v, u[-1]",
                      integral, omega, u=u, v=v, D=D)

DiffFluxT = DiffusionDGFluxTerm("diff_lf_flux(D.val, v, u)", "D.val, v,  u[-1]",
                                integral, omega, u=u, v=v, D=D)

DiffPen = DiffusionInteriorPenaltyTerm("diff_pen(Cw.val, v, u)", "Cw.val, v, u[-1]",
                                       integral, omega, u=u, v=v, Cw=Cwmat)

# ┌----------------------------┐
# |    Define source term      |
# └----------------------------┘
@Functionize
def source_fun(ts, coors, mode="qp", **kwargs):
    if mode == "qp":
        t = ts.dt * ts.step
        x_1 = coors[..., 0]
        x_2 = coors[..., 1]
        sin = nm.sin
        cos = nm.cos
        exp = nm.exp
        res = (
                + (5 * x_1 * cos(5 * x_1 * x_2) - 4 * (x_1 - 1) * cos(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1) ** 2 * (sin(5 * x_1 * x_2) - sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2))
                + (5 * x_2 * cos(5 * x_1 * x_2) - 4 * (x_2 - 1) * cos(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1) ** 2 * (sin(5 * x_1 * x_2) - sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2))
                - ((25 * x_1 ** 2 * sin(5 * x_1 * x_2) - 16 * (x_1 - 1) ** 2 * sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1)
                 + (25 * x_2 ** 2 * sin(5 * x_1 * x_2) - 16 * (x_2 - 1) ** 2 * sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * (exp(-t) - 1)) * diffusion_coef
                + (sin(5 * x_1 * x_2) - sin(4 * x_1 * x_2 - 4 * x_1 - 4 * x_2)) * exp(-t)
        )
        return {"val": res[..., None, None]}


g = Material('g', function=source_fun)
SourceTerm = LinearVolumeForceTerm("source_g(g.val, v)", "g.val, v", integral, omega, v=v, g=g)


# ┌----------------------------┐
# |===== Create equation ======|
# └----------------------------┘
eq = Equation('balance',
              MassT
              + BurgStiffT - BurgFluxT
              - (+ DivGrad - DiffFluxT) - diffusion_coef * DiffPen
              + SourceTerm
              )
eqs = Equations([eq])

# ┌----------------------------┐
# | Create boundary conditions |
# └----------------------------┘
@Functionize
def bc_funs(ts, coors, bc, problem):
    # return 2*coors[..., 1]
    t = ts.dt*ts.step
    x_1 = coors[..., 0]
    x_2 = coors[..., 1]
    sin = nm.sin
    cos = nm.cos
    exp = nm.exp
    if bc.diff == 0:
        if "left" in bc.name:
            res = -(exp(-t) - 1)*(sin(-5*x_2) + sin(8*x_2 - 4))
        elif "bottom" in bc.name:
            res = -(exp(-t) - 1) * (sin(-5 * x_1) + sin(8 * x_1 - 4))
        elif "right" in bc.name:
            res = -(exp(-t) - 1)*(sin(4) + sin(5*x_2))
        elif "top" in bc.name:
            res = -(exp(-t) - 1)*(sin(4) + sin(5*x_1))

    elif bc.diff == 1:
        if "left" in bc.name:
            res = nm.stack(((4*(x_2 - 1)*cos(4) - 5*x_2*cos(5*x_2))*(exp(-t) - 1),
                                -5*(exp(-t) - 1)*cos(5*x_2)),
                           axis=-2)
        elif "bottom" in bc.name:
            res = nm.stack(((5*cos(-5*x_1) - 8*cos(8*x_1 - 4))*(exp(-t) - 1),
                            -(5*x_1*cos(-5*x_1) - 4*(x_1 - 1)*cos(8*x_1 - 4))*(exp(-t) - 1)),
                           axis=-2)

        elif "right" in bc.name:
            res = nm.stack(((4*(x_2 - 1)*cos(4) - 5*x_2*cos(5*x_2))*(exp(-t) - 1),
                            -5*(exp(-t) - 1)*cos(5*x_2)),
                           axis=-2)
        elif "top" in bc.name:
            res = nm.stack((-5*(exp(-t) - 1)*cos(5*x_1),
                            (4*(x_1 - 1)*cos(4) - 5*x_1*cos(5*x_1))*(exp(-t) - 1)),
                           axis=-2)

    return res


left_bc_u = DGEssentialBC('left_fix_u', left, {'u.all': bc_funs})
right_bc_u = DGEssentialBC('right_fix_u', right, {'u.all': bc_funs})
top_bc_u = DGEssentialBC('top_fix_u', top, {'u.all': bc_funs})
bottom_bc_u = DGEssentialBC('bottom_fix_u', bottom, {'u.all': bc_funs})

left_bc_du = DGEssentialBC('left_fix_du', left, {'u.all': bc_funs}, diff=1)
right_bc_du = DGEssentialBC('right_fix_du', right, {'u.all': bc_funs}, diff=1)
top_bc_du = DGEssentialBC('top_fix_du', top, {'u.all': bc_funs}, diff=1)
bottom_bc_du = DGEssentialBC('bottom_fix_du', bottom, {'u.all': bc_funs}, diff=1)


# ┌--------------------------┐
# | Create initial condition |
# └--------------------------┘
@Functionize
def ic_fun(x, ic=None):
    return 0*gsmooth(x[..., 0:1] - .3) * gsmooth(x[..., 1:] - .3)


ics = InitialCondition('ic', omega, {'u.0': ic_fun})

# ┌----------------┐
# | Create problem |
# └----------------┘
pb = Problem(problem_name, equations=eqs, conf=Struct(options={"save_times": save_timestn}, ics={},
                                                      ebcs={}, epbcs={}, lcbcs={}, materials={},
                                                      ),
             active_only=False)

pb.set_ics(Conditions([ics]))
pb.set_bcs(
        ebcs=Conditions([left_bc_u, bottom_bc_u,
                          left_bc_du, bottom_bc_du,
                          right_bc_u, top_bc_u,
                          right_bc_du, top_bc_du
                          ]))

pb.setup_output(output_dir=output_folder, output_format=output_format)

# ┌----------------┐
# | Choose limiter |
# └----------------┘
limiter = IdentityLimiter

# ┌-------------------------┐
# | Set time discretization |
# └-------------------------┘
max_velo = nm.max(nm.linalg.norm(velo))
dx = nm.min(mesh.cmesh.get_volumes(2))
# dt = dx / max_velo * CFL / (2 * approx_order + 1)
dt = 1e-5
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
print("and EBCs: {}".format(pb.ebcs.names))
print("and EPBCS: {}".format(pb.epbcs.names))
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

print("Equation \n\n\t\t u_t - div(au)) = 0\n")
print("With IC: {}".format(ic_fun.name))
print("and EBCs: {}".format(pb.ebcs.names))
print("and EPBCS: {}".format(pb.epbcs.names))
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
