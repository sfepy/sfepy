import numpy as nm
import matplotlib.pyplot as plt

# sfepy imports
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.base.base import Struct
from sfepy.base.base import IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete.conditions import InitialCondition, EssentialBC, Conditions
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton
from sfepy.solvers.ts_solvers import SimpleTimeSteppingSolver

from sfepy.terms.terms_dot import ScalarDotMGradScalarTerm, DotProductVolumeTerm


from sfepy.base.conf import ProblemConf


# local imports
from dg_terms import AdvFluxDGTerm, ScalarDotMGradScalarDGTerm
# # from dg_equation import Equation
from dg_tssolver import EulerStepSolver, TVDRK3StepSolver
from dg_field import DGField

from my_utils.inits_consts import left_par_q, gsmooth, const_u, ghump, superic
from my_utils.visualizer import load_1D_vtks, plot1D_DG_sol

X1 = 0.
XN = 2*nm.pi
n_nod = 100
n_el = n_nod - 1
coors = nm.linspace(X1, XN, n_nod).reshape((n_nod, 1))
conn = nm.arange(n_nod, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
descs = ['1_2']
mesh = Mesh.from_data('advection_book_1d', coors, None,
                      [conn], [mat_ids], descs)

velo = 2*nm.pi
max_velo = nm.max(nm.abs(velo))

t0 = 0
t1 = 0.8
dx = (XN - X1) / n_nod
dt = dx / nm.abs(velo) *.4
# time_steps_N = int((tf - t0) / dt) * 2
tn = int(nm.ceil((t1 - t0) / dt))
dtdx = dt / dx
print("Space divided into {0} cells, {1} steps, step size is {2}".format(mesh.n_el, len(mesh.coors), dx))
print("Time divided into {0} nodes, {1} steps, step size is {2}".format(tn - 1, tn, dt))
print("Courant number c = max(abs(u)) * dt/dx = {0}".format(max_velo * dtdx))

approx_order = 0

integral = Integral('i', order=7)
domain = FEDomain('adv_sin_1D', mesh)
omega = domain.create_region('Omega', 'all')
left = domain.create_region('Gamma1',
                              'vertices in x == %.10f' % X1,
                              'vertex')
right = domain.create_region('Gamma2',
                              'vertices in x > %.10f' % (XN - dx),
                              'vertex')
field = DGField('dgfu', nm.float64, 'scalar', omega,
                approx_order=approx_order)

u = FieldVariable('u', 'unknown', field, history=1)
v = FieldVariable('v', 'test', field, primary_var_name='u')


MassT = DotProductVolumeTerm("adv_vol(v, u)", "v, u", integral, omega, u=u, v=v)

a = Material('a', val=[velo])
StiffT = ScalarDotMGradScalarDGTerm("adv_stiff(a.val, v, u)", "a.val, u, v", integral, omega,
                                    u=u, v=v, a=a)

FluxT = AdvFluxDGTerm(integral, omega, u=u, v=v, a=a)

eq = Equation('balance', MassT - StiffT + FluxT)
eqs = Equations([eq])

left_fix_u = EssentialBC('left_fix_u', left, {'u.all' : 0.0})
#
right_fix_u = EssentialBC('right_fix_u', right, {'u.all' : 0.0})

def ic_wrap(x, ic=None):
    return nm.sin(x)


ic_fun = Function('ic_fun', ic_wrap)
ics = InitialCondition('ic', omega, {'u.0': ic_fun})

pb = Problem('advection', equations=eqs)
pb.setup_output(output_dir="./output/") #, output_format="msh")
pb.set_bcs(ebcs=Conditions([left_fix_u, right_fix_u]))
pb.set_ics(Conditions([ics]))

# create post stage hook with limiter
from dg_field import get_unraveler, get_raveler
from dg_limiters import moment_limiter_1D
def limiter(vec):
    # TODO unify shapes for limiter
    u = get_unraveler(field.n_el_nod, field.n_cell)(vec).swapaxes(0, 1)
    u = moment_limiter_1D(u)
    rvec = get_raveler(field.n_el_nod, field.n_cell)(u.swapaxes(0, 1))
    return rvec[:, 0]

ls = ScipyDirect({})
nls_status = IndexedStruct()
# nls = Newton({'is_linear' : True}, lin_solver=ls, status=nls_status)
# nls = EulerStepSolver({}, lin_solver=ls, status=nls_status)
nls = Newton({}, lin_solver=ls, status=nls_status, post_stage_hook=limiter)

tss = TVDRK3StepSolver({'t0' : t0, 't1' : t1, 'n_step': tn},
                                nls=nls, context=pb, verbose=True)
pb.set_solver(tss)
pb.solve()


#--------
#| Plot |
#--------
lmesh, u = load_1D_vtks("./output/", "adv_sin_1D", tn, order=approx_order)
plot1D_DG_sol(lmesh, t0, t1, u, dt=dt, ic=ic_wrap,
              delay=100, polar=False)