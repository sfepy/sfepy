import numpy as nm

from sfepy.discrete.equations import Equation, Equations
from sfepy.discrete.variables import FieldVariable
from sfepy.discrete.fem import FEDomain
from sfepy.discrete.fem.meshio import Msh2MeshIO
from sfepy.discrete.fem.mesh import Mesh
from sfepy.discrete.functions import make_sfepy_function, Function
from sfepy.discrete.integrals import Integral, Integrals
from sfepy.discrete.materials import Material
from sfepy.discrete.problem import Problem
from sfepy.terms.terms_basic import SurfaceTerm


from sfepy.discrete.dg.dg_field import DGField


example_name = "kucera1"
mesh_name = "mesh_simp_2D_11_750.vtk"
approx_order = 2


data_path = "output\\{example_name}\\{approx_order}\\{example_name}{approx_order}.0000.msh"\
    .format(example_name=example_name, approx_order=approx_order)
mesh_path = "mesh\\" + mesh_name


gmsh_loader = Msh2MeshIO(data_path)


mesh = Mesh.from_file(mesh_path)
dmesh = Mesh()
gmsh_loader.read(dmesh, drop_z=True)


data, times, times_n, scheme = gmsh_loader.read_data(step="all")

assert(approx_order == nm.max(scheme.P))

domain = FEDomain(example_name, dmesh)
omega = domain.create_region('Omega', 'all')

field = DGField('dgfu', nm.float64, 'scalar', omega,
                approx_order=approx_order)

def u_fun(coors, t):
    x_1 = coors[..., 0]
    x_2 = coors[..., 1]
    sin = nm.sin
    exp = nm.exp

    res = (sin(4 * (x_1 + x_2 - x_1 * x_2)) + sin(5 * x_1 * x_2)) * (1 - exp(-t))
    return res


@make_sfepy_function
def sol_fun(ts, coors, mode="qp", **kwargs):
    t = ts.time
    if mode == "qp":
        return {"u":  u_fun(coors, times[-1])[..., None, None]}


nodes, nodal_vals = field.get_nodal_values(data[-1], omega)

exact_nodal_vals = u_fun(nodes, times[-1])

# nod, val = field.set_dofs(lambda coors: solution_fun(coors, times[-1]), region=omega)
#
# nodes, exact_nodal_vals = field.get_nodal_values(val, omega)

err_abs = nm.abs(exact_nodal_vals - nodal_vals)

from sfepy.discrete.common.mappings import get_jacobian

# Sufficient quadrature order for the analytical expression.
idiff = Integral('idiff', max(approx_order, 10))

u = FieldVariable("u", "unknown", field)
u.set_data(data[-1])

eqs = Equations([Equation('balance', SurfaceTerm("s()", "u", idiff, omega, u=u))])
pb = Problem("err_est", equations=eqs)

num_qp = pb.evaluate('ev_volume_integrate.idiff.Omega(u)',
                      u=u,
                      integrals=Integrals([idiff]), mode='qp')
aux = Material('aux', function=sol_fun)

ana_qp = pb.evaluate('ev_volume_integrate_mat.idiff.Omega(aux.u, u)',
                      aux=aux, u=u,
                      integrals=Integrals([idiff]), mode='qp')

det = get_jacobian(field, idiff)

diff_l2 = nm.sqrt((((num_qp - ana_qp)**2) * det).sum())
ana_l2 = nm.sqrt((((ana_qp)**2) * det).sum())
error = diff_l2 / ana_l2

pass



