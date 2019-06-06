import numpy as nm


from sfepy.discrete.fem import FEDomain
from sfepy.discrete.fem.meshio import Msh2MeshIO
from sfepy.discrete.fem.mesh import Mesh


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

nodes, nodal_vals = field.get_nodal_values(data[-1], omega)


def solution_fun(coors, t):
    x_1 = coors[..., 0]
    x_2 = coors[..., 1]
    sin = nm.sin
    exp = nm.exp

    res = (sin(4*(x_1 + x_2 - x_1 * x_2)) + sin(5*x_1*x_2))*(1 - exp(-t))
    return res


exact_nodal_vals = solution_fun(nodes, times[-1])

# nod, val = field.set_dofs(lambda coors: solution_fun(coors, times[-1]), region=omega)
#
# nodes, exact_nodal_vals = field.get_nodal_values(val, omega)

err_abs = nm.abs(exact_nodal_vals - nodal_vals)

pass



