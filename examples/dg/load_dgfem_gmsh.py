from sfepy.discrete.fem.meshio import Msh2MeshIO

gmsh_loader = Msh2MeshIO("output/ikucera1_tri/2/domain_2D.0000.msh")

mesh, data, times, times_n = gmsh_loader.read_data(setp=2)

pass