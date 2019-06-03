import gmsh
import sys

gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("test");

# add discrete surface with tag 1
gmsh.model.addDiscreteEntity(2, 1)

# add 4 mesh nodes
gmsh.model.mesh.setNodes(2, 1,
                         [1, 2, 3, 4], # node tags: 1, 2, 3, and 4
                         [0., 0., 0., # coordinates of node 1
                          1., 0., 0., # coordinates of node 2
                          1., 1., 0., # ...
                          0., 1., 0.])

# add 2 triangles
gmsh.model.mesh.setElements(2, 1,
                            [2], # single type : 3-node triangle
                            [[1, 2]], # triangle tags: 1 and 2
                            [[1, 2, 3, # triangle 1: nodes 1, 2, 3
                              1, 3, 4]]) # triangle 2: nodes 1, 3, 4

# export the mesh ; use explore.py to read and examine the mesh
gmsh.write("test.msh")

gmsh.finalize()

