import gmsh
import sys

gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 1)

# Copied from discrete.py...
gmsh.model.add("test");
gmsh.model.addDiscreteEntity(2, 1)
gmsh.model.mesh.setNodes(2, 1, [1, 2, 3, 4],
                         [0., 0., 0.,
                          1., 0., 0.,
                          1., 1., 0.,
                          0., 1., 0.])
gmsh.model.mesh.setElements(2, 1, [2], [[1, 2]],
                            [[1, 2, 3,
                              1, 3, 4]])
# ... end of copy

# Create a new post-processing view
t = gmsh.view.add("some data")

# add 10 steps of model-based data, on the nodes of the mesh
for step in range(0, 10):
    gmsh.view.addModelData(t, step, "test", "NodeData",
                           [1, 2, 3, 4], # tags of nodes
                           [[10.],[10.],[12.+step],[13.+step]]) # data, per node
    
gmsh.view.write(t, "data.msh")
    
gmsh.finalize()
