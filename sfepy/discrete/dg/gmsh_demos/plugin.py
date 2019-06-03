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

# create a view with some data
t = gmsh.view.add("some data")
gmsh.view.addModelData(t, 0, "test", "NodeData",
                       [1, 2, 3, 4],
                       [[1.],[10.],[20.],[1.]])

# test getting data back
dataType, tags, data, time, numComp = gmsh.view.getModelData(t, 0)
print dataType, tags

# compute the iso-curve at value 11
gmsh.plugin.setNumber("Isosurface", "Value", 11.)
gmsh.plugin.run("Isosurface")

# delete the source view
gmsh.view.remove(t)

# check how many views the plugin created (a priori, a single list-based one)
tags = gmsh.view.getTags()
if len(tags) == 1:
    gmsh.view.write(tags[0], "iso.msh")
    # test getting data back
    dataTypes, numElements, data = gmsh.view.getListData(tags[0])
    print dataTypes, numElements

gmsh.finalize()
