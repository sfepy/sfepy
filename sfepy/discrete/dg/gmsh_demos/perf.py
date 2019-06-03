import gmsh
import sys
gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("square")
gmsh.model.occ.addRectangle(0, 0, 0, 1, 1, 100)
gmsh.model.occ.synchronize()
gmsh.model.mesh.setTransfiniteSurface(100)
gmsh.model.mesh.generate(2)
gmsh.plugin.run("NewView")
print "before get"
type, tags, data, time, numComp = gmsh.view.getModelData(0, 0)
print "after get"
gmsh.finalize()
