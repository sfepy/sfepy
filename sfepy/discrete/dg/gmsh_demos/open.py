import gmsh
import sys

if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " file.geo [options]"
    exit(0)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.open(sys.argv[1])
gmsh.model.mesh.generate(3)
gmsh.write("test.msh")
gmsh.finalize()

