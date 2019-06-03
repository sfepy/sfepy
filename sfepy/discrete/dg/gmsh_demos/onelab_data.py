import gmsh
import sys

if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " file [options]"
    exit(0)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.open(sys.argv[1])

# attempts to run a client selected when opening the file (e.g. a .pro file)
gmsh.onelab.run();

json = gmsh.onelab.get()
print json

gmsh.finalize()
