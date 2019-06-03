import gmsh
import sys
import math

gmsh.initialize(sys.argv)
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("2d surface mesh with purely discrete boundary");

# create a discrete curve with N nodes and N line elements
gmsh.model.addDiscreteEntity(1, 100)
N = 50
dt = 2*math.pi/N
pts = [[math.cos(i * dt), math.sin(i * dt), 0] for i in range(N)]
flat_pts = [item for sublist in pts for item in sublist]
gmsh.model.mesh.setNodes(1, 100, range(1, N+1), flat_pts)
n = [item for sublist in [[i, i+1] for i in range(1, N+1)] for item in sublist]
n[-1] = 1
gmsh.model.mesh.setElements(1, 100, [1], [range(1, N+1)], [n])

# create a plane surface from the discrete curve
gmsh.model.geo.addCurveLoop([100], 101)
gmsh.model.geo.addPlaneSurface([101], 102)
gmsh.model.geo.synchronize()

# mesh the surface
gmsh.model.mesh.generate(2)
gmsh.fltk.run()

gmsh.finalize()
