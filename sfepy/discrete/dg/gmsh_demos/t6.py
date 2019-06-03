# This file reimplements gmsh/tutorial/t6.geo in Python.

import gmsh
import math

model = gmsh.model
factory = model.geo

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

model.add("t6")

# Copied from t1.py...
lc = 1e-2
factory.addPoint(0, 0, 0, lc, 1)
factory.addPoint(.1, 0,  0, lc, 2)
factory.addPoint(.1, .3, 0, lc, 3)
factory.addPoint(0, .3, 0, lc, 4)
factory.addLine(1, 2, 1)
factory.addLine(3, 2, 2)
factory.addLine(3, 4, 3)
factory.addLine(4, 1, 4)
factory.addCurveLoop([4, 1, -2, 3], 1)
factory.addPlaneSurface([1], 1)
model.addPhysicalGroup(0, [1, 2], 1)
model.addPhysicalGroup(1, [1, 2], 2)
model.addPhysicalGroup(2, [1], 6)
model.setPhysicalName(2, 6, "My surface")
# ...end of copy

# Delete surface 1 and left boundary (curve 4)
factory.remove([[2,1], [1,4]])

# Replace left boundary with 3 new lines
p1 = factory.addPoint(-0.05, 0.05, 0, lc)
p2 = factory.addPoint(-0.05, 0.1, 0, lc)
l1 = factory.addLine(1, p1)
l2 = factory.addLine(p1, p2)
l3 = factory.addLine(p2, 4)

# Recreate surface
factory.addCurveLoop([2, -1, l1, l2, l3, -3], 2)
factory.addPlaneSurface([-2], 1)

# Put 20 points with a refinement toward the extremities on curve 2
factory.mesh.setTransfiniteCurve(2, 20, "Bump", 0.05)

# Put 20 points total on combination of curves l1, l2 and l3 (beware that the
# points p1 and p2 are shared by the curves, so we do not create 6 + 6 + 10 = 22
# points, but 20!)
factory.mesh.setTransfiniteCurve(l1, 6)
factory.mesh.setTransfiniteCurve(l2, 6)
factory.mesh.setTransfiniteCurve(l3, 10)

# Put 30 points following a geometric progression on curve 1 (reversed) and on
# curve 3
factory.mesh.setTransfiniteCurve(1, 30, "Progression", -1.2)
factory.mesh.setTransfiniteCurve(3, 30, "Progression", 1.2)

# Define the Surface as transfinite, by specifying the four corners of the
# transfinite interpolation
factory.mesh.setTransfiniteSurface(1, "Left", [1,2,3,4])

# Recombine the triangles into quads
factory.mesh.setRecombine(2, 1)

# Apply an elliptic smoother to the grid
gmsh.option.setNumber("Mesh.Smoothing", 100)
model.addPhysicalGroup(2, [1], 1)

# When the surface has only 3 or 4 control points, the transfinite constraint
# can be applied automatically (without specifying the corners explictly).
factory.addPoint(0.2, 0.2, 0, 1.0, 7)
factory.addPoint(0.2, 0.1, 0, 1.0, 8)
factory.addPoint(0, 0.3, 0, 1.0, 9)
factory.addPoint(0.25, 0.2, 0, 1.0, 10)
factory.addPoint(0.3, 0.1, 0, 1.0, 11)
factory.addLine(8, 11, 10)
factory.addLine(11, 10, 11)
factory.addLine(10, 7, 12)
factory.addLine(7, 8, 13)
factory.addCurveLoop([13, 10, 11, 12], 14)
factory.addPlaneSurface([14], 15)
for i in range(10,14):
    factory.mesh.setTransfiniteCurve(i, 10)
factory.mesh.setTransfiniteSurface(15)
model.addPhysicalGroup(2, [15], 2)

model.mesh.generate(2)
gmsh.write("t6.msh")
gmsh.finalize()
