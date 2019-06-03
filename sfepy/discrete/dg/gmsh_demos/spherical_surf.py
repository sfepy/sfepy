import gmsh
import sys
import math

gmsh.initialize(sys.argv)
gmsh.model.add("sphere_cut")
R = 1
R1 = 0.95
sph = gmsh.model.occ.addSphere(0,0,0, R, -1, 0, math.pi/2, math.pi/2)
b1 = gmsh.model.occ.addBox(R1,0,0, R,R,R)
b2 = gmsh.model.occ.addBox(0,R1,0, R,R,R)
b3 = gmsh.model.occ.addBox(0,0,R1, R,R,R)
gmsh.model.occ.cut([(3,sph)], [(3,b1), (3,b2), (3,b3)])
gmsh.model.occ.synchronize()

gmsh.model.removeEntities([(3,sph)])
gmsh.model.removeEntities([(2,2), (2,4), (2,6)], True)
gmsh.fltk.run()
gmsh.finalize()
