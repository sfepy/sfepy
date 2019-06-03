# This file reimplements gmsh/tutorial/t5.geo in Python.

import gmsh
import math

model = gmsh.model
factory = model.geo

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

model.add("t5")

lcar1 = .1
lcar2 = .0005
lcar3 = .055

factory.addPoint(0.5,0.5,0.5, lcar2, 1)
factory.addPoint(0.5,0.5,0, lcar1, 2)
factory.addPoint(0,0.5,0.5, lcar1, 3)
factory.addPoint(0,0,0.5, lcar1, 4)
factory.addPoint(0.5,0,0.5, lcar1, 5)
factory.addPoint(0.5,0,0, lcar1, 6)
factory.addPoint(0,0.5,0, lcar1, 7)
factory.addPoint(0,1,0, lcar1, 8)
factory.addPoint(1,1,0, lcar1, 9)
factory.addPoint(0,0,1, lcar1, 10)
factory.addPoint(0,1,1, lcar1, 11)
factory.addPoint(1,1,1, lcar1, 12)
factory.addPoint(1,0,1, lcar1, 13)
factory.addPoint(1,0,0, lcar1, 14)

factory.addLine(8,9, 1);   factory.addLine(9,12, 2)
factory.addLine(12,11, 3); factory.addLine(11,8, 4)
factory.addLine(9,14, 5);  factory.addLine(14,13, 6)
factory.addLine(13,12, 7); factory.addLine(11,10, 8)
factory.addLine(10,13, 9); factory.addLine(10,4, 10)
factory.addLine(4,5, 11);  factory.addLine(5,6, 12)
factory.addLine(6,2, 13);  factory.addLine(2,1, 14)
factory.addLine(1,3, 15);  factory.addLine(3,7, 16)
factory.addLine(7,2, 17);  factory.addLine(3,4, 18)
factory.addLine(5,1, 19);  factory.addLine(7,8, 20)
factory.addLine(6,14, 21)

factory.addCurveLoop([-11,-19,-15,-18], 22)
factory.addPlaneSurface([22], 23)
factory.addCurveLoop([16,17,14,15], 24)
factory.addPlaneSurface([24], 25)
factory.addCurveLoop([-17,20,1,5,-21,13], 26)
factory.addPlaneSurface([26], 27)
factory.addCurveLoop([-4,-1,-2,-3], 28)
factory.addPlaneSurface([28], 29)
factory.addCurveLoop([-7,2,-5,-6], 30)
factory.addPlaneSurface([30], 31)
factory.addCurveLoop([6,-9,10,11,12,21], 32)
factory.addPlaneSurface([32], 33)
factory.addCurveLoop([7,3,8,9], 34)
factory.addPlaneSurface([34], 35)
factory.addCurveLoop([-10,18,-16,-20,4,-8], 36)
factory.addPlaneSurface([36], 37)
factory.addCurveLoop([-14,-13,-12,19], 38)
factory.addPlaneSurface([38], 39)

shells = []

# When the tag is not specified, a new one is automatically provided
sl = factory.addSurfaceLoop([35,31,29,37,33,23,39,25,27])
shells.append(sl)

def cheeseHole(x, y, z, r, lc, shells):
    p1 = factory.addPoint(x,  y,  z,   lc)
    p2 = factory.addPoint(x+r,y,  z,   lc)
    p3 = factory.addPoint(x,  y+r,z,   lc)
    p4 = factory.addPoint(x,  y,  z+r, lc)
    p5 = factory.addPoint(x-r,y,  z,   lc)
    p6 = factory.addPoint(x,  y-r,z,   lc)
    p7 = factory.addPoint(x,  y,  z-r, lc)

    c1 = factory.addCircleArc(p2,p1,p7)
    c2 = factory.addCircleArc(p7,p1,p5)
    c3 = factory.addCircleArc(p5,p1,p4)
    c4 = factory.addCircleArc(p4,p1,p2)
    c5 = factory.addCircleArc(p2,p1,p3)
    c6 = factory.addCircleArc(p3,p1,p5)
    c7 = factory.addCircleArc(p5,p1,p6)
    c8 = factory.addCircleArc(p6,p1,p2)
    c9 = factory.addCircleArc(p7,p1,p3)
    c10 = factory.addCircleArc(p3,p1,p4)
    c11 = factory.addCircleArc(p4,p1,p6)
    c12 = factory.addCircleArc(p6,p1,p7)

    l1 = factory.addCurveLoop([c5,c10,c4])
    l2 = factory.addCurveLoop([c9,-c5,c1])
    l3 = factory.addCurveLoop([c12,-c8,-c1])
    l4 = factory.addCurveLoop([c8,-c4,c11])
    l5 = factory.addCurveLoop([-c10,c6,c3])
    l6 = factory.addCurveLoop([-c11,-c3,c7])
    l7 = factory.addCurveLoop([-c2,-c7,-c12])
    l8 = factory.addCurveLoop([-c6,-c9,c2])

    s1 = factory.addSurfaceFilling([l1])
    s2 = factory.addSurfaceFilling([l2])
    s3 = factory.addSurfaceFilling([l3])
    s4 = factory.addSurfaceFilling([l4])
    s5 = factory.addSurfaceFilling([l5])
    s6 = factory.addSurfaceFilling([l6])
    s7 = factory.addSurfaceFilling([l7])
    s8 = factory.addSurfaceFilling([l8])

    sl = factory.addSurfaceLoop([s1, s2, s3, s4, s5, s6, s7, s8])
    v = factory.addVolume([sl])
    shells.append(sl)
    return v

x = 0; y = 0.75; z = 0; r = 0.09
for t in range(1, 6):
    x += 0.166
    z += 0.166
    v = cheeseHole(x, y, z, r, lcar3, shells)
    model.addPhysicalGroup(3, [v], t)

factory.addVolume(shells, 186)

model.addPhysicalGroup(3, [186], 10)
factory.synchronize()

model.mesh.generate(3)
gmsh.write("t5.msh")

gmsh.finalize()
