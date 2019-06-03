# This file reimplements gmsh/tutorial/t4.geo in Python.

import gmsh
import math

model = gmsh.model
factory = model.geo

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

model.add("t4")

cm = 1e-02
e1 = 4.5 * cm; e2 = 6 * cm / 2; e3 =  5 * cm / 2
h1 = 5 * cm; h2 = 10 * cm; h3 = 5 * cm; h4 = 2 * cm; h5 = 4.5 * cm
R1 = 1 * cm; R2 = 1.5 * cm; r = 1 * cm
Lc1 = 0.01
Lc2 = 0.003

def hypot(a, b):
    return math.sqrt(a * a + b * b)

ccos = (-h5*R1 + e2 * hypot(h5, hypot(e2, R1))) / (h5*h5 + e2*e2)
ssin = math.sqrt(1 - ccos*ccos)

factory.addPoint(-e1-e2, 0    , 0, Lc1, 1)
factory.addPoint(-e1-e2, h1   , 0, Lc1, 2)
factory.addPoint(-e3-r , h1   , 0, Lc2, 3)
factory.addPoint(-e3-r , h1+r , 0, Lc2, 4)
factory.addPoint(-e3   , h1+r , 0, Lc2, 5)
factory.addPoint(-e3   , h1+h2, 0, Lc1, 6)
factory.addPoint( e3   , h1+h2, 0, Lc1, 7)
factory.addPoint( e3   , h1+r , 0, Lc2, 8)
factory.addPoint( e3+r , h1+r , 0, Lc2, 9)
factory.addPoint( e3+r , h1   , 0, Lc2, 10)
factory.addPoint( e1+e2, h1   , 0, Lc1, 11)
factory.addPoint( e1+e2, 0    , 0, Lc1, 12)
factory.addPoint( e2   , 0    , 0, Lc1, 13)

factory.addPoint( R1 / ssin, h5+R1*ccos, 0, Lc2, 14)
factory.addPoint( 0        , h5        , 0, Lc2, 15)
factory.addPoint(-R1 / ssin, h5+R1*ccos, 0, Lc2, 16)
factory.addPoint(-e2       , 0.0       , 0, Lc1, 17)

factory.addPoint(-R2 , h1+h3   , 0, Lc2, 18)
factory.addPoint(-R2 , h1+h3+h4, 0, Lc2, 19)
factory.addPoint( 0  , h1+h3+h4, 0, Lc2, 20)
factory.addPoint( R2 , h1+h3+h4, 0, Lc2, 21)
factory.addPoint( R2 , h1+h3   , 0, Lc2, 22)
factory.addPoint( 0  , h1+h3   , 0, Lc2, 23)
                                                
factory.addPoint( 0, h1+h3+h4+R2, 0, Lc2, 24)
factory.addPoint( 0, h1+h3-R2,    0, Lc2, 25)

factory.addLine(1 , 17, 1)
factory.addLine(17, 16, 2)

factory.addCircleArc(14,15,16, 3)
factory.addLine(14,13, 4)
factory.addLine(13,12, 5)
factory.addLine(12,11, 6)
factory.addLine(11,10, 7)
factory.addCircleArc(8,9,10, 8)
factory.addLine(8,7, 9)
factory.addLine(7,6, 10)
factory.addLine(6,5, 11)
factory.addCircleArc(3,4,5, 12)
factory.addLine(3,2, 13)
factory.addLine(2,1, 14)
factory.addLine(18,19, 15)
factory.addCircleArc(21,20,24, 16)
factory.addCircleArc(24,20,19, 17)
factory.addCircleArc(18,23,25, 18)
factory.addCircleArc(25,23,22, 19)
factory.addLine(21,22, 20)

factory.addCurveLoop([17,-15,18,19,-20,16], 21)
factory.addPlaneSurface([21], 22)
factory.addCurveLoop([11,-12,13,14,1,2,-3,4,5,6,7,-8,9,10], 23)

# A surface with one hole is specified using 2 curve loops:
factory.addPlaneSurface([23,21], 24)

# FIXME: this will be implemented through the gmshView API
#  View "comments" {
#    T2(10, -10, 0){ StrCat("Created on ", Today, " with Gmsh") };
#    T3(0, 0.11, 0, TextAttributes("Align", "Center", "Font", "Helvetica")){ "Hole" };
#    T3(0, 0.09, 0, TextAttributes("Align", "Center")){ "file://image.png@0.01x0" };
#    T3(-0.01, 0.09, 0, 0){ "file://image.png@0.01x0,0,0,1,0,1,0" };
#    T3(0, 0.12, 0, TextAttributes("Align", "Center")){ "file://image.png@0.01x0#" };
#    T2(350, -7, 0){ "file://image.png@20x0" };
# };

factory.synchronize()

model.mesh.generate(2)

gmsh.write("t4.msh")

gmsh.finalize()
