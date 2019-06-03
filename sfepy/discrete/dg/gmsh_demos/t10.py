# This file reimplements gmsh/tutorial/t10.geo in Python.

import gmsh
import math

model = gmsh.model
factory = model.geo

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

model.add("t10")

lc = .15
factory.addPoint(0.0,0.0,0, lc, 1)
factory.addPoint(1,0.0,0, lc, 2)
factory.addPoint(1,1,0, lc, 3)
factory.addPoint(0,1,0, lc, 4)
factory.addPoint(0.2,.5,0, lc, 5)

factory.addLine(1,2, 1);
factory.addLine(2,3, 2);
factory.addLine(3,4, 3);
factory.addLine(4,1, 4);

factory.addCurveLoop([1,2,3,4], 5)
factory.addPlaneSurface([5], 6)

model.mesh.field.add("Attractor", 1)
model.mesh.field.setNumbers(1, "NodesList", [5])
model.mesh.field.setNumber(1, "NNodesByEdge", 100)
model.mesh.field.setNumbers(1, "EdgesList", [2])

model.mesh.field.add("Threshold", 2);
model.mesh.field.setNumber(2, "IField", 1);
model.mesh.field.setNumber(2, "LcMin", lc / 30)
model.mesh.field.setNumber(2, "LcMax", lc)
model.mesh.field.setNumber(2, "DistMin", 0.15)
model.mesh.field.setNumber(2, "DistMax", 0.5)

model.mesh.field.add("MathEval", 3)
model.mesh.field.setString(3, "F", "Cos(4*3.14*x) * Sin(4*3.14*y) / 10 + 0.101")

model.mesh.field.add("Attractor", 4)
model.mesh.field.setNumbers(4, "NodesList", [1])

model.mesh.field.add("MathEval", 5);
model.mesh.field.setString(5, "F", "F4^3 + " + str(lc / 100))

model.mesh.field.add("Box", 6)
model.mesh.field.setNumber(6, "VIn", lc / 15)
model.mesh.field.setNumber(6, "VOut", lc)
model.mesh.field.setNumber(6, "XMin", 0.3)
model.mesh.field.setNumber(6, "XMax", 0.6)
model.mesh.field.setNumber(6, "YMin", 0.3)
model.mesh.field.setNumber(6, "YMax", 0.6)

model.mesh.field.add("Min", 7)
model.mesh.field.setNumbers(7, "FieldsList", [2, 3, 5, 6])

model.mesh.field.setAsBackgroundMesh(7)

factory.synchronize()
model.mesh.generate(2)
gmsh.write("t10.msh")
gmsh.finalize()

