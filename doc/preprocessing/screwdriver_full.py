import sys
FREECADPATH = '/usr/lib/freecad/lib/'
sys.path.append(FREECADPATH)

from FreeCAD import Base, newDocument
import Part
import Draft
import ProfileLib.RegularPolygon as Poly
import MeshPart


doc = newDocument()

radius = 0.01
height = 0.1

cyl = doc.addObject("Part::Cylinder", "cyl")
cyl.Radius = radius
cyl.Height = height

sph = doc.addObject("Part::Sphere", "sph")
sph.Radius = radius

uni = doc.addObject("Part::MultiFuse", "uni")
uni.Shapes = [cyl, sph]

ske = doc.addObject('Sketcher::SketchObject', 'Sketch')
ske.Placement = Base.Placement(Base.Vector(0, 0, 0),
                               Base.Rotation(-0.707107, 0, 0, -0.707107))
Poly.makeRegularPolygon('Sketch', 5,
                        Base.Vector(-1.2 * radius, 0.9 * height, 0),
                        Base.Vector(-0.8 * radius, 0.9 * height, 0))

cut = doc.addObject("PartDesign::Revolution", "Revolution")
cut.Sketch = ske
cut.ReferenceAxis = (ske, ['V_Axis'])
cut.Angle = 360.0

dif = doc.addObject("Part::Cut", "dif")
dif.Base = uni
dif.Tool = cut

cyl1 = doc.addObject("Part::Cylinder", "cyl1")
cyl1.Radius = 0.2 * radius
cyl1.Height = 1.1 * height
cyl1.Placement = Base.Placement(Base.Vector(-1.1 * radius, 0, -0.2 * height),
                                Base.Rotation(0, 0, 0, 1))

arr = Draft.makeArray(cyl1, Base.Vector(1, 0, 0), Base.Vector(0, 1, 0), 2, 2)
arr.ArrayType = "polar"
arr.NumberPolar = 6

dif2 = doc.addObject("Part::Cut", "dif2")
dif2.Base = dif
dif2.Tool = arr

cyl2 = doc.addObject("Part::Cylinder", "cyl2")
cyl2.Radius = 0.3 * radius
cyl2.Height = height

dif3 = doc.addObject("Part::Cut", "dif3")
dif3.Base = dif2
dif3.Tool = cyl2

cyl3 = doc.addObject("Part::Cylinder", "cyl3")
cyl3.Radius = 0.3 * radius
cyl3.Height = height
cyl3.Placement = Base.Placement(Base.Vector(0, 0, height),
                                Base.Rotation(0, 0, 0, 1))

tip1 = doc.addObject("Part::Box", "tip1")
tip1.Length = radius
tip1.Width = 2 * radius
tip1.Height = 3 * radius
tip1.Placement = Base.Placement(Base.Vector(0, -radius, 1.71 * height),
                                Base.Rotation(Base.Vector(0, 1, 0), -10),
                                Base.Vector(0, 0, 3 * radius))

tip2 = doc.addObject("Part::Mirroring", "tip2")
tip2.Source = tip1
tip2.Normal = (1, 0, 0)

tip3 = doc.addObject("Part::MultiFuse", "tip3")
tip3.Shapes = [tip1, tip2]

dif4 = doc.addObject("Part::Cut", "dif4")
dif4.Base = cyl3
dif4.Tool = tip3

uni2 = doc.addObject("Part::MultiFuse", "uni2")
uni2.Shapes = [cyl2, dif4]

doc.recompute()

Part.export([dif3, uni2], 'screwdriver_full.step')
doc.saveAs('screwdriver_full.FCStd')
