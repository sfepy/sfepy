import gmsh
import sys

gmsh.initialize(sys.argv)
gmsh.fltk.initialize()

gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("square with cracks")

surf1 = 1
gmsh.model.occ.addRectangle(0, 0, 0, 1, 1, surf1)

pt1 = gmsh.model.occ.addPoint(0.2, 0.2, 0)
pt2 = gmsh.model.occ.addPoint(0.4, 0.4, 0)
line1 = gmsh.model.occ.addLine(pt1, pt2)
pt3 = gmsh.model.occ.addPoint(0.6, 0.1, 0)
pt4 = gmsh.model.occ.addPoint(0.1, 0.3, 0)
line2 = gmsh.model.occ.addLine(pt3, pt4)

o, m = gmsh.model.occ.fragment([(2,surf1)], [(1,line1), (1,line2)])
gmsh.model.occ.synchronize()

# m contains, for each input entity (surf1, line1 and line2), the child entities
# (if any) after the fragmentation, as lists of tuples. To apply the crack
# plugin we group all the intersecting lines in a physical group

new_surf = m[0][0][1]
new_lines = [item[1] for sublist in m[1:] for item in sublist]

gmsh.model.addPhysicalGroup(2, [new_surf], 100)
gmsh.model.addPhysicalGroup(1, new_lines, 101)

gmsh.model.mesh.generate(2)

gmsh.plugin.setNumber("Crack", "PhysicalGroup", 101)
gmsh.plugin.run("Crack")

gmsh.fltk.run()

gmsh.finalize()
