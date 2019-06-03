import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("periodic")
R = 2
gmsh.model.occ.addBox(0, 0, 0, R, R, R)
gmsh.model.occ.synchronize()

ent = gmsh.model.getEntities(0);
gmsh.model.mesh.setSize(ent, 1);
gmsh.model.mesh.setSize([(0,1)], 0.01);
gmsh.model.mesh.setPeriodic(2, [2], [1], [1,0,0,R, 0,1,0,0, 0,0,1,0, 0,0,0,1])

gmsh.model.mesh.generate(2)
#gmsh.model.mesh.setOrder(2)

masterTag, nodes, tfo = gmsh.model.mesh.getPeriodicNodes(2, 2)
# print masterTag, nodes, tfo

gmsh.write("periodic.msh")
