import gmsh
import sys

gmsh.initialize(sys.argv)
gmsh.model.add("normals")
gmsh.model.occ.addSphere(0, 0, 0, 1)
gmsh.model.occ.addBox(2, 0, 0, 1, 1, 1)
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)

nn = []
cc = []

# get all surfaces
ent = gmsh.model.getEntities(2)

for e in ent:
    surf = e[1]
    # get nodes on surf, including those on the boundary (contrary to internal
    # nodes, which store their parametric coordinates, boundary nodes will be
    # reparametrized on surf in order to compute their parametric coordinates -
    # which will be different when reparametrized on another adjacent surface)
    tags, coord, param = gmsh.model.mesh.getNodes(2, surf, True)
    # get surface normal on all nodes, i.e. including on the geometrical
    # singularities (edges/points)
    normals = gmsh.model.getNormal(surf, param)
    # get surface curvature
    curv = gmsh.model.getCurvature(2, surf, param)
    for i in range(0,len(coord),3):
        nn.append(coord[i])
        nn.append(coord[i+1])
        nn.append(coord[i+2])
        nn.append(normals[i])
        nn.append(normals[i+1])
        nn.append(normals[i+2])
        cc.append(coord[i])
        cc.append(coord[i+1])
        cc.append(coord[i+2])
        cc.append(curv[i/3])

t = gmsh.view.add("normals")
gmsh.view.addListData(t, "VP", len(nn)/6, nn)
gmsh.view.write(t, "normals.pos")

t = gmsh.view.add("curvatures")
gmsh.view.addListData(t, "SP", len(cc)/4, cc)
gmsh.view.write(t, "curvatures.pos")

gmsh.finalize()
