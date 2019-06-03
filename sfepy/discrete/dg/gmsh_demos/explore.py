import gmsh
import sys

if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " file.msh [options]"
    exit(0)

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.open(sys.argv[1])

# get all elementary entities in the model
entities = gmsh.model.getEntities()

for e in entities:
    # get the mesh nodes for each elementary entity
    nodeTags, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(e[0], e[1])
    # get the mesh elements for each elementary entity
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(e[0], e[1])
    # report some statistics
    numElem = sum(len(i) for i in elemTags)
    print str(len(nodeTags)) + " mesh nodes and " + str(numElem),\
          "mesh elements on entity " + str(e) + " " + gmsh.model.getType(e[0], e[1])
    for t in elemTypes:
        name, dim, order, numv, parv = gmsh.model.mesh.getElementProperties(t)
        print " - Element type: " + name + ", order " + str(order)
        print "   with " + str(numv) + " nodes in param coord: ", parv

# all mesh node coordinates
nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
x = dict(zip(nodeTags, nodeCoords[0::3]))
y = dict(zip(nodeTags, nodeCoords[1::3]))
z = dict(zip(nodeTags, nodeCoords[2::3]))

gmsh.finalize()
