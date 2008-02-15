import os
import pexpect

import geom

def mesh():
    meshgeometry="database/t1.geo.sphere-tri"
    pexpect.run("gmsh -0 %s -o tmp/x.geo"%(meshgeometry))
    g=geom.read_gmsh("tmp/x.geo")
    g.printinfo()
    geom.write_tetgen(g,"tmp/t.poly")
    geom.runtetgen("tmp/t.poly",a=0.3,Q=1.0,quadratic=False)
    m=geom.read_tetgen("tmp/t.1")
    m.printinfo()
    m.writemsh("tmp/t12.msh")

def mesh2():
    from sfe.fem.mesh import Mesh
    m = Mesh.fromFile("tmp/t.1.node")
    m.write("tmp/t.1.vtk", io = "auto")


try:
    os.makedirs( "tmp" )
except OSError, e:
    if e.errno != 17: # [Errno 17] File exists
        raise
#mesh()
mesh2()
