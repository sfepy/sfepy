#!/usr/bin/env python

import os
import pexpect

import geom
from sfe.fem.mesh import Mesh

def mesh():
    meshgeometry="database/box.geo"
    pexpect.run("gmsh -0 %s -o tmp/x.geo"%(meshgeometry))
    g=geom.read_gmsh("tmp/x.geo")
    g.printinfo()
    geom.write_tetgen(g,"tmp/t.poly")
    geom.runtetgen("tmp/t.poly",a=100.0,Q=1.0,quadratic=False)

    m = Mesh.fromFile("tmp/t.1.node")
    m.write("tmp/t.1.vtk", io = "auto")


try:
    os.makedirs( "tmp" )
except OSError, e:
    if e.errno != 17: # [Errno 17] File exists
        raise
mesh()
