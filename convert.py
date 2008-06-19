#!/usr/bin/env python
import os, sys
import pexpect

import geom
from sfe.fem.mesh import Mesh
try:
    from site_cfg import tetgen_path
except ImportError:
    tetgen_path = '/usr/bin/tetgen'

def mesh():
    if len( sys.argv ) == 3:
        geomFileName = sys.argv[1]
        vtkFileName = sys.argv[2]
    if len( sys.argv ) == 2:
        geomFileName = sys.argv[1]
        vtkFileName = "tmp/t.1.vtk"
    else:
        geomFileName = "database/box.geo"
        vtkFileName = "tmp/t.1.vtk"

    pexpect.run( "gmsh -0 %s -o tmp/x.geo" % geomFileName )
    g=geom.read_gmsh("tmp/x.geo")
    g.printinfo()
    geom.write_tetgen(g,"tmp/t.poly")
    geom.runtetgen("tmp/t.poly",a=0.0003,Q=1.0,quadratic=False,
                   tetgenpath = tetgen_path)

    m = Mesh.fromFile("tmp/t.1.node")
    m.write( vtkFileName, io = "auto" )

try:
    os.makedirs( "tmp" )
except OSError, e:
    if e.errno != 17: # [Errno 17] File exists
        raise

mesh()
