#!/usr/bin/env python
import os, sys

import geom
from sfepy.fem.mesh import Mesh
try:
    from site_cfg import tetgen_path
except ImportError:
    tetgen_path = '/usr/bin/tetgen'

def mesh():
    if len( sys.argv ) == 3:
        geom_file_name = sys.argv[1]
        vtk_file_name = sys.argv[2]
    if len( sys.argv ) == 2:
        geom_file_name = sys.argv[1]
        vtk_file_name = "tmp/t.1.vtk"
    else:
        geom_file_name = "database/box.geo"
        vtk_file_name = "tmp/t.1.vtk"

    os.system( "gmsh -0 %s -o tmp/x.geo" % geom_file_name )
    g=geom.read_gmsh("tmp/x.geo")
    g.printinfo()
    geom.write_tetgen(g,"tmp/t.poly")
    geom.runtetgen("tmp/t.poly",a=0.03,Q=1.0,quadratic=False,
                   tetgenpath = tetgen_path)

    m = Mesh.from_file("tmp/t.1.node")
    m.write( vtk_file_name, io = "auto" )

try:
    os.makedirs( "tmp" )
except OSError, e:
    if e.errno != 17: # [Errno 17] File exists
        raise

mesh()
