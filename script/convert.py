#!/usr/bin/env python
import os, sys

import sfepy.geom as geom
from sfepy.fem import Mesh
try:
    from site_cfg import tetgen_path
except ImportError:
    tetgen_path = '/usr/bin/tetgen'

def mesh():
    if len( sys.argv ) == 3:
        geom_filename = sys.argv[1]
        vtk_filename = sys.argv[2]
    if len( sys.argv ) == 2:
        geom_filename = sys.argv[1]
        vtk_filename = "tmp/t.1.vtk"
    else:
        geom_filename = "database/box.geo"
        vtk_filename = "tmp/t.1.vtk"

    os.system( "gmsh -0 %s -o tmp/x.geo" % geom_filename )
    g=geom.read_gmsh("tmp/x.geo")
    g.printinfo()
    geom.write_tetgen(g,"tmp/t.poly")
    geom.runtetgen("tmp/t.poly",a=0.03,Q=1.0,quadratic=False,
                   tetgenpath = tetgen_path)

    m = Mesh.from_file("tmp/t.1.node")
    m.write( vtk_filename, io = "auto" )

try:
    os.makedirs( "tmp" )
except OSError, e:
    if e.errno != 17: # [Errno 17] File exists
        raise

mesh()
