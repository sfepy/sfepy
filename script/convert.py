#!/usr/bin/env python
import os, sys

import sfepy.mesh as geom
from sfepy.fem import Mesh
try:
    from site_cfg import tetgen_path
except ImportError:
    tetgen_path = '/usr/bin/tetgen'

def mesh():
    if len( sys.argv ) == 3:
        geom_filename = sys.argv[1]
        vtk_filename = sys.argv[2]
    else:
        print 'Usage: %s <gmsh_filename> <mesh_filename>' % sys.argv[0]
        return

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
