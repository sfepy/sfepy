#!/usr/bin/env python
import os.path as op
from optparse import OptionParser

from sfepy.base.base import *
from sfepy.fem import Mesh, Domain
from sfepy.mechanics.friction import DualMesh

usage = """%prog [options] mesh_filename"""

def main():
    parser = OptionParser(usage=usage, version="%prog")

    options, args = parser.parse_args()

    if (len(args) == 1):
        mesh_filename = args[0];
    else:
        parser.print_help(),
        return

    mesh = Mesh('mesh', mesh_filename)
    print mesh

    domain = Domain('domain', mesh)
    print domain

    reg = domain.create_region('Surface',
                               'nodes of surface',
                               {'can_cells' : True})

    dual_mesh = DualMesh(reg)
    dual_mesh.save('dual_mesh.mesh',)
    dual_mesh.save_axes('axes.vtk',)

    print dual_mesh

if __name__ == '__main__':
    main()
