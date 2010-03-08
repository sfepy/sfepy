#!/usr/bin/env python
import os.path as op
from optparse import OptionParser

from sfepy.base.base import *
from sfepy.fem import Mesh, Domain
from sfepy.mechanics.friction import define_dual_mesh

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

    reg_omega = domain.create_region('Omega', 'all')
    reg = domain.create_region('Surface',
                               'nodes of surface',
                               {'can_cells' : True})

    dual_mesh = define_dual_mesh(reg, reg_omega)

if __name__ == '__main__':
    main()
