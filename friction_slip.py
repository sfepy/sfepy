#!/usr/bin/env python
import sys
import os.path as op

from sfepy.base.base import *
from sfepy.fem import Mesh, Domain
from sfepy.mechanics.friction import define_dual_mesh

def main():

    mesh = Mesh('mesh', 'meshes/3d/block.mesh')
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
