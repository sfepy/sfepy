#!/usr/bin/env python
"""
Plot mesh connectivities, facet orientations, global and local DOF ids etc.
"""
from __future__ import absolute_import
import sys
sys.path.append('.')
from argparse import ArgumentParser

import matplotlib.pyplot as plt

from sfepy.base.base import output
from sfepy.discrete.fem import Mesh, FEDomain
import sfepy.postprocess.plot_cmesh as pc

helps = {
}

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('filename')
    options = parser.parse_args()

    filename = options.filename

    mesh = Mesh.from_file(filename)
    output('Mesh:')
    output('  dimension: %d, vertices: %d, elements: %d'
           % (mesh.dim, mesh.n_nod, mesh.n_el))

    domain = FEDomain('domain', mesh)
    output(domain.cmesh)
    domain.cmesh.cprint(1)
    dim = domain.cmesh.dim

    entities_opts = [
        {'color' : 'k', 'label_global' : 12, 'label_local' : 8},
        {'color' : 'b', 'label_global' : 12, 'label_local' : 8},
        {'color' : 'g', 'label_global' : 12, 'label_local' : 8},
        {'color' : 'r', 'label_global' : 12},
    ]
    if dim == 2: entities_opts.pop(2)

    pc.plot_cmesh(None, domain.cmesh,
                  wireframe_opts = {'color' : 'k'},
                  entities_opts=entities_opts)

    plt.show()

if __name__ == '__main__':
    main()
