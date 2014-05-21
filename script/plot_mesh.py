#!/usr/bin/env python
"""
Plot mesh connectivities, facet orientations, global and local DOF ids etc.
"""
from optparse import OptionParser

from sfepy.base.base import output
from sfepy.discrete.fem import Mesh, FEDomain
import sfepy.postprocess.plot_cmesh as pc

usage = '%prog [options] filename\n' + __doc__.rstrip()

helps = {
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    options, args = parser.parse_args()

    if len(args) == 1:
        filename = args[0]
    else:
        parser.print_help(),
        return

    mesh = Mesh.from_file(filename)
    output('Mesh:')
    output('  dimension: %d, vertices: %d, elements: %d'
           % (mesh.dim, mesh.n_nod, mesh.n_el))

    domain = FEDomain('domain', mesh)
    output(domain.cmesh)
    domain.cmesh.cprint(1)
    dim = domain.cmesh.dim

    ax = pc.plot_wireframe(None, domain.cmesh)

    ax = pc.plot_entities(ax, domain.cmesh, 0, 'k')
    ax = pc.label_global_entities(ax, domain.cmesh, 0, 'k', 12)
    ax = pc.label_local_entities(ax, domain.cmesh, 0, 'k', 8)

    ax = pc.plot_entities(ax, domain.cmesh, 1, 'b')
    ax = pc.label_global_entities(ax, domain.cmesh, 1, 'b', 12)
    ax = pc.label_local_entities(ax, domain.cmesh, 1, 'b', 8)

    if dim == 3:
        ax = pc.plot_entities(ax, domain.cmesh, 2, 'g')
        ax = pc.label_global_entities(ax, domain.cmesh, 2, 'g', 12)
        ax = pc.label_local_entities(ax, domain.cmesh, 2, 'g', 8)

    ax = pc.plot_entities(ax, domain.cmesh, dim, 'r')
    ax = pc.label_global_entities(ax, domain.cmesh, dim, 'r', 12)

    pc.plt.show()

if __name__ == '__main__':
    main()
