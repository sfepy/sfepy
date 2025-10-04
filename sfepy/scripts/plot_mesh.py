#!/usr/bin/env python
"""
Plot mesh connectivities, facet orientations, global and local DOF ids etc.

To switch off plotting some mesh entities, set the corresponding color to
`None`.
"""
import sys
sys.path.append('.')
from argparse import ArgumentParser

import matplotlib.pyplot as plt

from sfepy.base.base import output
from sfepy.base.conf import dict_from_string
from sfepy.discrete.fem import Mesh, FEDomain
import sfepy.postprocess.plot_cmesh as pc

helps = {
    'vertex_opts' : 'plotting options for mesh vertices'
    ' [default: %(default)s]',
    'edge_opts' : 'plotting options for mesh edges'
    ' [default: %(default)s]',
    'face_opts' : 'plotting options for mesh faces'
    ' [default: %(default)s]',
    'cell_opts' : 'plotting options for mesh cells'
    ' [default: %(default)s]',
    'wireframe_opts' : 'plotting options for mesh wireframe'
    ' [default: %(default)s]',
    'no_axes' :
    'do not show the figure axes',
    'no_show' :
    'do not show the mesh plot figure',
}

def main():
    default_vertex_opts = """color='k', label_global=12,
                             label_local=8"""
    default_edge_opts = """color='b', label_global=12,
                           label_local=8"""
    default_face_opts = """color='g', label_global=12,
                           label_local=8"""
    default_cell_opts = """color='r', label_global=12"""
    default_wireframe_opts = "color='k'"

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('--vertex-opts', metavar='dict-like',
                        action='store', dest='vertex_opts',
                        default=default_vertex_opts,
                        help=helps['vertex_opts'])
    parser.add_argument('--edge-opts', metavar='dict-like',
                        action='store', dest='edge_opts',
                        default=default_edge_opts,
                        help=helps['edge_opts'])
    parser.add_argument('--face-opts', metavar='dict-like',
                        action='store', dest='face_opts',
                        default=default_face_opts,
                        help=helps['face_opts'])
    parser.add_argument('--cell-opts', metavar='dict-like',
                        action='store', dest='cell_opts',
                        default=default_cell_opts,
                        help=helps['cell_opts'])
    parser.add_argument('--wireframe-opts', metavar='dict-like',
                        action='store', dest='wireframe_opts',
                        default=default_wireframe_opts,
                        help=helps['wireframe_opts'])
    parser.add_argument('--no-axes',
                        action='store_false', dest='axes',
                        help=helps['no_axes'])
    parser.add_argument('-n', '--no-show',
                        action='store_false', dest='show',
                        help=helps['no_show'])
    parser.add_argument('filename')
    parser.add_argument('figname', nargs='?')
    options = parser.parse_args()

    entities_opts = [
        dict_from_string(options.vertex_opts),
        dict_from_string(options.edge_opts),
        dict_from_string(options.face_opts),
        dict_from_string(options.cell_opts),
    ]
    wireframe_opts = dict_from_string(options.wireframe_opts)

    filename = options.filename

    mesh = Mesh.from_file(filename)
    output('Mesh:')
    output('  dimension: %d, vertices: %d, elements: %d'
           % (mesh.dim, mesh.n_nod, mesh.n_el))

    domain = FEDomain('domain', mesh)
    output(domain.cmesh)
    domain.cmesh.cprint(1)
    dim = domain.cmesh.dim

    if dim == 2: entities_opts.pop(2)

    ax = pc.plot_cmesh(None, domain.cmesh,
                       wireframe_opts=wireframe_opts,
                       entities_opts=entities_opts)
    if dim == 2:
        ax.axis('image')

    if not options.axes:
        ax.axis('off')
    plt.tight_layout()

    if options.figname:
        fig = ax.figure
        fig.savefig(options.figname, bbox_inches='tight')

    if options.show:
        plt.show()

if __name__ == '__main__':
    main()
