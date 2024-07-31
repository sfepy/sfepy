#!/usr/bin/env python
"""
Print various information about a mesh.
"""
import sys
sys.path.append('.')
from argparse import RawDescriptionHelpFormatter, ArgumentParser

import numpy as nm

from sfepy.base.base import output
from sfepy.discrete.fem import Mesh, MeshIO, FEDomain
from sfepy.discrete.common.extmods.cmesh import graph_components

def show_mesh_info(options):
    mesh = Mesh.from_file(options.filename)

    output(mesh.cmesh)
    output('element types:', mesh.descs)
    output('nodal BCs:', sorted(mesh.nodal_bcs.keys()))
    output('vertex groups:', nm.unique(mesh.cmesh.vertex_groups))
    output('cell groups:', nm.unique(mesh.cmesh.cell_groups))

    bbox = mesh.get_bounding_box()
    output('bounding box:\n%s'
           % '\n'.join('%s: [%14.7e, %14.7e]' % (name, bbox[0, ii], bbox[1, ii])
                       for ii, name in enumerate('xyz'[:mesh.dim])))
    output('box dimensions:   [%s]'
           % ', '.join('%14.7e' % ii for ii in (bbox[1] - bbox[0])))
    output('centre:           [%s]'
           % ', '.join('%14.7e' % ii for ii in 0.5 * (bbox[0] + bbox[1])))
    output('coordinates mean: [%s]'
           % ', '.join('%14.7e' % ii for ii in mesh.coors.mean(0)))

    io = MeshIO.any_from_filename(options.filename)
    try:
        data = io.read_data(0)

    except ValueError:
        pass

    else:
        from sfepy.linalg.utils import output_array_stats

        output('data...')
        for key, val in data.items():
            if key in ('node_groups', 'mat_id'): continue
            output_array_stats(val.data, f'{key} ({val.mode}):')
        output('...done')

    if not options.detailed: return

    domain = FEDomain(mesh.name, mesh)

    for dim in range(1, mesh.cmesh.tdim + 1):
        volumes = mesh.cmesh.get_volumes(dim)
        output('volumes of %d %dD entities:\nmin: %.7e mean: %.7e median:'
               ' %.7e max: %.7e'
               % (mesh.cmesh.num[dim], dim, volumes.min(), volumes.mean(),
                  nm.median(volumes), volumes.max()))

    euler = lambda mesh: nm.dot(mesh.cmesh.num, [1, -1, 1, -1])
    ec = euler(mesh)
    output('Euler characteristic:', ec)

    graph = mesh.create_conn_graph(verbose=False)
    n_comp, _ = graph_components(graph.shape[0], graph.indptr, graph.indices)
    output('number of connected components:', n_comp)

    if mesh.dim > 1:
        region = domain.create_region('surf', 'vertices of surface', 'facet')
        surf_mesh = Mesh.from_region(region, mesh,
                                     localize=True, is_surface=True)
        FEDomain(surf_mesh.name, surf_mesh) # Calls CMesh.setup_entities().

        sec = euler(surf_mesh)
        output('surface Euler characteristic:', sec)
        if mesh.dim == 3:
            output('surface genus:', (2.0 - sec) / 2.0)

        surf_graph = surf_mesh.create_conn_graph(verbose=False)
        n_comp, _ = graph_components(surf_graph.shape[0],
                                     surf_graph.indptr, surf_graph.indices)
        output('number of connected surface components:', n_comp)

helps = {
    'filename' :
    'mesh file name',
    'detailed' :
    'show additional information (entity volume statistics)',
}

def add_args(parser):
    parser.add_argument('filename', help=helps['filename'])
    parser.add_argument('-d', '--detailed',
                        action='store_true', dest='detailed',
                        default=False, help=helps['detailed'])

def main():
    parser = ArgumentParser(description=__doc__.rstrip(),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s')
    add_args(parser)

    options = parser.parse_args()
    show_mesh_info(options)

if __name__ == '__main__':
    main()
