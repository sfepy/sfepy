#!/usr/bin/env python
"""
Print various information about a mesh.
"""
from __future__ import absolute_import
from argparse import RawDescriptionHelpFormatter, ArgumentParser

import numpy as nm

from sfepy.base.base import output
from sfepy.discrete.fem import Mesh

helps = {
    'filename' :
    'mesh file name',
    'detailed' :
    'show additional information (entity volume statistics)',
}

def main():
    parser = ArgumentParser(description=__doc__.rstrip(),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('filename', help=helps['filename'])
    parser.add_argument('-d', '--detailed',
                        action='store_true', dest='detailed',
                        default=False, help=helps['detailed'])
    options = parser.parse_args()

    mesh = Mesh.from_file(options.filename)

    output(mesh.cmesh)
    output('element types:', mesh.descs)
    output('nodal BCs:', sorted(mesh.nodal_bcs.keys()))

    bbox = mesh.get_bounding_box()
    output('bounding box: %s'
           % ', '.join('%s: [%s, %s]' % (name, bbox[0, ii], bbox[1, ii])
                       for ii, name in enumerate('xyz'[:mesh.dim])))

    output('centre:', 0.5 * (bbox[0] + bbox[1]))
    output('coordinates mean:', mesh.coors.mean(0))

    if not options.detailed: return

    from sfepy.discrete.fem.geometry_element import create_geometry_elements
    gels = create_geometry_elements()
    mesh.cmesh.set_local_entities(gels)
    mesh.cmesh.setup_entities()

    for dim in range(1, mesh.cmesh.tdim + 1):
        volumes = mesh.cmesh.get_volumes(dim)
        output('volumes of %d %dD entities: min: %s mean: %s max: %s'
               % (mesh.cmesh.num[dim],
                  dim, volumes.min(), volumes.mean(), volumes.max()))

    output('Euler characteristic:', nm.dot(mesh.cmesh.num, [1, -1, 1, -1]))

if __name__ == '__main__':
    main()
