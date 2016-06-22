#!/usr/bin/env python
"""
Print various information about a mesh.
"""
from __future__ import absolute_import
from argparse import RawDescriptionHelpFormatter, ArgumentParser

from sfepy.base.base import output
from sfepy.discrete.fem import Mesh

helps = {
    'filename' :
    'mesh file name',
}

def main():
    parser = ArgumentParser(description=__doc__.rstrip(),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('filename', help=helps['filename'])
    options = parser.parse_args()

    mesh = Mesh.from_file(options.filename)

    output(mesh.cmesh)
    output('element types:', mesh.descs)
    output('nodal BCs:', mesh.nodal_bcs)

    bbox = mesh.get_bounding_box()
    output('bounding box: %s'
           % ', '.join('%s: [%s, %s]' % (name, bbox[0, ii], bbox[1, ii])
                       for ii, name in enumerate('xyz'[:mesh.dim])))

    output('centre:', mesh.coors.mean(0))

if __name__ == '__main__':
    main()
