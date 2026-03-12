#!/usr/bin/env python
"""
Simple mesh generators and statistics.
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from . import (blockgen, cylindergen, gen_iga_patch, combine_meshes,
               show_mesh_info)

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='valid subcommands',
                                       help='additional help',
                                       dest='mesh_kind',
                                       required=True)

    sub = subparsers.add_parser('block', help='generate a block mesh',
                                description=blockgen.__doc__,
                                formatter_class=RawDescriptionHelpFormatter)
    sub.set_defaults(fun=blockgen.gen_block)
    blockgen.add_args(sub)

    sub = subparsers.add_parser('cylinder', help='generate a cylinder mesh',
                                description=cylindergen.__doc__,
                                formatter_class=RawDescriptionHelpFormatter)
    sub.set_defaults(fun=cylindergen.gen_cylinder)
    cylindergen.add_args(sub)

    sub = subparsers.add_parser('iga-patch', help='generate an IGA patch',
                                description=gen_iga_patch.__doc__,
                                formatter_class=RawDescriptionHelpFormatter)
    sub.set_defaults(fun=gen_iga_patch.gen_iga_patch)
    gen_iga_patch.add_args(sub)

    sub = subparsers.add_parser('combine', help='combine meshes',
                                description=combine_meshes.__doc__,
                                formatter_class=RawDescriptionHelpFormatter)
    sub.set_defaults(fun=combine_meshes.combine_meshes)
    combine_meshes.add_args(sub)

    sub = subparsers.add_parser('info', help='show various mesh statistics',
                                description=show_mesh_info.__doc__,
                                formatter_class=RawDescriptionHelpFormatter)
    sub.set_defaults(fun=show_mesh_info.show_mesh_info)
    show_mesh_info.add_args(sub)

    options = parser.parse_args()

    options.fun(options)

if __name__ == '__main__':
    main()
