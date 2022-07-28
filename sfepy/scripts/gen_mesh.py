#!/usr/bin/env python
"""
Simple mesh generators.
"""
from argparse import ArgumentParser

from . import blockgen, cylindergen, gen_iga_patch

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='valid subcommands',
                                       help='additional help',
                                       dest='mesh_kind',
                                       required=True)

    parser0 = subparsers.add_parser('block', help='generate a block mesh')
    parser0.set_defaults(fun=blockgen.gen_block)
    blockgen.add_args(parser0)

    parser1 = subparsers.add_parser('cylinder', help='generate a cylinder mesh')
    parser1.set_defaults(fun=cylindergen.gen_cylinder)
    cylindergen.add_args(parser1)

    parser2 = subparsers.add_parser('iga-patch', help='generate an IGA patch')
    parser2.set_defaults(fun=gen_iga_patch.gen_iga_patch)
    gen_iga_patch.add_args(parser2)

    options = parser.parse_args()

    options.fun(options)

if __name__ == '__main__':
    main()
