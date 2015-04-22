#!/usr/bin/env python
"""
Simple main SfePy scripts wrapper.
"""

from __future__ import print_function

import os.path as op

import glob

import argparse
import subprocess

import sfepy


def get_commands():
    """
    Get available commands (AKA scripts) to run via main SfePy wrapper.

    Commands are dynamically defined by presence in pre-defined directories:
        - <sfeppy.data_dir>/scripts-common
        - user defined directory (not implemented)

    TBD:
        - Windows UNC-paths may be not working correctly (according os.path doc)
        - Check for mandatory/default scripts existence?

    :rtype : dict { command: path_to_script }
    """

    bin_dir = 'scripts-common'  # TBD: Get actual values from SfePy
    if not sfepy.in_source_tree:
        bin_dir = op.normpath(op.join(sfepy.data_dir, bin_dir))

    scripts = glob.glob(op.normpath(op.join(bin_dir, '*.py')))
    cmd = [op.splitext(op.basename(i))[0] for i in scripts]

    commands = dict(zip(cmd, scripts))

    return commands


def main():
    """
    TBD
    :rtype : object
    """

    cmd_list = get_commands()

    parser = argparse.ArgumentParser(
        description='Main Sfepy commands wrapper.',
        version='%(prog)s' + sfepy.__version__,
        usage='%(prog)s [command] [options]'
    )

    parser.add_argument(
        '--dir',
        '-d',
        help='Optional scripts directory (not implemented yet).'
    )

    parser.add_argument(
        '--window',
        '-w',
        help='Use alternative pythonw interpreter.',
        action='store_false',
        dest='py_cmd'
    )

    parser.add_argument(
        'command',
        choices=cmd_list.keys(),
        help='Available SfePy command(s).')

    parser.add_argument(
        'options',
        nargs=argparse.REMAINDER,
        help='Additional options passed directly to selected [command].')

    py_cmd = 'python' if parser.parse_args().py_cmd else 'pythonw'
    command = parser.parse_args().command
    options = ', '.join(parser.parse_args().options)

    args = [py_cmd, cmd_list[command], options]

    print(options)
    print(py_cmd + cmd_list[command], options)
    print(args)

    subprocess.call(args)

if __name__ == '__main__':
    main()
