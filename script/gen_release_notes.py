#!/usr/bin/env python
"""
Generate release notes using git log starting from the given version.
"""
from __future__ import print_function
from argparse import ArgumentParser

import sys
import subprocess
from textwrap import wrap

def main():
    parser = ArgumentParser(description=__doc__.rstrip())
    parser.add_argument('version')
    options = parser.parse_args()

    print(options.version)

    cmd = """git log --pretty=format:"%%s%%n%%b%%n" --topo-order --reverse release_%s..HEAD""" % options.version

    print(cmd)

    if sys.version_info > (3, 0):
        raw = subprocess.check_output(cmd.split(), encoding='utf-8')

    else:
        raw = subprocess.check_output(cmd.split())

    msgs = raw.split('\n"\n')

    merges = []
    ims = [0]
    for msg in msgs:
        msg = msg.replace('\n"', '')
        if 'Merge' in msg:
            merges.insert(ims[-1], msg)
            ims.append(len(merges))

        else:
            merges.append(msg)

    for ii in range(len(ims) - 1):
        merge = merges[ims[ii]:ims[ii+1]]

        msg = merge[0][1:].strip().lower()
        wmsg = '\n  '.join(wrap(msg, 77))

        print('- ' + wmsg)
        print()

        for msg in merge[1:]:
            submsgs = msg[1:].split('\n- ')

            wmsg = '\n    '.join(wrap(submsgs[0].strip(), 75))
            print('  - ' + wmsg)
            for submsg in submsgs[1:]:
                wmsg = '\n    '.join(wrap(submsg.strip(), 73))
                print('    - ' + wmsg)

        print()

if __name__ == '__main__':
    main()
