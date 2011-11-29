#!/usr/bin/env python
import shlex
import subprocess

def main():
    args = shlex.split('git log --pretty=format:"%an <%ae>"')
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    out = p.communicate()[0].split('\n')

    done = set()
    unique = []
    for line in reversed(out):
        line = line.strip()
        if len(line) and not line in done:
            done.add(line)
            unique.append(line)

    for line in unique:
        print line

if __name__ == '__main__':
    main()
