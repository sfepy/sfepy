#!/usr/bin/env python
from __future__ import print_function
import shlex
import subprocess

def main():
    args = shlex.split('git log --pretty=format:"%an <%ae>"')
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    out = p.communicate()[0].split('\n')

    done = set()
    unique = []
    counts = {}
    for line in reversed(out):
        line = line.strip()
        if not len(line): continue

        name, email = line.split(' <')
        name = name.strip()
        email = '<' + email.strip()

        if not name in done:
            done.add(name)
            unique.append(name)

        record = counts.setdefault(name, [0, email])
        record[0] += 1
        if not email in record[1:]:
            record.append(email)


    print('List of contributors ordered by date of the first commit,'
          ' with commit counts:')

    for line in unique:
        record = counts[line]
        print(('%6d %s %s' % (record[0], line, ', '.join(record[1:]))))

if __name__ == '__main__':
    main()
