#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
import shlex
import subprocess

def main():
    args = shlex.split('git log --pretty=format:"%an <%ae>"')
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    out = p.communicate()[0].decode('utf-8').split('\n')

    done = set()
    e2n = {}
    unique = []
    counts = {}
    for line in reversed(out):
        line = line.strip()
        if not len(line): continue

        name, email = line.split(' <')
        name = name.strip()
        email = '<' + email.strip()

        if not name in done:
            if email not in e2n:
                done.add(name)
                e2n[email] = name

                unique.append(name)

            elif name not in counts[e2n[email]]:
                counts[e2n[email]].append(name)

        name = e2n.get(email, name)

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
