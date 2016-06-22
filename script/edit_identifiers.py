#!/usr/bin/env python
"""
Convert mixedCase identifiers to under_scores.
"""
from __future__ import print_function
from __future__ import absolute_import
import sys, re, os

##
# Taken from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66009.
def cw2us(x): # capwords to underscore notation
#    return re.sub(r'(?<=[a-z])[A-Z]', r"_\g<0>", x).lower()
    return re.sub(r'(?<=[a-z])[A-Z]|(?<!^)[A-Z](?=[a-z])', r"_\g<0>", x).lower()

def mc2us(x): # mixed case to underscore notation
    return cw2us(x)

def us2mc(x): # underscore to mixed case notation
    return re.sub(r'_([a-z])', lambda m: (m.group(1).upper()), x)

def us2cw(x): # underscore to capwords notation
    s = us2mc(x)
    return s[0].upper()+s[1:]
##

misc = '[\'\"\[\]\(\)\{\}]'
mixed = '^%s*(?=[^A-Z])[a-z]+' % misc
match_candidate = re.compile(mixed).match

def split_on(token, chars):
    if len(chars) == 1:
        return token.split(chars[0])
    else:
        aux = token.split(chars[0])
        out = []
        for item in aux:
            out.extend(split_on(item, chars[1:]))
        return out

def edit(line):
    count = 0
    aux = line.split()
    if len(aux) > 2:
        if aux[0] == 'from' and aux[2] == 'import':
            aux = aux[3:]
        elif aux[0] == 'import':
            aux = aux[3:]

    for token in aux:
        for item in split_on(token, '.,=*/+-_'):
            if match_candidate(item):
                line = line.replace(item, cw2us(item), 1)
                count += 1
    return line, count

def main():

    write = True

    for name in sys.argv[1:]:
        print(name)

        rfd = open(name, 'r')
        path = os.path.dirname(name)
        base = os.path.basename(name)
        new_name = os.path.join(path, 'new_' + base)
        if write:
            wfd = open(new_name, 'w')

        n_edit = 0
        for line in rfd:
            eline, count = edit(line)
            if write:
                wfd.write(eline)
            n_edit += count

        rfd.close()

        print('%d edit candidates' % n_edit)

        if write:
            wfd.close()
            os.rename(new_name, name)

if __name__ == '__main__':
    main()
