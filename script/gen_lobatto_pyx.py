#!/usr/bin/env python
"""
Generate lobatto.pyx file.
"""
import sys
sys.path.append('.')
import os
from optparse import OptionParser

import sympy as sm
import numpy as nm
import matplotlib.pyplot as plt

from sfepy import top_dir
from sfepy.base.ioutils import InDir

cdef = """
cdef float64 %s(float64 x):
    return %s
"""

fun_list = """
cdef fun %s[%d]

%s[:] = [%s]
"""

def gen_lobbato(max_order):
    assert max_order > 2

    x = sm.symbols('x')

    lobs = [0, 1]
    lobs[0] = (1 - x) / 2
    lobs[1] = (1 + x) / 2

    legs = [sm.legendre(0, 'y')]
    clegs = [sm.ccode(legs[0])]
    clobs = []
    denoms = [] # for lobs.

    for ii, lob in enumerate(lobs):
        clob = sm.ccode(lob)

        clobs.append(clob)

    for ii in range(2, max_order + 1):
        coef = sm.sympify('sqrt(2 * (2 * %s - 1)) / 2' % ii)
        leg = sm.legendre(ii - 1, 'y')

        pleg = leg.as_poly()
        coefs = pleg.all_coeffs()
        denom = max(sm.denom(val) for val in coefs)

        cleg = sm.ccode(sm.horner(leg*denom)/denom)

        lob = sm.simplify(coef * sm.integrate(leg, ('y', -1, x)))
        lobnc = sm.simplify(sm.integrate(leg, ('y', -1, x)))

        plobnc = lobnc.as_poly()
        coefs = plobnc.all_coeffs()
        denom = sm.denom(coef) * max(sm.denom(val) for val in coefs)

        clob = sm.ccode(sm.horner(lob*denom)/denom)

        legs.append(leg)
        clegs.append(cleg)
        lobs.append(lob)
        clobs.append(clob)
        denoms.append(denom)

    coef = sm.sympify('sqrt(2 * (2 * %s - 1)) / 2' % (max_order + 1))
    leg = sm.legendre(max_order, 'y')

    pleg = leg.as_poly()
    coefs = pleg.all_coeffs()
    denom = max(sm.denom(val) for val in coefs)

    cleg = sm.ccode(sm.horner(leg*denom)/denom)
    legs.append(leg)
    clegs.append(cleg)

    kerns = []
    ckerns = []

    for ii, lob in enumerate(lobs[2:]):
        kern = sm.simplify(lob / (lobs[0] * lobs[1]))

        denom = denoms[ii] / 4
        ckern = sm.ccode(sm.horner(kern*denom)/denom)

        kerns.append(kern)
        ckerns.append(ckern)

    return legs, clegs, lobs, clobs, kerns, ckerns, denoms

usage = """%prog [options]

Generate lobatto.pyx file.
"""
help = {
    'max_order' :
    'maximum order of polynomials [default: %default]',
    'plot' :
    'plot polynomials',
}

def main():
    parser = OptionParser(usage=usage, version='%prog')
    parser.add_option('-m', '--max-order', metavar='order', type=int,
                      action='store', dest='max_order',
                      default=10, help=help['max_order'])
    parser.add_option('', '--plot',
                      action='store_true', dest='plot',
                      default=False, help=help['plot'])
    options, args = parser.parse_args()

    max_order = options.max_order

    legs, clegs, lobs, clobs, kerns, ckerns, denoms = gen_lobbato(max_order)

    if options.plot:
        x, y = sm.symbols('x,y')

        plt.figure(1)
        plt.clf()

        vx = nm.linspace(-1, 1, 100)

        for ii, lob in enumerate(lobs):
            print ii
            print lob
            # print sm.integrate(lob, (x, -1, 1))
            print lob.as_poly(x).all_coeffs()

            vy = [float(lob.subs(x, xx)) for xx in vx]
            plt.plot(vx, vy)

        plt.figure(2)
        plt.clf()

        for ii, kern in enumerate(kerns):
            print ii + 2
            print kern
            print kern.as_poly(x).all_coeffs()

            vy = [float(kern.subs(x, xx)) for xx in vx]
            plt.plot(vx, vy)


        plt.figure(3)
        plt.clf()

        for ii, leg in enumerate(legs):
            print ii
            print sm.simplify(leg)
            print leg.as_poly(y).all_coeffs()
            # print sm.integrate(leg, ('y', -1, 1))

            vy = [float(leg.subs(y, xx)) for xx in vx]
            plt.plot(vx, vy)

        plt.show()

    indir = InDir(os.path.join(top_dir, 'sfepy/fem/extmods/'))
    fd = open(indir('lobatto_template.pyx'), 'r')
    template = fd.read()
    fd.close()

    filename = indir('lobatto.pyx')

    fd = open(filename, 'w')

    out = []

    names_lobatto = []
    out.append('\n# Lobatto functions.\n')
    for ii, clob in enumerate(clobs):
        name = 'lobatto_%03d' % ii
        function = cdef % (name, clob)
        out.append(function)
        names_lobatto.append(name)

    names_kernel = []
    out.append('\n# Kernel functions.\n')
    for ii, ckern in enumerate(ckerns):
        name = 'kernel_%03d' % (ii + 2)
        function = cdef % (name, ckern)
        out.append(function)
        names_kernel.append(name)

    names_legendre = []
    out.append('\n# Legendre functions.\n')
    for ii, cleg in enumerate(clegs):
        name = 'legendre_%03d' % ii
        function = cdef % (name, cleg.replace('y', 'x'))
        out.append(function)
        names_legendre.append(name)

    out.append('\n# Lists of functions.\n')

    args = ', '.join(['&%s' % name for name in names_lobatto])
    list_lobatto = fun_list % ('lobatto', max_order + 1, 'lobatto', args)
    out.append(list_lobatto)

    args = ', '.join(['&%s' % name for name in names_kernel])
    list_kernel = fun_list % ('kernel', max_order - 1, 'kernel', args)
    out.append(list_kernel)

    args = ', '.join(['&%s' % name for name in names_legendre])
    list_legendre = fun_list % ('legendre', max_order + 1, 'legendre', args)
    out.append(list_legendre)

    fd.write(template % ''.join(out))

    fd.close()

if __name__ == '__main__':
    main()
