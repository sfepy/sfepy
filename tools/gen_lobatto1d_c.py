#!/usr/bin/env python
"""
Generate lobatto1d.c and lobatto1h.c files.
"""
import sys
sys.path.append('.')
import os
from argparse import ArgumentParser

import sympy as sm
import numpy as nm
import matplotlib.pyplot as plt

from sfepy import top_dir
from sfepy.base.ioutils import InDir

hdef = 'float64 %s(float64 x);\n'

cdef = """
float64 %s(float64 x)
{
  return(%s);
}
"""

fun_list = """
const fun %s[%d] = {%s};
"""

def gen_lobatto(max_order):
    assert max_order > 2

    x = sm.symbols('x')

    lobs = [0, 1]
    lobs[0] = (1 - x) / 2
    lobs[1] = (1 + x) / 2

    dlobs = [lob.diff('x') for lob in lobs]

    legs = [sm.legendre(0, 'y')]
    clegs = [sm.ccode(legs[0])]
    dlegs = [sm.legendre(0, 'y').diff('y')]
    cdlegs = [sm.ccode(dlegs[0])]

    clobs = [sm.ccode(lob) for lob in lobs]
    cdlobs = [sm.ccode(dlob) for dlob in dlobs]

    denoms = [] # for lobs.

    for ii in range(2, max_order + 1):
        coef = sm.sympify('sqrt(2 * (2 * %s - 1)) / 2' % ii)
        leg = sm.legendre(ii - 1, 'y')

        pleg = leg.as_poly()
        coefs = pleg.all_coeffs()
        denom = max(sm.denom(val) for val in coefs)

        cleg = sm.ccode(sm.horner(leg*denom)/denom)

        dleg = leg.diff('y')
        cdleg = sm.ccode(sm.horner(dleg*denom)/denom)

        lob = sm.simplify(coef * sm.integrate(leg, ('y', -1, x)))
        lobnc = sm.simplify(sm.integrate(leg, ('y', -1, x)))

        plobnc = lobnc.as_poly()
        coefs = plobnc.all_coeffs()
        denom = sm.denom(coef) * max(sm.denom(val) for val in coefs)

        clob = sm.ccode(sm.horner(lob*denom)/denom)

        dlob = lob.diff('x')
        cdlob = sm.ccode(sm.horner(dlob*denom)/denom)

        legs.append(leg)
        clegs.append(cleg)
        dlegs.append(dleg)
        cdlegs.append(cdleg)
        lobs.append(lob)
        clobs.append(clob)
        dlobs.append(dlob)
        cdlobs.append(cdlob)
        denoms.append(denom)

    coef = sm.sympify('sqrt(2 * (2 * %s - 1)) / 2' % (max_order + 1))
    leg = sm.legendre(max_order, 'y')

    pleg = leg.as_poly()
    coefs = pleg.all_coeffs()
    denom = max(sm.denom(val) for val in coefs)

    cleg = sm.ccode(sm.horner(leg*denom)/denom)

    dleg = leg.diff('y')
    cdleg = sm.ccode(sm.horner(dleg*denom)/denom)

    legs.append(leg)
    clegs.append(cleg)
    dlegs.append(dleg)
    cdlegs.append(cdleg)

    kerns = []
    ckerns = []
    dkerns = []
    cdkerns = []
    for ii, lob in enumerate(lobs[2:]):
        kern = sm.simplify(lob / (lobs[0] * lobs[1]))
        dkern = kern.diff('x')

        denom = denoms[ii] / 4
        ckern = sm.ccode(sm.horner(kern*denom)/denom)
        cdkern = sm.ccode(sm.horner(dkern*denom)/denom)

        kerns.append(kern)
        ckerns.append(ckern)
        dkerns.append(dkern)
        cdkerns.append(cdkern)

    return (legs, clegs, dlegs, cdlegs,
            lobs, clobs, dlobs, cdlobs,
            kerns, ckerns, dkerns, cdkerns,
            denoms)

def plot_polys(fig, polys, var_name='x'):
    plt.figure(fig)
    plt.clf()

    x = sm.symbols(var_name)
    vx = nm.linspace(-1, 1, 100)

    for ii, poly in enumerate(polys):
        print(ii)
        print(poly)
        print(poly.as_poly(x).all_coeffs())

        vy = [float(poly.subs(x, xx)) for xx in vx]
        plt.plot(vx, vy)

def append_declarations(out, cpolys, comment, cvar_name, shift=0):
    names = []
    out.append('\n// %s functions.\n' % comment)
    for ii, cpoly in enumerate(cpolys):
        name = '%s_%03d' % (cvar_name, ii + shift)
        function = hdef % name
        out.append(function)
        names.append(name)

    return names

def append_polys(out, cpolys, comment, cvar_name, var_name='x', shift=0):
    names = []
    out.append('\n// %s functions.\n' % comment)
    for ii, cpoly in enumerate(cpolys):
        name = '%s_%03d' % (cvar_name, ii + shift)
        function = cdef % (name, cpoly.replace(var_name, 'x'))
        out.append(function)
        names.append(name)

    return names

def append_lists(out, names, length):
    args = ', '.join(['&%s' % name for name in names])
    name = names[0][:-4]
    _list = fun_list % (name, length, args)
    out.append(_list)

helps = {
    'max_order' :
    'maximum order of polynomials [default: %(default)s]',
    'plot' :
    'plot polynomials',
}

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('-m', '--max-order', metavar='order', type=int,
                        action='store', dest='max_order',
                        default=10, help=helps['max_order'])
    parser.add_argument('--plot',
                        action='store_true', dest='plot',
                        default=False, help=helps['plot'])
    options = parser.parse_args()

    max_order = options.max_order

    (legs, clegs, dlegs, cdlegs,
     lobs, clobs, dlobs, cdlobs,
     kerns, ckerns, dkerns, cdkerns,
     denoms) = gen_lobatto(max_order)

    if options.plot:
        plot_polys(1, lobs)
        plot_polys(11, dlobs)

        plot_polys(2, kerns)
        plot_polys(21, dkerns)

        plot_polys(3, legs, var_name='y')
        plot_polys(31, dlegs, var_name='y')

        plt.show()

    indir = InDir(os.path.join(top_dir, 'sfepy/discrete/fem/extmods/'))

    fd = open(indir('lobatto1d_template.h'), 'r')
    template = fd.read()
    fd.close

    fd = open(indir('lobatto1d.h'), 'w')

    out = []

    append_declarations(out, clobs, 'Lobatto', 'lobatto')
    append_declarations(out, cdlobs, 'Derivatives of Lobatto', 'd_lobatto')

    append_declarations(out, ckerns, 'Kernel', 'kernel',
                        shift=2)
    append_declarations(out, cdkerns, 'Derivatives of kernel', 'd_kernel',
                        shift=2)

    append_declarations(out, clegs, 'Legendre', 'legendre')
    append_declarations(out, cdlegs, 'Derivatives of Legendre', 'd_legendre')

    fd.write(template.replace('// REPLACE_TEXT', ''.join(out)))

    fd.close()


    fd = open(indir('lobatto1d_template.c'), 'r')
    template = fd.read()
    fd.close()

    fd = open(indir('lobatto1d.c'), 'w')

    out = []

    names_lobatto = append_polys(out, clobs,
                                 'Lobatto', 'lobatto')
    names_d_lobatto = append_polys(out, cdlobs,
                                   'Derivatives of Lobatto', 'd_lobatto')

    names_kernel = append_polys(out, ckerns,
                                'Kernel', 'kernel',
                                shift=2)
    names_d_kernel = append_polys(out, cdkerns,
                                  'Derivatives of kernel', 'd_kernel',
                                  shift=2)

    names_legendre = append_polys(out, clegs,
                                  'Legendre', 'legendre',
                                  var_name='y')
    names_d_legendre = append_polys(out, cdlegs,
                                    'Derivatives of Legendre', 'd_legendre',
                                    var_name='y')

    out.append('\n// Lists of functions.\n')

    out.append('\nconst int32 max_order = %d;\n' % max_order)

    append_lists(out, names_lobatto, max_order + 1)
    append_lists(out, names_d_lobatto, max_order + 1)

    append_lists(out, names_kernel, max_order - 1)
    append_lists(out, names_d_kernel, max_order - 1)

    append_lists(out, names_legendre, max_order + 1)
    append_lists(out, names_d_legendre, max_order + 1)

    fd.write(template.replace('// REPLACE_TEXT', ''.join(out)))

    fd.close()

if __name__ == '__main__':
    main()
