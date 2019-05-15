#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser
import numpy as nm

import os
import sys

sys.path.append('.')

from sfepy.base.base import output, Struct, IndexedStruct
from sfepy.discrete import (FieldVariable, Material, Integral, Function,
                            Equation, Equations, Problem)
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.terms import Term
from sfepy.discrete.conditions import Conditions, EssentialBC, InitialCondition
from sfepy.solvers import Solver

helps = {
    'tss_name'     : 'time stepping solver name [default: %(default)s]',
    'implicit'     : 'use implicit time stepping',
    'order'        : 'field approximation order [default: %(default)s]',
    'refine'       : 'uniform mesh refinement level [default: %(default)s]',
    'mesh_filename':
        'mesh file name [default: %(default)s]',
    'output_dir'   :
        'output directory [default: %(default)s]',
}


def main():
    from sfepy import data_dir

    parser = ArgumentParser()
    parser.add_argument('--version', action='version', version='%(prog)s')
    parser.add_argument('--tss',
                        action="store", dest='tss_name',
                        default='simple', help=helps['tss_name'])
    parser.add_argument('-i', '--implicit',
                        action="store_true", dest='implicit',
                        default=False, help=helps['implicit'])
    parser.add_argument('--order', metavar='int', type=int,
                        action='store', dest='order',
                        default=2, help=helps['order'])
    parser.add_argument('-r', '--refine', metavar='int', type=int,
                        action='store', dest='refine',
                        default=0, help=helps['refine'])
    parser.add_argument('mesh_filename', nargs='?',
                        default=data_dir + '/meshes/2d/square_quad.mesh',
                        help=helps['mesh_filename'])
    parser.add_argument('output_dir', nargs='?',
                        default='output/fem_adv_1D',
                        help=helps['output_dir'])
    options = parser.parse_args()

    mesh = Mesh.from_file(options.mesh_filename)
    domain = FEDomain('domain', mesh)

    if options.refine > 0:
        for ii in range(options.refine):
            output('refine %d...' % ii)
            domain = domain.refine()
            output('... %d nodes %d elements'
                   % (domain.shape.n_nod, domain.shape.n_el))

    omega = domain.create_region('Omega', 'all')
    gamma = domain.create_region('Gamma', 'vertices of surface', 'facet')

    field = Field.from_args('fu', nm.float64, 1, omega,
                            approx_order=options.order)

    u = FieldVariable('u', 'unknown', field, history=1)
    v = FieldVariable('v', 'test', field, primary_var_name='u')

    m = Material('m', vel=[[1], [0]])

    integral = Integral('i', order=2 * options.order)
    # integral = Integral('i', order=3)

    t1 = Term.new('dw_volume_dot(v, du/dt)',
                  integral, omega, v=v, u=u)
    if options.implicit:
        t2 = Term.new('dw_advect_div_free(m.vel, v, u)',
                      integral, omega, m=m, v=v, u=u)

    else:
        t2 = Term.new('dw_advect_div_free(m.vel, v, u[-1])',
                      integral, omega, m=m, v=v, u=u)

    eq = Equation('eq', t1 + t2)
    eqs = Equations([eq])

    # Initial conditions.
    def get_ic(coors, ic):
        x, y = coors.T
        x = (x - x.min()) / (x.max() - x.min())
        val = .3 * nm.piecewise(
                x, [x <= 0.1, x >= 0.1, .3 < x],
                [0, lambda x: nm.exp(1 / ((10 * (x - .2)) ** 2 - 1)) / nm.exp(1 / (- 1)), 0]
        )
        return val

    ic_fun = Function('ic_fun', get_ic)
    ic = InitialCondition('ic', omega, {'u.all': ic_fun})

    fix = EssentialBC('fix', gamma, {'u.all': 0.0})

    ls = Solver.any_from_conf(Struct(**{
        'name': 'ls', 'kind': 'ls.scipy_direct',
    }))
    nls_status = IndexedStruct()
    nls = Solver.any_from_conf(Struct(**{
        'name'     : 'newton',
        'kind'     : 'nls.newton',
        'i_max'    : 1,
        'is_linear': True,
    }), lin_solver=ls, status=nls_status)
    tss = Solver.any_from_conf(Struct(**{
        'name'  : 'tss',
        'kind'  : 'ts.%s' % options.tss_name,
        't0'    : 0.0,
        't1'    : 1.0,
        'dt'    : 0.01,
        'n_step': None,
    }), nls=nls, verbose=True)

    pb = Problem('advection', equations=eqs, active_only=False)
    pb.set_output_dir(options.output_dir)
    pb.save_regions_as_groups(os.path.join(options.output_dir, 'regions'))

    pb.set_bcs(ebcs=Conditions([fix]))
    pb.set_ics(Conditions([ic]))

    pb.set_solver(tss)

    status = IndexedStruct()
    pb.solve(status=status)
    output(status)


if __name__ == '__main__':
    main()
