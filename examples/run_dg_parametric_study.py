#!python
# -*- coding:utf-8 -*-
""""
Script for running parametric studies for DG from conf files
"""

import sys
from os import path
from os.path import join as pjoin
from copy import copy, deepcopy
import argparse
import numpy as nm
import pandas as pd

from sfepy.applications.pde_solver_app import PDESolverApp

from sfepy.base.conf import ProblemConf
from sfepy.base.ioutils import ensure_path
from sfepy.base.base import (get_default, output, assert_,
                             Struct, basestr, IndexedStruct)
from sfepy.discrete import Integral, Integrals, Material
from sfepy.discrete.common.mappings import get_jacobian


from examples.dg.example_dg_common import get_gen_block_mesh_hook
from sfepy.discrete.dg.my_utils.plot_1D_dg import clear_folder


parser = argparse.ArgumentParser(description='Run SfePy DG example conf python files',
                                 epilog='(c) 2019 by T. Zitka , Man-machine Interaction at NTC UWB')
parser.add_argument("conf_file", help="""File containing problem configuration""")
parser.add_argument('-p', '--plot', help="To plot 1D case", action="store_true")
parser.add_argument('-o', '--output', help="""Root output folder""", default="output")


def iter_h_tens(base_pc, n_refirements, dimensions=(1., 1.), center=(.5, .5)):

    n_don = 2
    for i in range(n_refirements):
        pc = copy(base_pc)

        pc.h = nm.sqrt(2) / (n_don - 1)
        pc.h_coef = 1 / (n_don - 1)
        pc.filename_mesh = get_gen_block_mesh_hook(dimensions, (n_don, n_don), center)

        pc.example_name += "_h" + str(n_don - 1)
        pc.output_folder = pjoin(pc.output_folder, "h" + str(n_don - 1) + "/")

        n_don = n_don + n_don - 1

        yield pc


def iter_order(base_pc, approx_orders=range(1, 7)):
    for approx_order in approx_orders:
        pc = copy(base_pc)
        pc.fields = deepcopy(pc.fields)
        pc.integrals = deepcopy(pc.integrals)

        update_approx_order(pc, approx_order)

        pc.example_name += "_o" + str(pc.approx_order)
        pc.output_folder = pjoin(pc.output_folder, "o" + str(pc.approx_order) + "/")

        yield pc


def update_approx_order(pc, approx_order):
    pc.approx_order = approx_order
    for val in pc.integrals.values():
        val.order = 2 * approx_order
    for val in pc.fields.values():
        val.approx_order = str(approx_order) + 'd'
    return pc


def flatten(it, map_iter='values', max_depth=128):
    from collections import Iterable, Mapping
    from operator import methodcaller
    if max_depth < 0:
        raise RecursionError('maximum recursion depth exceeded in flatten')
    elif isinstance(it, str):
        yield it
    elif isinstance(it, Mapping):
        for item in methodcaller(map_iter)(it):
            yield from flatten(item, map_iter=map_iter, max_depth=max_depth-1)
    elif isinstance(it, Iterable):
        for item in it:
            yield from flatten(item, map_iter=map_iter, max_depth=max_depth-1)
    else:
        yield it


def iterate_all(iterators, base):
    result = [base]
    for iterator in iterators:
        result = [list(iterator(init)) for init in result]
        result = list(flatten(result))

    for prod in result:
        yield prod


def main_from_file(argv):
    if argv is None:
        argv = sys.argv[1:]

    args = parser.parse_args(argv)
    conf_file_name = args.conf_file
    output_folder = args.output

    output("Using conf file {} as base".format(conf_file_name))
    pc_base = ProblemConf.from_file(conf_file_name)
    pc_base.output_folder = pjoin(output_folder, pc_base.example_name)

    err_list = []

    for pc in iterate_all([iter_order, iter_h_tens], pc_base):
        pc.output_name_trunk = pjoin(pc.output_folder, pc.example_name)
        output("Running {}".format(pc.example_name))
        ensure_path(pc.output_name_trunk)
        output_format = "{}.*.{}".format(pc.output_name_trunk,
                                          pc.options.output_format
                                          if hasattr(pc.options, "output_format") else "vtk")
        output("Output set to {}, clearing ...".format(output_format))
        clear_folder(output_format, confirm=False)

        sa = PDESolverApp(pc, Struct(output_filename_trunk=pc.output_name_trunk,
                                     save_ebc=False,
                                     save_ebc_nodes=False,
                                     save_region=False,
                                     save_regions=False,
                                     save_regions_as_groups=False,
                                     save_field_meshes=False,
                                     solve_not=False), "sfepy")
        pb, state = sa()

        ts = pb.ts
        ts.set_step(step=ts.n_step - 1)

        idiff = Integral('idiff', max(pc.approx_order, 10))
        field = pb.fields['density']
        u = pb.equations.variables["u"]
        num_qp = pb.evaluate('ev_volume_integrate.idiff.Omega(u)',
                             u=u, ts=pb.ts,
                             integrals=Integrals([idiff]), mode='qp')
        aux = Material('aux', function=pc.sol_fun)

        ana_qp = pb.evaluate('ev_volume_integrate_mat.idiff.Omega(aux.u, u)',
                             aux=aux, u=u, ts=pb.ts,
                             integrals=Integrals([idiff]), mode='qp')

        det = get_jacobian(field, idiff)

        diff_l2 = nm.sqrt((((num_qp - ana_qp) ** 2) * det).sum())
        ana_l2 = nm.sqrt((((ana_qp) ** 2) * det).sum())
        error = diff_l2 / ana_l2

        avrg_vol = nm.mean(pb.fields["density"].domain.cmesh.get_volumes(pb.get_dim()))
        var_vol = nm.var(pb.fields["density"].domain.cmesh.get_volumes(pb.get_dim()))
        err_list.append({"h": pc.h, "h_coef": pc.h_coef, "order": pc.approx_order,
                         "avrg_vol": avrg_vol, "var_vol": var_vol,
                         "ana_l2": ana_l2, "err_rel": error, "err_l2": diff_l2})

    err_df = pd.DataFrame(err_list, columns=["h", "h_coef", "order", "avrg_vol", "var_vol", "ana_l2", "err_rel", "err_l2"])
    err_df.to_csv(pjoin(output_folder, pc_base.example_name, "h-order_error.csv"))


if __name__ == '__main__':
    main_from_file(sys.argv[1:])
