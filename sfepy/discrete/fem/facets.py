"""
Helper functions related to mesh facets and Lagrange FE approximation.

Line: ori - iter:

0 - iter0
1 - iter1

Triangle: ori - iter:

0 - iter21
1 - iter12
3 - iter02
4 - iter20
6 - iter10
7 - iter01

Possible couples:

1, 4, 7 <-> 0, 3, 6

Square: ori - iter:

 0 - iter10x01y
 7 - iter10y01x
11 - iter01y01x
30 - iter01x10y
33 - iter10x10y
52 - iter01y10x
56 - iter10y10x
63 - iter01x01y

Possible couples:

7, 33, 52, 63 <-> 0, 11, 30, 56

_quad_ori_groups:

i < j < k < l

all faces are permuted to

l --- k
|     |
|     |
i --- j

ijkl

which is the same as

l --- j
|     |
|     |
i --- k

ikjl

k --- l
|     |
|     |
i --- j

ijlk

- start at one vertex and go around clock-wise or anticlock-wise
-> 8 groups of 3
-> same face nodes order in
ijkl (63), ikjl (59), ijlk (31)
ilkj (11), iklj (15), iljk (43)
jkli ( 7), jlki ( 3), kjli ( 6)
kjil (56), jkil (57), ljik (48)
lijk (52), likj (20), kijl (60)
lkji ( 0), ljki ( 4), klji ( 1)
klij (33), lkij (32), jlik (41)
jilk (30), kilj (22), jikl (62)
"""
from __future__ import print_function
from __future__ import absolute_import
import numpy as nm
import six
from six.moves import range

_quad_ori_groups = {
    0 : 0,
    1 : 0,
    3 : 7,
    4 : 0,
    6 : 7,
    7 : 7,
    11 : 11,
    15 : 11,
    20 : 52,
    22 : 30,
    30 : 30,
    31 : 63,
    32 : 33,
    33 : 33,
    41 : 33,
    43 : 11,
    48 : 56,
    52 : 52,
    56 : 56,
    57 : 56,
    59 : 63,
    60 : 52,
    62 : 30,
    63 : 63,
}

def build_orientation_map(n_fp):
    """
    The keys are binary masks of the lexicographical ordering of facet
    vertices. A bit i set to one means `v[i] < v[i+1]`.

    The values are `[original_order, permutation]`, where `permutation` can be
    used to sort facet vertices lexicographically. Hence `permuted_facet =
    facet[permutation]`.
    """
    from sfepy.linalg import permutations

    indices = list(range(n_fp))

    cmps = [(i1, i2) for i2 in indices for i1 in indices[:i2]]
    powers = [2**ii for ii in range(len(cmps))]

    ori_map = {}
    for indx in permutations(indices):
        key = 0
        sign = 1
        for ip, power in enumerate(powers):
            i1, i2 = cmps[ip]
            less = (indx[i1] < indx[i2])
            key += power * less
            if not less:
                sign *= -1

        isort = nm.argsort(indx)
        ori_map[key] = [indx, isort]

    return ori_map, cmps, powers

def iter0(num):
    for ir in range(num - 1, -1, -1):
        yield ir

def iter1(num):
    for ir in range(num):
        yield ir

ori_line_to_iter = {
    0 : iter0,
    1 : iter1,
}

def make_line_matrix(order):
    if (order < 2):
        return nm.zeros((0, 0), dtype=nm.int32)

    oo = order - 1
    mtx = nm.arange(oo, dtype=nm.int32)

    return mtx

def iter01(num):
    for ir in range(num - 1, -1, -1):
        for ic in range(ir + 1):
            yield ir, ic

def iter10(num):
    for ir in range(num - 1, -1, -1):
        for ic in range(ir, -1, -1):
            yield ir, ic

def iter02(num):
    for ic in range(num):
        for ir in range(num - 1, ic - 1, -1):
            yield ir, ic

def iter20(num):
    for ic in range(num):
        for ir in range(ic, num):
            yield ir, ic

def iter12(num):
    for idiag in range(num):
        irs, ics = nm.diag_indices(num - idiag)
        for ii in range(irs.shape[0] - 1, -1, -1):
            yield irs[ii] + idiag, ics[ii]

def iter21(num):
    for idiag in range(num):
        irs, ics = nm.diag_indices(num - idiag)
        for ii in range(irs.shape[0]):
            yield irs[ii] + idiag, ics[ii]

ori_triangle_to_iter = {
    0 : iter21,
    1 : iter12,
    3 : iter02,
    4 : iter20,
    6 : iter10,
    7 : iter01,
}

def make_triangle_matrix(order):
    if (order < 3):
        return nm.zeros((0, 0), dtype=nm.int32)

    oo = order - 2
    mtx = nm.zeros((oo, oo), dtype=nm.int32)
    for ii, (ir, ic) in enumerate(iter01(oo)):
        mtx[ir, ic] = ii

    return mtx

def iter01x01y(num):
    for ir in range(num):
        for ic in range(num):
            yield ir, ic

def iter01y01x(num):
    for ir, ic in iter01x01y(num):
        yield ic, ir

def iter10x01y(num):
    for ir in range(num - 1, -1, -1):
        for ic in range(num):
            yield ir, ic

def iter10y01x(num):
    for ir, ic in iter10x01y(num):
        yield ic, ir

def iter01x10y(num):
    for ir in range(num):
        for ic in range(num - 1, -1, -1):
            yield ir, ic

def iter01y10x(num):
    for ir, ic in iter01x10y(num):
        yield ic, ir

def iter10x10y(num):
    for ir in range(num - 1, -1, -1):
        for ic in range(num - 1, -1, -1):
            yield ir, ic

def iter10y10x(num):
    for ir, ic in iter10x10y(num):
        yield ic, ir

ori_square_to_iter = {
    0 : iter10x01y,
    7 : iter10y01x,
    11 : iter01y01x,
    30 : iter01x10y,
    33 : iter10x10y,
    52 : iter01y10x,
    56 : iter10y10x,
    63 : iter01x01y,
}

def make_square_matrix(order):
    if (order < 2):
        return nm.zeros((0, 0), dtype=nm.int32)

    oo = order - 1
    mtx = nm.arange(oo * oo, dtype=nm.int32)
    mtx.shape = (oo, oo)

    return mtx

def get_facet_dof_permutations(n_fp, order):
    """
    Prepare DOF permutation vector for each possible facet orientation.
    """
    from sfepy.base.base import dict_to_array

    if n_fp == 2:
        mtx = make_line_matrix(order)
        ori_map = ori_line_to_iter
        fo = order - 1

    elif n_fp == 3:
        mtx = make_triangle_matrix(order)
        ori_map = ori_triangle_to_iter
        fo = order - 2

    elif n_fp == 4:
        mtx = make_square_matrix(order)
        ori_map = {}
        for key, val in six.iteritems(_quad_ori_groups):
            ori_map[key] = ori_square_to_iter[val]
        fo = order - 1

    else:
        raise ValueError('unsupported number of facet points! (%d)' % n_fp)

    dof_perms = {}
    for key, itfun in six.iteritems(ori_map):
        dof_perms[key] = [mtx[ii] for ii in itfun(fo)]

    dof_perms = dict_to_array(dof_perms)

    return dof_perms

if __name__ == '__main__':
    order = 5
    mtx = make_triangle_matrix(order)
    print(mtx)

    oo = order - 2
    print([mtx[ir, ic] for ir, ic in ori_triangle_to_iter[0](oo)])
    print([mtx[ir, ic] for ir, ic in ori_triangle_to_iter[1](oo)])
    print([mtx[ir, ic] for ir, ic in ori_triangle_to_iter[3](oo)])
    print([mtx[ir, ic] for ir, ic in ori_triangle_to_iter[4](oo)])
    print([mtx[ir, ic] for ir, ic in ori_triangle_to_iter[6](oo)])
    print([mtx[ir, ic] for ir, ic in ori_triangle_to_iter[7](oo)])

    order = 4
    mtx = make_square_matrix(order)
    print(mtx)

    oo = order - 1
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[0](oo)])
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[7](oo)])
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[11](oo)])
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[30](oo)])
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[33](oo)])
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[52](oo)])
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[56](oo)])
    print([mtx[ir, ic] for ir, ic in ori_square_to_iter[63](oo)])
