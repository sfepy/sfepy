"""
Helper functions related to mesh facets.
"""
import numpy as nm

from sfepy.linalg import permutations

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

_quad_orientations = {
    0  : [0, 1, 3, 2],
    7  : [3, 2, 0, 1],
    11 : [2, 3, 0, 1],
    30 : [1, 0, 2, 3],
    33 : [1, 0, 3, 2],
    52 : [2, 3, 1, 0],
    56 : [3, 2, 1, 0],
    63 : [0, 1, 2, 3],
}

def build_orientation_map(n_fp):
    """
    The keys are binary masks of the lexicographical ordering of facet
    vertices. A bit i set to one means `v[i] < v[i+1]`.

    The values are `[original_order, permutation]`, where `permutation` can be
    used to sort facet vertices lexicographically. Hence `permuted_facet =
    facet[permutation]`.
    """
    indices = range(n_fp)

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
