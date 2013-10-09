"""
Helper functions related to mesh facets.
"""
import numpy as nm

from sfepy.base.base import dict_to_array, assert_
from sfepy.base.compat import in1d
from sfepy.linalg import permutations, map_permutations

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

def get_dof_orientation_maps(n_fp, raw_ori_maps):
    """
    Given description of facet DOF nodes, return the corresponding
    integer coordinates and orientation maps.

    Notes
    -----
    Assumes single facet type in all groups.
    """
    if n_fp <= 3: # Simplex facet.
        ori_maps = raw_ori_maps

    else: # Tensor product facet.
        ori_maps = {}
        for ig, ori_map in raw_ori_maps.iteritems():
            ori_maps[ig] = {}
            for key in ori_map.iterkeys():
                new_key = _quad_ori_groups[key]
                ori_maps[ig][key] = [None, _quad_orientations[new_key]]

    return ori_maps

def permute_facets(facets, ori, ori_map):
    """
    Return a copy of `facets` array with vertices sorted lexicographically.
    """
    assert_((in1d(nm.unique(ori), ori_map.keys())).all())

    permuted_facets = facets.copy()

    for key, ori_map in ori_map.iteritems():
        perm = ori_map[1]
        ip = nm.where(ori == key)[0]
        for ic0, ic1 in enumerate(perm):
            permuted_facets[ip,ic0] = facets[ip,ic1]

    return permuted_facets

def get_facet_dof_permutations(n_fp, igs, nodes):
    """
    Prepare DOF permutation vector for each possible facet orientation.
    """
    ori_map, _, _ = build_orientation_map(n_fp)
    raw_ori_maps = {}
    for ig in igs:
        raw_ori_maps[ig] = ori_map
    ori_maps = get_dof_orientation_maps(n_fp, raw_ori_maps)

    int_coors = nodes
    ori = nm.empty((int_coors.shape[0],), dtype=nm.int32)
    dof_perms = {}
    for ig, ori_map in ori_maps.iteritems():
        dof_perms[ig] = {}
        for key in ori_map.iterkeys():
            ori.fill(key)
            permuted_int_coors = permute_facets(int_coors, ori, ori_map)

            perm = map_permutations(int_coors, permuted_int_coors,
                                    check_same_items=True)

            dof_perms[ig][key] = perm

    adof_perms = {}
    for ig, val in dof_perms.iteritems():
        adof_perms[ig] = dict_to_array(val)

    return adof_perms
