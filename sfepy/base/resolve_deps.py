"""
Functions for resolving dependencies.
"""
from __future__ import absolute_import
import itertools as it

from sfepy.base.base import basestr
import six
from six.moves import range

def get_nums(deps):
    """
    Get number of prerequisite names for each name in dependencies.
    """
    nums = {}
    for key, val in six.iteritems(deps):
        nums[key] = len(val)

    return nums

def solvable(deps, names):
    """
    Return True if `names` form a solvable block, i.e. the set of names equals
    to the set of their prerequisites.
    """
    if not names: return False # Missing self-reference.

    dep_names = set()
    for name in names:
        dep_names.update(deps[name])

    return dep_names == set(names)

def remove_known(deps, known):
    """
    Remove known names from dependencies.
    """
    if isinstance(known, basestr):
        out = {}
        for key, val in six.iteritems(deps):
            if key == known: continue
            out[key] = [ii for ii in val if ii != known]

        return out

    else:
        out = deps
        for ii in known:
            out = remove_known(out, ii)

        return out

def try_block(deps, num):
    """
    Return generator of lists of solvable blocks of the length `num`.
    """
    keys = list(deps.keys())
    for ic in it.combinations(keys, num):
        if solvable(deps, ic):
            yield sorted(ic)

def resolve(deps):
    """
    Resolve dependencies among equations so that smaller blocks are solved
    first.

    The dependencies are given in terms of variable names.

    Parameters
    ----------
    deps : dict
        The dependencies as a dictionary with names as keys and sets of
        prerequisite names as values.

    Returns
    -------
    order : list
        The list of blocks in the order of solving. Each block is a list of
        names.
    """
    order = []
    if not(len(deps)): return order

    nums = get_nums(deps)
    ib0 = min(nums.values())
    for ib in range(ib0, len(deps) + 1):
        blocks = [ii for ii in try_block(deps, ib)]
        if len(blocks):
            new_deps = remove_known(deps, blocks[0])
            order.extend([blocks[0]] + resolve(new_deps))

        if len(list(it.chain(*order))) == len(deps):
            break

    return order
