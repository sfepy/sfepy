"""
This module contains functions that have different names or behavior
depending on NumPy and Scipy versions.
"""
import numpy as nm
import scipy as sc

__all__ = ['in1d', 'unique']

try:
    in1d = nm.in1d

except AttributeError:
    in1d = nm.setmember1d

unique = nm.unique
try:
    nm.unique([0], return_index=True, return_inverse=True)

except TypeError:
    unique = nm.unique1d


try:
    factorial = sc.factorial

except AttributeError:
    import scipy.special as scm

    factorial = scm.factorial

if nm.lib.NumpyVersion(nm.__version__) >= '1.14.0':
    def lstsq(a, b, **kwargs):
        return nm.linalg.lstsq(a, b, rcond=None, **kwargs)
else:
    lstsq = nm.linalg.lstsq


try:
    block = nm.block

except AttributeError:
    from numpy.core import numeric as _nx
    from numpy.core.numeric import array

    def _block_check_depths_match(arrays, parent_index=[]):
        """
        Copied from NumPy 1.14.1.
        """
        def format_index(index):
            idx_str = ''.join('[{}]'.format(i) for i in index if i is not None)
            return 'arrays' + idx_str
        if type(arrays) is tuple:
            # not strictly necessary, but saves us from:
            #  - more than one way to do things - no point treating tuples like
            #    lists
            #  - horribly confusing behaviour that results when tuples are
            #    treated like ndarray
            raise TypeError(
                '{} is a tuple. '
                'Only lists can be used to arrange blocks, and np.block does '
                'not allow implicit conversion from tuple to ndarray.'.format(
                    format_index(parent_index)
                )
            )
        elif type(arrays) is list and len(arrays) > 0:
            idxs_ndims = (_block_check_depths_match(arr, parent_index + [i])
                          for i, arr in enumerate(arrays))

            first_index, max_arr_ndim = next(idxs_ndims)
            for index, ndim in idxs_ndims:
                if ndim > max_arr_ndim:
                    max_arr_ndim = ndim
                if len(index) != len(first_index):
                    raise ValueError(
                        "List depths are mismatched. First element was at depth "
                        "{}, but there is an element at depth {} ({})".format(
                            len(first_index),
                            len(index),
                            format_index(index)
                        )
                    )
            return first_index, max_arr_ndim
        elif type(arrays) is list and len(arrays) == 0:
            # We've 'bottomed out' on an empty list
            return parent_index + [None], 0
        else:
            # We've 'bottomed out' - arrays is either a scalar or an array
            return parent_index, _nx.ndim(arrays)


    def _block(arrays, max_depth, result_ndim):
        """
        Copied from NumPy 1.14.1.
        """
        def atleast_nd(a, ndim):
            # Ensures `a` has at least `ndim` dimensions by prepending
            # ones to `a.shape` as necessary
            return array(a, ndmin=ndim, copy=False, subok=True)

        def block_recursion(arrays, depth=0):
            if depth < max_depth:
                if len(arrays) == 0:
                    raise ValueError('Lists cannot be empty')
                arrs = [block_recursion(arr, depth+1) for arr in arrays]
                return _nx.concatenate(arrs, axis=-(max_depth-depth))
            else:
                # We've 'bottomed out' - arrays is either a scalar or an array
                # type(arrays) is not list
                return atleast_nd(arrays, result_ndim)

        try:
            return block_recursion(arrays)
        finally:
            # recursive closures have a cyclic reference to themselves, which
            # requires gc to collect (gh-10620). To avoid this problem, for
            # performance and PyPy friendliness, we break the cycle:
            block_recursion = None


    def block(arrays):
        """
        Copied from NumPy 1.14.1.
        """
        bottom_index, arr_ndim = _block_check_depths_match(arrays)
        list_ndim = len(bottom_index)
        return _block(arrays, list_ndim, max(arr_ndim, list_ndim))
