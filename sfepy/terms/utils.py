from __future__ import print_function
from __future__ import absolute_import
import numpy as nm
from six.moves import range

def check_finiteness(data, info):
    is_finite = nm.isfinite(data)
    if not is_finite.all():
        ii = nm.where(is_finite == False)
        print(ii)
        print(data[ii])
        msg = 'infinite %s!, see above' % info
        raise ValueError(msg)

def get_range_indices(num):
    """
    Return indices and slices in given range.

    Returns
    -------
    indx : list of tuples
        The list of `(ii, slice(ii, ii + 1))` of the indices. The first
        item is the index itself, the second item is a convenience slice
        to index components of material parameters.
    """
    indx = [(ii, slice(ii, ii + 1)) for ii in range(num)]

    return indx
