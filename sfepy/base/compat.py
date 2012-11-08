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
    import scipy.misc as scm

    factorial = scm.factorial
