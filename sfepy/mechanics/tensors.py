"""
Functions to compute some tensor-related quantities usual in continuum mechanics.
"""
from sfepy.base.base import *

def dim2sym(dim):
    """
    Given the space dimension, return the symmetric storage size.
    """
    return (dim + 1) * dim / 2

def sym2dim(sym):
    """
    Given the symmetric storage size, return the space dimension.
    """
    return int((sym / 3) + 1)

def get_trace(tensor, sym_storage=True):
    """
    The trace of a tensor.
    """
    if sym_storage:
        dim = sym2dim(tensor.shape[1])
        trace = nm.sum(tensor[:,:dim], axis=1)

    else:
        trace = nm.trace(tensor, axis1=1, axis2=2)

    return trace

def get_volumetric_tensor(tensor, sym_storage=True):
    """
    The volumetric part of a tensor.
    """
    dim = tensor.shape[1]
    if sym_storage:
        dim = sym2dim(dim)
    
    trace = get_trace(tensor, sym_storage=sym_storage)
    val = trace / float(dim)

    if sym_storage:
        vt = nm.zeros_like(tensor)
        vt[:,:dim] = val[:,None]

    else:
        vt = val[:,None,None] * nm.eye(dim, dtype=nm.float64)

    return vt

def get_deviator(tensor, sym_storage=True):
    """
    The deviatoric part (deviator) of a tensor.
    """
    vt = get_volumetric_tensor(tensor, sym_storage=sym_storage)
    dev = tensor - vt

    return dev

def get_von_mises_stress(stress, sym_storage=True):
    r"""
    Given a symmetric stress tensor, compute the von Mises stress (also known
    as Equivalent tensile stress).

    Notes
    -----
    .. math::
        \sigma_V = \sqrt{\frac{(\sigma_{11} - \sigma_{22})^2 +
        (\sigma_{22} - \sigma_{33})^2 + (\sigma_{11} - \sigma_{33})^2 + 6
        (\sigma_{12}^2 + \sigma_{13}^2 + \sigma_{23}^2)}{2}}
    """
    dim = stress.shape[1]
    if sym_storage:
        dim = sym2dim(dim)

    assert_(dim == 3)

    if sym_storage:
        s11 = stress[:,0]
        s22 = stress[:,1]
        s33 = stress[:,2]
        s12 = stress[:,3]
        s13 = stress[:,4]
        s23 = stress[:,5]

    else:
        s11 = stress[:,0,0]
        s22 = stress[:,1,1]
        s33 = stress[:,2,2]
        s12 = stress[:,0,1]
        s13 = stress[:,0,2]
        s23 = stress[:,1,2]

    vms = nm.sqrt(0.5 * ((s11 - s22)**2 + (s22 - s33)**2 + (s11 - s33)**2
                         + 6.0 * (s12**2 + s13**2 + s23**2)))[:,None]

    return vms
