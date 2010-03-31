"""
Functions to compute some tensor-related quantities usual in continuum mechanics.
"""
from sfepy.base.base import *
from sfepy.base.la import dot_sequences, make_axis_rotation_matrix

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

def prepare_cylindrical_transform(coors, origin):
    """
    Prepare matrices for transforming tensors into cylindrical coordinates with
    the axis 'z' in a given origin.

    Parameters
    ----------
    coors : array
        The Cartesian coordinates.
    origin : array of length 3
        The origin.

    Returns
    -------
    mtx : array
        The array of transformation matrices for each coordinate in `coors`.
    """
    x, y = coors[:,0] - origin[0], coors[:,1] - origin[1]
    theta = nm.arctan2(y, x)

    mtx = nm.zeros((coors.shape[0], 3, 3), dtype=nm.float64)
    for ii, th in enumerate(theta):
        mtx[ii] = make_axis_rotation_matrix([0.0, 0.0, 1.0], th)

    return mtx

def transform_data(data, coors=None, mode='cylindrical', mtx=None):
    """
    Transform vector or tensor data components into another coordinate system.

    For vectors:

    .. math::
        \ul{v}' = M \cdot \ul{v}

    For tensors (assuming orthogonal coordinates):

    .. math::
        \ull{t}' = M^{-T} \cdot \ull{t} \cdot M

    Parameters
    ----------
    data : array
        The vectors or tensors (symmetric storage) to be transformed.
    coors : array
        The Cartesian coordinates of the data. Not needed when `mtx` argument
        is given.
    mode : one of ['cylindrical']
        The requested coordinate system. Not needed when `mtx` argument
        is given.
    mtx : array
        The array of transformation matrices :math:`M` for each data row.

    Returns
    -------
    new_data : array
        The transformed data.
    """
    if (coors is None) and (mtx is None):
        raise ValueError('one of (coors, mtx) arguments must be set!')
    
    if mtx is None:
        if mode == 'cylindrical':
            mtx = prepare_cylindrical_transform(coors, [0.0, 0.0, 0.0])

        else:
            raise ValueError('transformation mode %s is not supported!' % mode)

    shape = data.shape

    data = nm.squeeze(data)
    if data.ndim == 1:
        data.shape = (1, data.shape[0])

    if data.shape[1] == 3: # Vectors.
        new_data = dot_sequences(mtx, data)

    elif data.shape[1] == 6: # Symmetric tensors.
        ii = [[0, 3, 4], [3, 1, 5], [4, 5, 2]]

        aux = data[:,ii]
        aux2 = dot_sequences(dot_sequences(mtx, aux), mtx.transpose((0, 2, 1)))
        assert nm.allclose(aux2[0], nm.dot(nm.dot(mtx[0], aux[0]), mtx[0].T))

        aux3 = aux2.reshape((aux2.shape[0], 9))

        new_data = aux3[:, [0, 4, 8, 1, 2, 5]]

    else:
        raise ValueError('unsupported data shape! (%s)' % str(data.shape))

    # Restore the correct shape.
    new_data.shape = shape

    return new_data
