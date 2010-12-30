"""
Functions to compute some tensor-related quantities usual in continuum mechanics.
"""
import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import assert_, Struct
from sfepy.linalg \
     import apply_to_sequence, dot_sequences, make_axis_rotation_matrix

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

def get_full_indices(dim):
    """
    The indices for converting the symmetric storage to the full storage.
    """
    return {
        2 : [[0, 2], [2, 1]],
        3 : [[0, 3, 4], [3, 1, 5], [4, 5, 2]],
    }[dim]

def get_sym_indices(dim):
    """
    The indices for converting the full storage to the symmetric storage.
    """
    return {
        2 : [0, 3, 1],
        3 : [0, 4, 8, 1, 2, 5],
    }[dim]

def get_non_diagonal_indices(dim):
    """
    The non_diagonal indices for the full vector storage.
    """
    return {
        2 : ([1], [2]),
        3 : ([1, 2, 5], [3, 6, 7]),
    }[dim]

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
        ii = get_full_indices(3)

        aux = data[:,ii]
        aux2 = dot_sequences(dot_sequences(mtx, aux), mtx.transpose((0, 2, 1)))
        assert nm.allclose(aux2[0], nm.dot(nm.dot(mtx[0], aux[0]), mtx[0].T))

        aux3 = aux2.reshape((aux2.shape[0], 9))

        new_data = aux3[:, get_sym_indices(3)]

    else:
        raise ValueError('unsupported data shape! (%s)' % str(data.shape))

    # Restore the correct shape.
    new_data.shape = shape

    return new_data

class StressTransform(Struct):
    """
    Encapsulates functions to convert various stress tensors in the symmetric
    storage given the deformation state.
    """

    def __init__(self, def_grad, jacobian=None):
        """
        Set :math:`\ull{F} = \pdiff{\ul{x}}{\ul{X}}` and optionally also
        :math:`J = \det(\ull{F})`.
        """
        self.def_grad = def_grad
        self.n_el, self.n_qp, self.dim = self.def_grad.shape[:3]

        self.s2f = get_full_indices(self.dim)
        self.f2s = get_sym_indices(self.dim)

        if jacobian is None:
            self.jacobian = apply_to_sequence(self.def_grad, nla.det,
                                              2, (1, 1))

        else:
            self.jacobian = jacobian

    def _assert_symmetry(self, stress):
        i1, i2 = get_non_diagonal_indices(self.dim)
        assert_(nm.allclose(stress[:,:,i1], stress[:,:,i2]))

    def get_cauchy_from_2pk(self, stress_in):
        """
        Get the Cauchy stress given the second Piola-Kirchhoff stress.
        
        .. math::

            \sigma_{ij} = J^{-1} F_{ik} S_{kl} F_{jl}
        """

        stress_in_full = stress_in[:,:,self.s2f,0]

        val_il = dot_sequences(self.def_grad, stress_in_full)
        val_ij = dot_sequences(val_il, self.def_grad, use_rows=True)

        stress_out_full = val_ij / self.jacobian

        ii = get_sym_indices(self.dim)
        sh = stress_out_full.shape
        stress_out_full.shape = (sh[0], sh[1], sh[2] * sh[3])

        self._assert_symmetry(stress_out_full)

        stress_out = nm.empty_like(stress_in)
        stress_out[...,0] = stress_out_full[:,:,self.f2s]
        return stress_out
        
