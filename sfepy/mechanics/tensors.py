"""
Functions to compute some tensor-related quantities usual in continuum mechanics.
"""
import numpy as nm
import numpy.linalg as nla

from sfepy.base.base import assert_, Struct
from sfepy.linalg import dot_sequences, make_axis_rotation_matrix

def dim2sym(dim):
    """
    Given the space dimension, return the symmetric storage size.
    """
    return (dim + 1) * dim // 2

def sym2dim(sym):
    """
    Given the symmetric storage size, return the space dimension.

    Notes
    -----
    This function works for any space dimension.
    """
    val = int(-0.5 + nm.sqrt(2 * sym + 0.25))
    assert_(dim2sym(val) == sym)

    return val

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


def get_cauchy_strain(grad):
    """
    Given a gradient, return the corresponding Cauchy strain (symmetric
    gradient).
    """
    nc, _, dim = grad.shape
    sym = dim2sym(dim)
    out = nm.empty((nc, sym, 1), dtype=grad.dtype)

    strain_tab = {
        1: (0, ),
        2: (0, 3, (1, 2)),
        3: (0, 4, 8, (1, 3), (2, 6), (5, 7)),
    }

    vals_ = grad.reshape((nc, -1))
    for ii, idx in enumerate(strain_tab[dim]):
        if isinstance(idx, tuple):
            out[:, ii, 0] = nm.sum(vals_[:, nm.array(idx)], axis=1)
        else:
            out[:, ii, 0] = vals_[:, idx]

    return out


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

    if dim == 2:

        if sym_storage:
            s11 = stress[:,0]
            s22 = stress[:,1]
            s12 = stress[:,2]

        else:
            s11 = stress[:,0,0]
            s22 = stress[:,1,1]
            s12 = stress[:,0,1]

        vms = nm.sqrt(s11**2 - s11*s22 + s22**2 + 3*s12**2)[:,None]

    else:

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

def get_t4_from_t2s(t2s):
    """
    Get the full 4D tensor with major/minor symmetries from its 2D matrix
    representation.

    Parameters
    ----------
    t2s : array
        The symmetrically-stored tensor of shape (S, S), where S it the
        symmetric storage size.

    Returns
    -------
    t4 : array
        The full 4D tensor of shape (D, D, D, D), where D is the space
        dimension.
    """
    dim = sym2dim(t2s.shape[0])
    iif = get_full_indices(dim)

    t4 = t2s[:, iif][iif, ...]

    return t4

def prepare_cylindrical_transform(coors, origin, mode='axes'):
    """
    Prepare matrices for transforming tensors into cylindrical coordinates with
    the axis 'z' in a given origin.

    Parameters
    ----------
    coors : array
        The Cartesian coordinates.
    origin : array of length 3
        The origin.
    mode : 'axes' or 'data'
        In 'axes' (default) mode the matrix transforms data to different
        coordinate system, while in 'data' mode the matrix transforms
        the data in the same coordinate system and is transpose of the
        matrix in the 'axes' mode.

    Returns
    -------
    mtx : array
        The array of transformation matrices for each coordinate in `coors`.
    """
    assert_(mode in ['axes', 'data'])

    x, y = coors[:,0] - origin[0], coors[:,1] - origin[1]
    theta = nm.arctan2(y, x)
    if mode == 'data':
        theta = -theta

    mtx = nm.zeros((coors.shape[0], 3, 3), dtype=nm.float64)
    for ii, th in enumerate(theta):
        mtx[ii] = make_axis_rotation_matrix([0.0, 0.0, 1.0], th)

    return mtx

def transform_data(data, coors=None, mode='cylindrical', mtx=None):
    r"""
    Transform vector or tensor data components between orthogonal
    coordinate systems in 3D using transformation matrix :math:`M`, that
    should express rotation of the original coordinate system to the new
    system denoted by :math:`\bullet'` below.

    For vectors:

    .. math::
        \ul{v}' = M \cdot \ul{v}

    For second order tensors:

    .. math::
        \ull{t}' = M \cdot \ull{t} \cdot M^T

        \mbox{or}

        t_{ij}' = M_{ip} M_{jq} t_{pq}

    For fourth order tensors:

    .. math::

        t_{ijkl}' = M_{ip} M_{jq} M_{kr} M_{ls} t_{pqrs}

    Parameters
    ----------
    data : array, shape (num, n_r) or (num, n_r, n_c)
        The vectors (`n_r` is 3) or tensors (symmetric storage, `n_r` is 6,
        `n_c`, if available, is 1 or 6) to be transformed.
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

    if shape[0] != mtx.shape[0]:
        raise ValueError('incompatible numbers of points! (data: %d, mtx: %d)'
                         % (shape[0], mtx.shape[0]))

    if shape[1] == 3: # Vectors.
        new_data = dot_sequences(mtx, data)

    elif shape[1] == 6: # Symmetric tensors.
        iif = get_full_indices(3)
        iis = get_sym_indices(3)

        if ((data.ndim == 2)
            or ((data.ndim == 3) and (shape[2] == 1))): # Second order.
            if data.ndim == 3:
                aux = data[:, iif, 0]

            else:
                aux = data[:, iif]

            aux2 = dot_sequences(dot_sequences(mtx, aux, 'AB'), mtx, 'ABT')
            assert nm.allclose(aux2[0],
                               nm.dot(nm.dot(mtx[0], aux[0]), mtx[0].T))

            aux3 = aux2.reshape((aux2.shape[0], 9))

            new_data = aux3[:, iis]
            if data.ndim == 3:
                new_data = new_data[..., None]

        elif (data.ndim == 3) and (shape[2] == 6): # Fourth order.
            # Note: nm.einsum() is much slower than dot_sequences().
            df = data[:, iif][..., iif]
            tdf = nm.einsum('apqrs,aip,ajq,akr,als->aijkl',
                            df, mtx, mtx, mtx, mtx)
            tdf2 = tdf.reshape(tdf.shape[0], 9, 9)
            new_data = tdf2[:, :, iis][:, iis]

        else:
            raise ValueError('unsupported data shape! (%s)' % str(shape))


    else:
        raise ValueError('unsupported data shape! (%s)' % str(shape))

    assert_(new_data.shape == shape)

    return new_data

class StressTransform(Struct):
    """
    Encapsulates functions to convert various stress tensors in the symmetric
    storage given the deformation state.
    """

    def __init__(self, def_grad, jacobian=None):
        r"""
        Set :math:`\ull{F} = \pdiff{\ul{x}}{\ul{X}}` and optionally also
        :math:`J = \det(\ull{F})`.
        """
        self.def_grad = nm.asarray(def_grad, dtype=nm.float64)
        self.n_el, self.n_qp, self.dim = self.def_grad.shape[:3]

        self.s2f = get_full_indices(self.dim)
        self.f2s = get_sym_indices(self.dim)

        if jacobian is None:
            self.jacobian = nla.det(self.def_grad)[..., None, None]

        else:
            self.jacobian = nm.asarray(jacobian, dtype=nm.float64)

    def _assert_symmetry(self, stress):
        i1, i2 = get_non_diagonal_indices(self.dim)
        assert_(nm.allclose(stress[:,:,i1], stress[:,:,i2]))

    def get_1pk_from_2pk(self, stress_in):
        r"""
        Get the first Piola-Kirchhoff stress given the second Piola-Kirchhoff
        stress.

        .. math::

            P_{ij} = F_{ik} S_{kj}

        Parameters
        ----------
        stress_in : array_like
            The second Piola-Kirchhoff stress in vector symmetric storage with
            the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`

        Returns
        -------
        stress_out_full : array
            The first Piola-Kirchhoff stress in matrix storage.
        """
        stress_in = nm.asarray(stress_in, dtype=nm.float64)

        stress_in_full = stress_in[:,:,self.s2f,0]
        stress_out_full = dot_sequences(self.def_grad, stress_in_full)
        return stress_out_full

    def get_cauchy_from_2pk(self, stress_in):
        r"""
        Get the Cauchy stress given the second Piola-Kirchhoff stress.

        .. math::

            \sigma_{ij} = J^{-1} F_{ik} S_{kl} F_{jl}

        Parameters
        ----------
        stress_in : array_like
            The second Piola-Kirchhoff stress in vector symmetric storage with
            the indices ordered as :math:`[11, 22, 33, 12, 13, 23]`

        Returns
        -------
        stress_out : array
            The Cauchy stress in vector symmetric storage with the indices
            ordered as :math:`[11, 22, 33, 12, 13, 23]`.
        """
        stress_in = nm.asarray(stress_in, dtype=nm.float64)

        stress_in_full = stress_in[:,:,self.s2f,0]

        val_il = dot_sequences(self.def_grad, stress_in_full)
        val_ij = dot_sequences(val_il, self.def_grad, mode='ABT')

        stress_out_full = val_ij / self.jacobian

        sh = stress_out_full.shape
        stress_out_full.shape = (sh[0], sh[1], sh[2] * sh[3])

        self._assert_symmetry(stress_out_full)

        stress_out = nm.empty_like(stress_in)
        stress_out[...,0] = stress_out_full[:,:,self.f2s]
        return stress_out
