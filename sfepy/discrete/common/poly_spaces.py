from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import load_classes, Struct
from sfepy import get_paths

def transform_basis(transform, bf):
    """
    Transform a basis `bf` using `transform` array of matrices.
    """
    if bf.ndim == 3:
        nbf = nm.einsum('cij,qdj->cqdi', transform, bf, order='C')

    elif bf.ndim == 4:
        if bf.shape[0] == 1:
            nbf = nm.einsum('cij,qdj->cqdi', transform, bf[0], order='C')

        else:
            nbf = nm.einsum('cij,cqdj->cqdi', transform, bf, order='C')

    # Note: the 2nd derivatives are not supported here.
    # Workaround for NumPy 1.14.0 - order is ignored(?)
    nbf = nm.ascontiguousarray(nbf)

    return nbf

class PolySpace(Struct):
    """Abstract polynomial space class."""
    _all = None

    keys = {
        (0, 1) : 'simplex',
        (1, 2) : 'simplex',
        (2, 3) : 'simplex',
        (3, 4) : 'simplex',
        (2, 4) : 'tensor_product',
        (3, 8) : 'tensor_product',
    }

    @staticmethod
    def any_from_args(name, geometry, order, base='lagrange',
                      force_bubble=False):
        """
        Construct a particular polynomial space classes according to the
        arguments passed in.
        """
        if name is None:
            name = PolySpace.suggest_name(geometry, order, base, force_bubble)

        if PolySpace._all is None:
            ps_files = get_paths('sfepy/discrete/fem/poly_spaces.py')
            ps_files += get_paths('sfepy/discrete/dg/dg_poly_spaces.py')
            PolySpace._all = load_classes(ps_files, [PolySpace],
                                          ignore_errors=True,
                                          name_attr='name')
        table = PolySpace._all

        key = '%s_%s' % (base, PolySpace.keys[(geometry.dim,
                                               geometry.n_vertex)])
        if (geometry.name == '1_2') and (key not in table):
            key = '%s_%s' % (base, 'tensor_product')

        if force_bubble:
            key += '_bubble'

        return table[key](name, geometry, order)

    @staticmethod
    def suggest_name(geometry, order, base='lagrange',
                     force_bubble=False):
        """
        Suggest the polynomial space name given its constructor parameters.
        """
        aux = geometry.get_interpolation_name()[:-1]
        if force_bubble:
            return aux + ('%dB' % order)
        else:
            return aux + ('%d' % order)

    def __init__(self, name, geometry, order):
        self.name = name
        self.geometry = geometry
        self.order = order

        self.bbox = nm.vstack((geometry.coors.min(0), geometry.coors.max(0)))

    def eval_base(self, coors, diff=0, ori=None, force_axis=False,
                  transform=None, suppress_errors=False, eps=1e-15):
        """
        Evaluate the basis or its first or second derivatives in points given
        by coordinates. The real work is done in _eval_base() implemented in
        subclasses.

        Note that the second derivative code is a work-in-progress and only
        `coors` and `transform` arguments are used.

        Parameters
        ----------
        coors : array_like
            The coordinates of points where the basis is evaluated. See Notes.
        diff : 0, 1 or 2
            If nonzero, return the given derivative.
        ori : array_like, optional
            Optional orientation of element facets for per element basis.
        force_axis : bool
            If True, force the resulting array shape to have one more axis even
            when `ori` is None.
        transform : array_like, optional
            The basis transform array.
        suppress_errors : bool
            If True, do not report points outside the reference domain.
        eps : float
            Accuracy for comparing coordinates.

        Returns
        -------
        base : array
            The basis (shape (n_coor, 1, n_base)) or its first derivative
            (shape (n_coor, dim, n_base)) or its second derivative (shape
            (n_coor, dim, dim, n_base)) evaluated in the given points. An
            additional axis is pre-pended of length n_cell, if `ori` is given,
            or of length 1, if `force_axis` is True.

        Notes
        -----
        If coors.ndim == 3, several point sets are assumed, with equal number
        of points in each of them. This is the case, for example, of the
        values of the volume base functions on the element facets. The indexing
        (of bf_b(g)) is then (ifa,iqp,:,n_ep), so that the facet can be set in
        C using FMF_SetCell.
        """
        coors = nm.asarray(coors)
        if not coors.ndim in (2, 3):
            raise ValueError('coordinates must have 2 or 3 dimensions! (%d)'
                             % coors.ndim)

        if (coors.ndim == 2):
            base = self._eval_base(coors, diff=diff, ori=ori,
                                   suppress_errors=suppress_errors,
                                   eps=eps)

            if (base.ndim == 3) and force_axis:
                base = base[None, ...]

            if not base.flags['C_CONTIGUOUS']:
                base = nm.ascontiguousarray(base)

        else: # Several point sets.
            if diff:
                bdim = self.geometry.dim
            else:
                bdim = 1

            base = nm.empty((coors.shape[0], coors.shape[1],
                             bdim, self.n_nod), dtype=nm.float64)

            for ii, _coors in enumerate(coors):
                base[ii] = self._eval_base(_coors, diff=diff, ori=ori,
                                           suppress_errors=suppress_errors,
                                           eps=eps)

        if transform is not None:
            base = transform_basis(transform, base)

        return base
