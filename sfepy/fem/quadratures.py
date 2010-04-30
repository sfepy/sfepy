"""
`quadrature_tables` are organized as follows::

    quadrature_tables = {
        '<geometry1>' : {
            order1 : QuadraturePoints(args1),
            order2 : QuadraturePoints(args2),
            ...
        },
        '<geometry2>' : {
            order1 : QuadraturePoints(args1),
            order2 : QuadraturePoints(args2),
            ...
        },
        ...
    }

Naming conventions in problem description files::

    `<family>_<order>_<dimension>`

Integral 'family' is just an arbitrary name given by user.

Quadrature coordinates and weights copied from "The Finite Element
Method Displayed" by Gouri Dhatt and Gilbert Touzat,
Wiley-Interscience Production, 1984. 

Examples
--------
* gauss_o2_d2 # second order, 2D
* gauss_o1_d3 # first order, 3D
* my_int_o1_d3 # same as above
"""
from sfepy.base.base import *

class QuadraturePoints(Struct):
    """
    Representation of a set of quadrature points.

    Parameters
    ----------
    data : array_like
        The array of shape `(n_point, dim + 1)` of quadrature point
        coordinates (first `dim` columns) and weights (the last column). 
    coors : array_like, optional
        Optionally, instead of using `data`, the coordinates and weights can
        be provided separately - `data` are then ignored.
    weights : array_like, optional
        Optionally, instead of using `data`, the coordinates and weights can
        be provided separately - `data` are then ignored.
    bounds : (float, float), optional
        The coordinates and weights should correspond to a reference
        element in `[0, 1]` x `dim`. Provide the correct bounds if this is
        not the case.
    tp_fix : float, optional
        The value that is used to multiply the tensor product element
        volume (= 1.0) to get the correct volume.
    """

    def __init__(self, data, coors=None, weights=None, bounds=None, tp_fix=1.0):
        if coors is None:
            data = nm.array(data, dtype=nm.float64, ndmin=2)
            self.coors = data[:,:-1].copy()
            self.weights = data[:,-1].copy()

        elif weights is not None:
            self.coors = nm.array(coors, dtype=nm.float64, ndmin=2)
            self.weights = nm.array(weights, dtype=nm.float64)

        else:
            raise ValueError('both "coors" and "weights" have to be provided!')

        self.n_point, self.dim = self.coors.shape

        self.bounds = (0, 1)
        bbox = nm.array([self.bounds] * self.dim, dtype=nm.float64)
        self.volume = nm.prod(bbox.sum(axis=1)) * tp_fix

        if bounds is not None:
            # Transform from given bounds to self.bounds.
            bbox = nm.array([bounds] * self.dim, dtype=nm.float64)
            volume = nm.prod(nm.diff(bbox, axis=1)) * tp_fix

            a, b = bounds
            c, d = self.bounds

            c1 = (d - c) / (b - a)
            c2 = ((b * c) - (a * d)) / (b - a)

            self.coors = c1 * self.coors + c2
            self.weights *= self.volume / volume

_QP = QuadraturePoints
quadrature_tables = {
    '1_2' : {

        1 : _QP( [ 0.5, 1.0]),

        3 : _QP([[-1.0/nm.sqrt(3), 1.0],
                 [ 1.0/nm.sqrt(3), 1.0]], bounds=(-1.0, 1.0)),

        5 : _QP([[-nm.sqrt(3.0/5.0), 5.0/9.0],
                 [ 0.0             , 8.0/9.0],
                 [ nm.sqrt(3.0/5.0), 5.0/9.0]], bounds=(-1.0, 1.0)),

        7 : _QP([[-nm.sqrt((3+2*nm.sqrt(6.0/5.0))/7), 0.5-1/(6*nm.sqrt(6.0/5.0))],
                 [-nm.sqrt((3-2*nm.sqrt(6.0/5.0))/7), 0.5+1/(6*nm.sqrt(6.0/5.0))],
                 [ nm.sqrt((3-2*nm.sqrt(6.0/5.0))/7), 0.5+1/(6*nm.sqrt(6.0/5.0))],
                 [ nm.sqrt((3+2*nm.sqrt(6.0/5.0))/7), 0.5-1/(6*nm.sqrt(6.0/5.0))]], bounds=(-1.0, 1.0)),

        9 : _QP([[-nm.sqrt(5+4*nm.sqrt(5.0/14.0))/3, 161.0/450.0-13/(180*nm.sqrt(5.0/14.0))],
                 [-nm.sqrt(5-4*nm.sqrt(5.0/14.0))/3, 161.0/450.0+13/(180*nm.sqrt(5.0/14.0))],
                 [ 0.0                             , 128.0/225.0                           ],
                 [ nm.sqrt(5-4*nm.sqrt(5.0/14.0))/3, 161.0/450.0+13/(180*nm.sqrt(5.0/14.0))],
                 [ nm.sqrt(5+4*nm.sqrt(5.0/14.0))/3, 161.0/450.0-13/(180*nm.sqrt(5.0/14.0))]], bounds=(-1.0, 1.0)),

    },

    '2_3' : {

        1 : _QP([[ 1.0/3.0, 1.0/3.0, 0.5]], tp_fix=0.5),

        2 : _QP([[ 1.0/6.0, 1.0/6.0, 1.0/6.0],
                 [ 2.0/3.0, 1.0/6.0, 1.0/6.0],
                 [ 1.0/6.0, 2.0/3.0, 1.0/6.0]], tp_fix=0.5),

        3 : _QP([[ 1.0/3.0, 1.0/3.0,-27.0/96.0],
                 [ 1.0/5.0, 1.0/5.0, 25.0/96.0],
                 [ 3.0/5.0, 1.0/5.0, 25.0/96.0],
                 [ 1.0/5.0, 3.0/5.0, 25.0/96.0]], tp_fix=0.5),

    },

    '2_4' : {

        2 : _QP([[ nm.sqrt(2.0/3.0), 0.0         , 4.0/3.0],
                 [-1/nm.sqrt(6)    , 1/nm.sqrt(2), 4.0/3.0],
                 [-1/nm.sqrt(6)    ,-1/nm.sqrt(2), 4.0/3.0]], bounds=(-1.0, 1.0)),

        3 : _QP([[-1/nm.sqrt(3),-1/nm.sqrt(3), 1.0],
                 [ 1/nm.sqrt(3),-1/nm.sqrt(3), 1.0],
                 [ 1/nm.sqrt(3), 1/nm.sqrt(3), 1.0],
                 [-1/nm.sqrt(3), 1/nm.sqrt(3), 1.0]], bounds=(-1.0, 1.0)),

        5 : _QP([[ 0.0             , 0.0               , 8.0/7.0],
                 [ 0.0             ,-nm.sqrt(14.0/15.0), 20.0/63.0],
                 [ 0.0             , nm.sqrt(14.0/15.0), 20.0/63.0],
                 [-nm.sqrt(3.0/5.0),-nm.sqrt(3.0/5.0)  , 20.0/36.0],
                 [-nm.sqrt(3.0/5.0), nm.sqrt(3.0/5.0)  , 20.0/36.0],
                 [ nm.sqrt(3.0/5.0),-nm.sqrt(3.0/5.0)  , 20.0/36.0],
                 [ nm.sqrt(3.0/5.0), nm.sqrt(3.0/5.0)  , 20.0/36.0]], bounds=(-1.0, 1.0)),

    },

    '3_4' : {
        1 : _QP([[ 1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/6.0]], tp_fix=1.0/6.0),

        2 : _QP([[ (5-nm.sqrt(5))/20  , (5-nm.sqrt(5))/20  , (5-nm.sqrt(5))/20  , 1.0/24.0],
                 [ (5-nm.sqrt(5))/20  , (5-nm.sqrt(5))/20  , (5+3*nm.sqrt(5))/20, 1.0/24.0],
                 [ (5-nm.sqrt(5))/20  , (5+3*nm.sqrt(5))/20, (5-nm.sqrt(5))/20  , 1.0/24.0],
                 [ (5+3*nm.sqrt(5))/20, (5-nm.sqrt(5))/20  , (5-nm.sqrt(5))/20  , 1.0/24.0]], tp_fix=1.0/6.0),

        3 : _QP([[ 1.0/4.0, 1.0/4.0, 1.0/4.0,-2.0/15.0],
                 [ 1.0/6.0, 1.0/6.0, 1.0/6.0, 3.0/40.0],
                 [ 1.0/6.0, 1.0/6.0, 1.0/2.0, 3.0/40.0],
                 [ 1.0/6.0, 1.0/2.0, 1.0/6.0, 3.0/40.0],
                 [ 1.0/2.0, 1.0/6.0, 1.0/6.0, 3.0/40.0]], tp_fix=1.0/6.0),

    },

    '3_8' : {

        2 : _QP([[ 0.0             , nm.sqrt(2.0/3.0),-1/nm.sqrt(3), 2.0],
                 [ 0.0             ,-nm.sqrt(2.0/3.0),-1/nm.sqrt(3), 2.0],
                 [ nm.sqrt(2.0/3.0), 0.0             , 1/nm.sqrt(3), 2.0],
                 [-nm.sqrt(2.0/3.0), 0.0             , 1/nm.sqrt(3), 2.0]], bounds=(-1.0, 1.0)),

        3 : _QP([[-1.0, 0.0, 0.0, 4.0/3.0],
                 [ 1.0, 0.0, 0.0, 4.0/3.0],
                 [ 0.0,-1.0, 0.0, 4.0/3.0],
                 [ 0.0, 1.0, 0.0, 4.0/3.0],
                 [ 0.0, 0.0,-1.0, 4.0/3.0],
                 [ 0.0, 0.0, 1.0, 4.0/3.0]], bounds=(-1.0, 1.0)),

        5 : _QP([[-nm.sqrt(19.0/30.0), 0.0               , 0.0               , 320.0/361.0],
                 [ nm.sqrt(19.0/30.0), 0.0               , 0.0               , 320.0/361.0],
                 [ 0.0               ,-nm.sqrt(19.0/30.0), 0.0               , 320.0/361.0],
                 [ 0.0               , nm.sqrt(19.0/30.0), 0.0               , 320.0/361.0],
                 [ 0.0               , 0.0               ,-nm.sqrt(19.0/30.0), 320.0/361.0],
                 [ 0.0               , 0.0               , nm.sqrt(19.0/30.0), 320.0/361.0],
                 [ nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0), 121.0/361.0],
                 [ nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0), 121.0/361.0],
                 [ nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0), 121.0/361.0],
                 [ nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0), 121.0/361.0],
                 [-nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0), 121.0/361.0],
                 [-nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0), 121.0/361.0],
                 [-nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0), nm.sqrt(19.0/33.0), 121.0/361.0],
                 [-nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0),-nm.sqrt(19.0/33.0), 121.0/361.0]], bounds=(-1.0, 1.0)),

    },
}
del _QP
