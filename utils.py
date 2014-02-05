"""
Utility functions based on igakit.
"""
from itertools import product

import numpy as nm

from sfepy.mesh.mesh_generators import get_tensor_product_conn

def create_linear_fe_mesh(nurbs, pars=None):
    """
    Convert a NURBS object into a nD-linear tensor product FE mesh.

    Parameters
    ----------
    nurbs : igakit.nurbs.NURBS instance
        The NURBS object.
    pars : sequence of array, optional
        The values of parameters in each parametric dimension. If not given,
        the values are set so that the resulting mesh has the same number of
        vertices as the number of control points/basis functions of the NURBS
        object.

    Returns
    -------
    coors : array
        The coordinates of mesh vertices.
    conn : array
        The vertex connectivity array.
    desc : str
        The cell kind.
    """
    knots = nurbs.knots
    shape = nurbs.weights.shape

    if pars is None:
        pars = []
        for ii, kv in enumerate(knots):
            par = nm.linspace(kv[0], kv[-1], shape[ii])
            pars.append(par)

    coors = nm.array([nurbs(*ii) for ii in product(*pars)])
    conn, desc = get_tensor_product_conn([len(ii) for ii in pars])

    if (coors[:, -1] == 0.0).all():
        coors = coors[:, :-1]

    return coors, conn, desc
