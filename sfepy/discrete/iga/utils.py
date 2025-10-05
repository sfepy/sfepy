"""
Utility functions based on igakit.
"""
import numpy as nm

from sfepy.base.base import Struct
from sfepy.discrete.fem import Mesh
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

    coors = nurbs(*pars)
    coors.shape = (-1, coors.shape[-1])

    conn, desc = get_tensor_product_conn([len(ii) for ii in pars])

    if (coors[:, -1] == 0.0).all():
        coors = coors[:, :-1]

    return coors, conn, desc

def create_mesh_and_output(nurbs, pars=None, **kwargs):
    """
    Create a nD-linear tensor product FE mesh using
    :func:`create_linear_fe_mesh()`, evaluate field variables given as keyword
    arguments in the mesh vertices and create a dictionary of output data
    usable by Mesh.write().

    Parameters
    ----------
    nurbs : igakit.nurbs.NURBS instance
        The NURBS object.
    pars : sequence of array, optional
        The values of parameters in each parametric dimension. If not given,
        the values are set so that the resulting mesh has the same number of
        vertices as the number of control points/basis functions of the NURBS
        object.
    **kwargs : kwargs
        The field variables as keyword arguments. Their names serve as keys in
        the output dictionary.

    Returns
    -------
    mesh : Mesh instance
        The finite element mesh.
    out : dict
        The output dictionary.
    """
    coors, conn, desc = create_linear_fe_mesh(nurbs, pars)
    mat_id = nm.zeros(conn.shape[0], dtype=nm.int32)
    mesh = Mesh.from_data('nurbs', coors, None, [conn], [mat_id], [desc])

    out = {}
    for key, variable in kwargs.items():
        if variable.ndim == 2:
            nc = variable.shape[1]
            field = variable.reshape(nurbs.weights.shape + (nc,))

        else:
            field = variable.reshape(nurbs.weights.shape)
            nc = 1

        vals = nurbs.evaluate(field, *pars)
        out[key] = Struct(name='output_data', mode='vertex',
                          data=vals.reshape((-1, nc)))

    return mesh, out

def save_basis(nurbs, pars):
    """
    Save a NURBS object basis on a FE mesh corresponding to the given
    parametrization in VTK files.

    Parameters
    ----------
    nurbs : igakit.nurbs.NURBS instance
        The NURBS object.
    pars : sequence of array, optional
        The values of parameters in each parametric dimension.
    """
    coors, conn, desc = create_linear_fe_mesh(nurbs, pars)
    mat_id = nm.zeros(conn.shape[0], dtype=nm.int32)
    mesh = Mesh.from_data('nurbs', coors, None, [conn], [mat_id], [desc])

    n_dof = nurbs.weights.ravel().shape[0]
    variable = nm.zeros(n_dof, dtype=nm.float64)
    field = variable.reshape(nurbs.weights.shape)
    for ic in range(n_dof):
        variable[ic - 1] = 0.0
        variable[ic] = 1.0

        vals = nurbs.evaluate(field, *pars).reshape((-1))
        out = {}
        out['bf'] = Struct(name='output_data', mode='vertex',
                           data=vals[:, None])
        mesh.write('iga_basis_%03d.vtk' % ic, io='auto', out=out)
