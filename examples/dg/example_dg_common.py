"""
Functions common to DG examples
"""
import os
from glob import glob

import numpy as nm

from sfepy.base.base import output
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
from sfepy.mesh.mesh_generators import gen_block_mesh


diffusion_schemes_implicit = {
    "symmetric":
          " + dw_dg_diffusion_flux.i.Omega(D.val, p, v)"
        + " + dw_dg_diffusion_flux.i.Omega(D.val, v, p)",
    "non-symmetric":
          " + dw_dg_diffusion_flux.i.Omega(D.val, p, v)"
        + " - dw_dg_diffusion_flux.i.Omega(D.val, v, p)",
    "incomplete":
        " + dw_dg_diffusion_flux.i.Omega(D.val, p, v)"}

diffusion_schemes_explicit = {
    "symmetric":
          " - dw_dg_diffusion_flux.i.Omega(D.val, p[-1], v)"
        + " - dw_dg_diffusion_flux.i.Omega(D.val, v, p[-1])",
    "non-symmetric":
          " - dw_dg_diffusion_flux.i.Omega(D.val, p[-1], v)"
        + " + dw_dg_diffusion_flux.i.Omega(D.val, v, p[-1])",
    "incomplete":
        " - dw_dg_diffusion_flux.i.Omega(D.val, p[-1], v)"}

functions = {}
def local_register_function(fun):
    try:
        functions.update({fun.__name__: (fun,)})

    except AttributeError:  # Already a sfepy Function.
        fun = fun.function
        functions.update({fun.__name__: (fun,)})

    return fun


def get_cfl_setup(CFL=None, dt=None):
    """
    Provide either CFL or dt to create preprocess hook that sets up
    Courant-Friedrichs-Levi stability condition for either advection or
    diffusion.

    Parameters
    ----------
    CFL : float, optional
    dt: float, optional

    Returns
    -------
    setup_cfl_condition : callable
        expects sfepy.discrete.problem as argument

    """

    if CFL is None and dt is None:
        raise ValueError("Specifiy either CFL or dt in CFL setup")

    def setup_cfl_condition(problem):
        """
        Sets up CFL condition for problem ts_conf in problem

        Parameters
        ----------
            problem : discrete.problem.Problem
        """
        ts_conf = problem.ts_conf
        mesh = problem.domain.mesh
        dim = mesh.dim
        first_field = list(problem.fields.values())[0]
        first_field_name = list(problem.fields.keys())[0]
        approx_order = first_field.approx_order
        mats = problem.create_materials(['a', 'D'])
        try:
            # make this more general?
            #  maybe require material name in parameter
            velo = problem.conf_materials['material_a__0'].values["val"]
            max_velo = nm.max(nm.linalg.norm(velo))
        except KeyError:
            max_velo = 1

        try:
            # make this more general?
            #  maybe require material name in parameter
            diffusion = problem.conf_materials['material_D__0'].values["val"]
            max_diffusion = nm.max(nm.linalg.norm(diffusion))
        except KeyError:
            max_diffusion = None

        dx = nm.min(problem.domain.mesh.cmesh.get_volumes(dim))

        output("Preprocess hook - setup_cfl_condition:...")
        output("Approximation order of field {}({}) is {}"
               .format(first_field_name, first_field.family_name, approx_order))
        output("Space divided into {0} cells, {1} steps, step size {2}"
               .format(mesh.n_el, len(mesh.coors), dx))

        if dt is None:
            adv_dt = get_cfl_advection(max_velo, dx, approx_order, CFL)
            diff_dt = get_cfl_diffusion(max_diffusion, dx, approx_order, CFL)
            _dt = min(adv_dt, diff_dt)
        else:
            output("CFL coefficient {0} ignored, dt specified directly"
                   .format(CFL))
            _dt = dt

        tn = int(nm.ceil((ts_conf.t1 - ts_conf.t0) / _dt))
        dtdx = _dt / dx

        ts_conf.dt = _dt
        ts_conf.n_step = tn
        ts_conf.cour = max_velo * dtdx

        output("Time divided into {0} nodes, {1} steps, step size is {2}"
               .format(tn - 1, tn, _dt))
        output("Courant number c = max(norm(a)) * dt/dx = {0}"
               .format(ts_conf.cour))
        output("Time stepping solver is {}".format(ts_conf.kind))
        output("... CFL setup done.")

    return setup_cfl_condition


def get_cfl_advection(max_velo, dx, approx_order, CFL):
    """

    Parameters
    ----------
    max_velo : float
    dx : float
    approx_order : int
    CFL : CFL

    Returns
    -------
    dt : float
    """
    order_corr = 1. / (2 * approx_order + 1)


    dt = dx / max_velo * CFL * order_corr

    if not (nm.isfinite(dt)):
        dt = 1
    output(("CFL advection: CFL coefficient was {0} " +
           "and order correction 1/{1} = {2}")
           .format(CFL, (2 * approx_order + 1), order_corr))
    output("CFL advection: resulting dt={}".format((dt)))
    return dt


def get_cfl_diffusion(max_diffusion, dx, approx_order, CFL,
                      do_order_corr=False):
    """

    Parameters
    ----------
    max_diffusion : float
    dx : float
    approx_order : int
    CFL : float
    do_order_corr : bool


    Returns
    -------
    dt : float
    """
    if max_diffusion is None:
        return 1

    if do_order_corr:
        order_corr = 1. / (2 * approx_order + 1)
    else:
        order_corr = 1

    dt = dx**2 / max_diffusion * CFL * order_corr

    if not (nm.isfinite(dt)):
        dt = 1
    output(("CFL diffusion: CFL coefficient was {0} " +
            "and order correction 1/{1} = {2}")
           .format(CFL, (2 * approx_order + 1), order_corr))
    output("CFL diffusion: resulting dt={}".format(dt))
    return dt


def get_gen_1D_mesh_hook(XS, XE, n_nod):
    """

    Parameters
    ----------
    XS : float
        leftmost coordinate
    XE : float
        rightmost coordinate
    n_nod : int
        number of nodes, number of cells is then n_nod - 1

    Returns
    -------
    mio : UserMeshIO instance
    """
    def mesh_hook(mesh, mode):
        """
        Generate the 1D mesh.
        """
        if mode == 'read':

            coors = nm.linspace(XS, XE, n_nod).reshape((n_nod, 1))
            conn = nm.arange(n_nod, dtype=nm.int32) \
                       .repeat(2)[1:-1].reshape((-1, 2))
            mat_ids = nm.zeros(n_nod - 1, dtype=nm.int32)
            descs = ['1_2']

            mesh = Mesh.from_data('uniform_1D{}'.format(n_nod), coors, None,
                                  [conn], [mat_ids], descs)
            return mesh

        elif mode == 'write':
            pass

    mio = UserMeshIO(mesh_hook)
    return mio

def get_gen_block_mesh_hook(dims, shape, centre, mat_id=0, name='block',
                            coors=None, verbose=True):
    """

    Parameters
    ----------
    dims : array of 2 or 3 floats
        Dimensions of the block.
    shape : array of 2 or 3 ints
        Shape (counts of nodes in x, y, z) of the block mesh.
    centre : array of 2 or 3 floats
        Centre of the block.
    mat_id : int, optional
        The material id of all elements.
    name : string
        Mesh name.
    verbose : bool
        If True, show progress of the mesh generation.

    Returns
    -------
    mio : UserMeshIO instance
    """
    def mesh_hook(mesh, mode):
        """
        Generate the 1D mesh.
        """
        if mode == 'read':

            mesh = gen_block_mesh(dims, shape, centre, mat_id=mat_id, name=name,
                                  coors=coors, verbose=verbose)
            return mesh

        elif mode == 'write':
            pass

    mio = UserMeshIO(mesh_hook)
    return mio

def clear_folder(clear_format, confirm=False, doit=True):
    """
    Deletes files matching the format

    Parameters
    ----------
    clear_format : str
    confirm : bool
    doit : bool
        if False do not delete anything no matter the confirmation

    Returns
    -------
    deleted_anything :
        True if there was something to delete
    """
    files = glob(clear_format)
    if confirm:
        for file in files:
            output("Will delete file {}".format(file))
        doit = input("--------------\nDelete files [Y/n]? ").strip() == "Y"

    if doit:
        for file in files:
            os.remove(file)
    return bool(files)
