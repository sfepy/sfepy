"""
Global interpolation functions.
"""
import time
import numpy as nm

from sfepy.base.base import output, get_default_attr
from sfepy.discrete.fem.geometry_element import create_geometry_elements
from sfepy.discrete.fem.extmods.crefcoors import find_ref_coors

def get_ref_coors(field, coors, strategy='kdtree', close_limit=0.1, cache=None,
                  verbose=False):
    """
    Get reference element coordinates and elements corresponding to given
    physical coordinates.

    Parameters
    ----------
    field : Field instance
        The field defining the approximation.
    coors : array
        The physical coordinates.
    strategy : str, optional
        The strategy for finding the elements that contain the
        coordinates. Only 'kdtree' is supported for the moment.
    close_limit : float, optional
        The maximum limit distance of a point from the closest
        element allowed for extrapolation.
    cache : Struct, optional
        To speed up a sequence of evaluations, the field mesh and other data
        can be cached. Optionally, the cache can also contain the reference
        element coordinates as `cache.ref_coors`, `cache.cells` and
        `cache.status`, if the evaluation occurs in the same coordinates
        repeatedly. In that case the mesh related data are ignored.
    verbose : bool
        If False, reduce verbosity.

    Returns
    -------
    ref_coors : array
        The reference coordinates.
    cells : array
        The cell indices corresponding to the reference coordinates.
    status : array
        The status: 0 is success, 1 is extrapolation within `close_limit`, 2 is
        extrapolation outside `close_limit`, 3 is failure, 4 is failure due to
        non-convergence of the Newton iteration in tensor product cells.

    Notes
    -----
    Outline of the algorithm for finding xi such that X(xi) = P:

    1. make inverse connectivity - for each vertex have cells it is in.
    2. find the closest vertex V.
    3. choose initial cell: i0 = first from cells incident to V.
    4. while not P in C_i, change C_i towards P, check if P in new C_i.
    """
    ref_coors = get_default_attr(cache, 'ref_coors', None)

    if ref_coors is None:
        ref_coors = nm.empty_like(coors)
        cells = nm.empty((coors.shape[0],), dtype=nm.int32)
        status = nm.empty((coors.shape[0],), dtype=nm.int32)

        cmesh = get_default_attr(cache, 'cmesh', None)
        if cmesh is None:
            mesh = field.create_mesh(extra_nodes=False)

            tt = time.clock()
            mesh = field.create_mesh(extra_nodes=False)
            cmesh = mesh.cmesh

            gels = create_geometry_elements()

            cmesh.set_local_entities(gels)
            cmesh.setup_entities()

            centroids = cmesh.get_centroids(cmesh.tdim)

            if field.gel.name != '3_8':
                normals0 = cmesh.get_facet_normals()
                normals1 = None

            else:
                normals0 = cmesh.get_facet_normals(0)
                normals1 = cmesh.get_facet_normals(1)

            output('cmesh setup: %f s' % (time.clock()-tt), verbose=verbose)

        else:
            centroids = cache.centroids
            normals0 = cache.normals0
            normals1 = cache.normals1

        kdtree = get_default_attr(cache, 'kdtree', None)
        if kdtree is None:
            from scipy.spatial import cKDTree as KDTree

            tt = time.clock()
            kdtree = KDTree(cmesh.coors)
            output('kdtree: %f s' % (time.clock()-tt), verbose=verbose)

        tt = time.clock()
        ics = kdtree.query(coors)[1]
        output('kdtree query: %f s' % (time.clock()-tt), verbose=verbose)

        ics = nm.asarray(ics, dtype=nm.int32)

        ap = field.ap
        ps = ap.interp.gel.interp.poly_spaces['v']
        mtx_i = ps.get_mtx_i()

        ac = nm.ascontiguousarray

        tt = time.clock()
        find_ref_coors(ref_coors, cells, status, ac(coors), cmesh, centroids,
                       normals0, normals1, ics,
                       ps.geometry.coors, ps.nodes, mtx_i,
                       1, close_limit, 1e-15, 100, 1e-8)
        output('ref. coordinates: %f s' % (time.clock()-tt), verbose=verbose)

    else:
        cells = cache.cells
        status = cache.status

    return ref_coors, cells, status
