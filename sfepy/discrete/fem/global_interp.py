"""
Global interpolation functions.
"""
import time
import numpy as nm

from sfepy.base.base import output, get_default_attr
from sfepy.discrete.fem.extmods.bases import find_ref_coors

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
        To speed up a sequence of evaluations, the field mesh, the inverse
        connectivity of the field mesh and the KDTree instance can be cached as
        `cache.mesh`, `cache.offsets`, `cache.iconn` and
        `cache.kdtree`. Optionally, the cache can also contain the reference
        element coordinates as `cache.ref_coors`, `cache.cells` and
        `cache.status`, if the evaluation occurs in the same coordinates
        repeatedly. In that case the KDTree related data are ignored.
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
        extrapolation outside `close_limit`, 3 is failure.
    """
    ref_coors = get_default_attr(cache, 'ref_coors', None)
    if ref_coors is None:
        mesh = get_default_attr(cache, 'mesh', None)
        if mesh is None:
            mesh = field.create_mesh(extra_nodes=False)

        cmesh = mesh.cmesh
        cconn = cmesh.get_conn(cmesh.tdim, 0)
        conn = cconn.indices.reshape((cmesh.n_el, -1)).astype(nm.int32)

        cmesh.setup_connectivity(0, cmesh.tdim)
        ciconn = cmesh.get_conn(0, cmesh.tdim)
        iconn = ciconn.indices.astype(nm.int32)
        offsets = ciconn.offsets.astype(nm.int32)

        scoors = mesh.coors
        output('reference field: %d vertices' % scoors.shape[0],
               verbose=verbose)

        if strategy == 'kdtree':
            kdtree = get_default_attr(cache, 'kdtree', None)
            if kdtree is None:
                from scipy.spatial import cKDTree as KDTree

                tt = time.clock()
                kdtree = KDTree(scoors)
                output('kdtree: %f s' % (time.clock()-tt), verbose=verbose)

            tt = time.clock()
            ics = kdtree.query(coors)[1]
            output('kdtree query: %f s' % (time.clock()-tt), verbose=verbose)

            tt = time.clock()
            ics = nm.asarray(ics, dtype=nm.int32)

            ap = field.ap
            ps = ap.interp.gel.interp.poly_spaces['v']
            mtx_i = ps.get_mtx_i()

            # Get reference element coordinates corresponding to
            # destination coordinates.
            ref_coors = nm.empty_like(coors)
            cells = nm.empty((coors.shape[0],), dtype=nm.int32)
            status = nm.empty((coors.shape[0],), dtype=nm.int32)

            find_ref_coors(ref_coors, cells, status, coors,
                        ics, offsets, iconn, scoors, conn,
                        ps.geometry.coors, ps.nodes, mtx_i,
                        1, close_limit, 1e-15, 100, 1e-8)
            output('ref. coordinates: %f s' % (time.clock()-tt),
                   verbose=verbose)

        elif strategy == 'crawl':
            raise NotImplementedError

        else:
            raise ValueError('unknown search strategy! (%s)' % strategy)

    else:
        ref_coors = cache.ref_coors
        cells = cache.cells
        status = cache.status

    return ref_coors, cells, status
