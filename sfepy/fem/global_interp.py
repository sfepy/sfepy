"""
Global interpolation functions.
"""
import time
import numpy as nm

from sfepy.base.base import output, get_default_attr
from sfepy.fem.mesh import make_inverse_connectivity
from sfepy.fem.extmods.bases import find_ref_coors

def get_ref_coors(field, coors, strategy='kdtree', close_limit=0.1, cache=None,
                  verbose=True):
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

        scoors = mesh.coors
        output('reference field: %d vertices' % scoors.shape[0],
               verbose=verbose)

        iconn = get_default_attr(cache, 'iconn', None)
        if iconn is None:
            offsets, iconn = make_inverse_connectivity(mesh.conns,
                                                       mesh.n_nod,
                                                       ret_offsets=True)

            ii = nm.where(offsets[1:] == offsets[:-1])[0]
            if len(ii):
                raise ValueError('some vertices not in any element! (%s)'
                                 % ii)

        else:
            offsets = cache.offsets

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

            vertex_coorss, nodess, mtx_is = [], [], []
            conns = []
            for ig, ap in field.aps.iteritems():
                ps = ap.interp.gel.interp.poly_spaces['v']

                vertex_coorss.append(ps.geometry.coors)
                nodess.append(ps.nodes)
                mtx_is.append(ps.get_mtx_i())

                conns.append(mesh.conns[ig].copy())

            # Get reference element coordinates corresponding to
            # destination coordinates.
            ref_coors = nm.empty_like(coors)
            cells = nm.empty((coors.shape[0], 2), dtype=nm.int32)
            status = nm.empty((coors.shape[0],), dtype=nm.int32)

            find_ref_coors(ref_coors, cells, status, coors,
                        ics, offsets, iconn,
                        scoors, conns,
                        vertex_coorss, nodess, mtx_is,
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
