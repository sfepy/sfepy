"""
Global interpolation functions.
"""
import numpy as nm

from sfepy.base.base import assert_, output, get_default_attr
from sfepy.base.timing import Timer
from sfepy.discrete.fem.geometry_element import create_geometry_elements
import sfepy.discrete.common.extmods.crefcoors as crc

def get_ref_coors_convex(field, coors, close_limit=0.1, cache=None,
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
    timer = Timer()

    ref_coors = get_default_attr(cache, 'ref_coors', None)
    if ref_coors is None:
        extrapolate = close_limit > 0.0

        ref_coors = nm.empty_like(coors)
        cells = nm.empty((coors.shape[0],), dtype=nm.int32)
        status = nm.empty((coors.shape[0],), dtype=nm.int32)

        cmesh = get_default_attr(cache, 'cmesh', None)
        if cmesh is None:
            timer.start()
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

            output('cmesh setup: %f s' % timer.stop(), verbose=verbose)

        else:
            centroids = cache.centroids
            normals0 = cache.normals0
            normals1 = cache.normals1

        kdtree = get_default_attr(cache, 'kdtree', None)
        if kdtree is None:
            from scipy.spatial import cKDTree as KDTree

            timer.start()
            kdtree = KDTree(cmesh.coors)
            output('kdtree: %f s' % timer.stop(), verbose=verbose)

        timer.start()
        ics = kdtree.query(coors)[1]
        output('kdtree query: %f s' % timer.stop(), verbose=verbose)

        ics = nm.asarray(ics, dtype=nm.int32)

        coors = nm.ascontiguousarray(coors)
        ctx = field.create_basis_context()

        timer.start()
        crc.find_ref_coors_convex(ref_coors, cells, status, coors, cmesh,
                                  centroids, normals0, normals1, ics,
                                  extrapolate, 1e-15, close_limit, ctx)
        output('ref. coordinates: %f s' % timer.stop(), verbose=verbose)

    else:
        cells = cache.cells
        status = cache.status

    return ref_coors, cells, status

def get_potential_cells(coors, cmesh, centroids=None, extrapolate=True):
    """
    Get cells that potentially contain points with the given physical
    coordinates.

    Parameters
    ----------
    coors : array
        The physical coordinates.
    cmesh : CMesh instance
        The cmesh defining the cells.
    centroids : array, optional
        The centroids of the cells.
    extrapolate : bool
        If True, even the points that are surely outside of the
        cmesh are considered and assigned potential cells.

    Returns
    -------
    potential_cells : array
        The indices of the cells that potentially contain the points.
    offsets : array
        The offsets into `potential_cells` for each point: a point ``ip`` is
        potentially in cells ``potential_cells[offsets[ip]:offsets[ip+1]]``.
    """
    from scipy.spatial import cKDTree as KDTree

    if centroids is None:
        centroids = cmesh.get_centroids(cmesh.tdim)

    kdtree = KDTree(coors)

    conn = cmesh.get_cell_conn()
    cc = conn.indices.reshape(cmesh.n_el, -1)
    cell_coors = cmesh.coors[cc]

    rays = cell_coors - centroids[:, None]
    radii = nm.linalg.norm(rays, ord=nm.inf, axis=2).max(axis=1)

    potential_cells = [[]] * coors.shape[0]
    for ic, centroid in enumerate(centroids):
        ips = kdtree.query_ball_point(centroid, radii[ic], p=nm.inf)
        if len(ips):
            for ip in ips:
                if not len(potential_cells[ip]):
                    potential_cells[ip] = []

                potential_cells[ip].append(ic)

    lens = nm.array([0] + [len(ii) for ii in potential_cells], dtype=nm.int32)

    if extrapolate:
        # Deal with the points outside of the field domain - insert elements
        # incident to the closest mesh vertex.
        iin = nm.where(lens[1:] == 0)[0]
        if len(iin):
            kdtree = KDTree(cmesh.coors)
            ics = kdtree.query(coors[iin])[1]
            cmesh.setup_connectivity(0, cmesh.tdim)
            conn = cmesh.get_conn(0, cmesh.tdim)

            oo = conn.offsets
            for ii, ip in enumerate(iin):
                ik = ics[ii]
                potential_cells[ip] = conn.indices[oo[ik]:oo[ik+1]]
                lens[ip+1] = len(potential_cells[ip])

    offsets = nm.cumsum(lens, dtype=nm.int32)
    potential_cells = nm.concatenate(potential_cells).astype(nm.int32)

    return potential_cells, offsets

def get_ref_coors_general(field, coors, close_limit=0.1, get_cells_fun=None,
                          cache=None, verbose=False):
    """
    Get reference element coordinates and elements corresponding to given
    physical coordinates.

    Parameters
    ----------
    field : Field instance
        The field defining the approximation.
    coors : array
        The physical coordinates.
    close_limit : float, optional
        The maximum limit distance of a point from the closest
        element allowed for extrapolation.
    get_cells_fun : callable, optional
        If given, a function with signature ``get_cells_fun(coors, cmesh,
        **kwargs)`` returning cells and offsets that potentially contain points
        with the coordinates `coors`. When not given,
        :func:`get_potential_cells()` is used.
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
        non-convergence of the Newton iteration in tensor product cells. If
        close_limit is 0, then status 5 indicates points outside of the field
        domain that had no potential cells.
    """
    timer = Timer()

    ref_coors = get_default_attr(cache, 'ref_coors', None)
    if ref_coors is None:
        extrapolate = close_limit > 0.0

        get = get_potential_cells if get_cells_fun is None else get_cells_fun

        ref_coors = nm.empty_like(coors)
        cells = nm.empty((coors.shape[0],), dtype=nm.int32)
        status = nm.empty((coors.shape[0],), dtype=nm.int32)

        cmesh = get_default_attr(cache, 'cmesh', None)
        if cmesh is None:
            timer.start()
            mesh = field.create_mesh(extra_nodes=False)
            cmesh = mesh.cmesh

            if get_cells_fun is None:
                centroids = cmesh.get_centroids(cmesh.tdim)

            else:
                centroids = None

            output('cmesh setup: %f s' % timer.stop(), verbose=verbose)

        else:
            centroids = cache.centroids

        timer.start()
        potential_cells, offsets = get(coors, cmesh, centroids=centroids,
                                       extrapolate=extrapolate)
        output('potential cells: %f s' % timer.stop(), verbose=verbose)

        coors = nm.ascontiguousarray(coors)
        ctx = field.create_basis_context()

        eval_cmesh = get_default_attr(cache, 'eval_cmesh', None)
        if eval_cmesh is None:
            timer.start()
            mesh = field.create_eval_mesh()
            if mesh is None:
                eval_cmesh = cmesh

            else:
                eval_cmesh = mesh.cmesh

            output('eval_cmesh setup: %f s'
                   % timer.stop(), verbose=verbose)

        timer.start()

        crc.find_ref_coors(ref_coors, cells, status, coors, eval_cmesh,
                           potential_cells, offsets, extrapolate,
                           1e-15, close_limit, ctx)
        if extrapolate:
            assert_(nm.all(status < 5))

        output('ref. coordinates: %f s' % timer.stop(), verbose=verbose)

    else:
        cells = cache.cells
        status = cache.status

    return ref_coors, cells, status

def get_ref_coors(field, coors, strategy='general', close_limit=0.1,
                  get_cells_fun=None, cache=None, verbose=False):
    """
    Get reference element coordinates and elements corresponding to given
    physical coordinates.

    Parameters
    ----------
    field : Field instance
        The field defining the approximation.
    coors : array
        The physical coordinates.
    strategy : {'general', 'convex'}, optional
        The strategy for finding the elements that contain the coordinates. For
        convex meshes, the 'convex' strategy might be faster than the 'general'
        one.
    close_limit : float, optional
        The maximum limit distance of a point from the closest
        element allowed for extrapolation.
    get_cells_fun : callable, optional
        If given, a function with signature ``get_cells_fun(coors, cmesh,
        **kwargs)`` returning cells and offsets that potentially contain points
        with the coordinates `coors`. Applicable only when `strategy` is
        'general'. When not given, :func:`get_potential_cells()` is used.
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
        non-convergence of the Newton iteration in tensor product cells. If
        close_limit is 0, then for the 'general' strategy the status 5
        indicates points outside of the field domain that had no potential
        cells.
    """
    if strategy == 'general':
        return get_ref_coors_general(field, coors, close_limit=close_limit,
                                     get_cells_fun=get_cells_fun,
                                     cache=cache, verbose=verbose)

    elif strategy == 'convex':
        return get_ref_coors_convex(field, coors, close_limit=close_limit,
                                    cache=cache, verbose=verbose)

    else:
        raise ValueError('unsupported strategy! (%s)' % strategy)
