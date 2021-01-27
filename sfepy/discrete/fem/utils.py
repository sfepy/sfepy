from __future__ import absolute_import
import numpy as nm

import sfepy.linalg as la
from sfepy.discrete.integrals import Integral
from sfepy.discrete import PolySpace
from six.moves import range

def prepare_remap(indices, n_full):
    """
    Prepare vector for remapping range `[0, n_full]` to its subset given
    by `indices`.
    """
    remap = nm.empty((n_full,), dtype=nm.int32)
    remap.fill(-1)
    remap[indices] = nm.arange(indices.shape[0], dtype=nm.int32)

    return remap

def invert_remap(remap):
    """
    Return the inverse of `remap`, i.e. a mapping from a sub-range
    indices to a full range, see :func:`prepare_remap()`.
    """
    if remap is not None:
        inverse = nm.where(remap >= 0)[0].astype(nm.int32)

    else:
        inverse = None

    return inverse

def prepare_translate(old_indices, new_indices):
    """
    Prepare vector for translating `old_indices` to `new_indices`.

    Returns
    -------
    translate : array
        The translation vector. Then `new_ar = translate[old_ar]`.
    """
    old_indices = nm.asarray(old_indices)
    new_indices = nm.asarray(new_indices)

    translate = nm.zeros(old_indices.max() + 1, dtype=new_indices.dtype)
    translate[old_indices] = new_indices

    return translate

def compute_nodal_normals(nodes, region, field, return_imap=False):
    """
    Nodal normals are computed by simple averaging of element normals of
    elements every node is contained in.
    """
    dim = region.dim

    field.domain.create_surface_group(region)
    field.setup_surface_data(region)

    # Custom integral with quadrature points in nodes.
    ps = PolySpace.any_from_args('', field.gel.surface_facet,
                                 field.approx_order)
    qp_coors = ps.node_coors
    # Unit normals -> weights = ones.
    qp_weights = nm.ones(qp_coors.shape[0], dtype=nm.float64)

    integral = Integral('aux', coors=qp_coors, weights=qp_weights)

    normals = nm.zeros((nodes.shape[0], dim), dtype=nm.float64)
    mask = nm.zeros((nodes.max() + 1,), dtype=nm.int32)
    imap = nm.empty_like(mask)
    imap.fill(nodes.shape[0]) # out-of-range index for normals.
    imap[nodes] = nm.arange(nodes.shape[0], dtype=nm.int32)

    cmap, _ = field.get_mapping(region, integral, 'surface')
    e_normals = cmap.normal[..., 0]

    sd = field.surface_data[region.name]
    econn = sd.get_connectivity()
    mask[econn] += 1

    # normals[imap[econn]] += e_normals
    im = imap[econn]
    for ii, en in enumerate(e_normals):
        normals[im[ii]] += en

    # All nodes must have a normal.
    if not nm.all(mask[nodes] > 0):
        raise ValueError('region %s has not complete faces!' % region.name)

    norm = la.norm_l2_along_axis(normals)[:, nm.newaxis]
    if (norm < 1e-15).any():
        raise ValueError('zero nodal normal! (a node in volume?)')
    normals /= norm

    if return_imap:
        return normals, imap

    else:
        return normals

def _get_edge_path(graph, seed, mask, cycle=False):
    """
    Get a path in an edge graph starting with seed. The mask is incremented by
    one at positions of the path vertices.
    """
    if mask[seed]:
        return []

    path = [seed]
    mask[seed] = 1

    row = graph[seed].indices
    nv = len(row)
    while nv:
        if nv == 2:
            if mask[row[0]]:
                if mask[row[1]]:
                    if cycle:
                        path.append(seed)
                    break

                else:
                    vert = row[1]

            else:
                vert = row[0]

        elif mask[row[0]]:
            break

        else:
            vert = row[0]

        path.append(vert)
        mask[vert] = 1

        row = graph[vert].indices
        nv = len(row)

    path = nm.array(path, dtype=nm.int32)

    return path

def get_edge_paths(graph, mask):
    """
    Get all edge paths in a graph with non-masked vertices. The mask is
    updated.
    """
    nodes = nm.unique(graph.indices)
    npv = nm.diff(graph.indptr)

    if npv.max() > 2:
        raise ValueError('more than 2 edges sharing a vertex!')

    seeds = nm.where(npv == 1)[0]

    # 1. get paths.
    paths = []
    for seed in seeds:
        path = _get_edge_path(graph, seed, mask)
        if len(path):
            paths.append(path)

    # 2. get possible remaing cycles.
    while 1:
        ii = nm.where(mask[nodes] == 0)[0]
        if not len(ii):
            break

        path = _get_edge_path(graph, nodes[ii[0]], mask, cycle=True)
        if len(path):
            paths.append(path)

    return paths

def compute_nodal_edge_dirs(nodes, region, field, return_imap=False):
    """
    Nodal edge directions are computed by simple averaging of direction vectors
    of edges a node is contained in. Edges are assumed to be straight and a
    node must be on a single edge (a border node) or shared by exactly two
    edges.
    """
    coors = region.domain.mesh.coors
    dim = coors.shape[1]

    graph = region.get_edge_graph()

    imap = prepare_remap(nodes, nodes.max() + 1)
    mask = nm.zeros_like(imap)

    try:
        paths = get_edge_paths(graph, mask)

    except ValueError:
        raise ValueError('more than 2 edges sharing a vertex in region %s!'
                         % region.name)

    # All nodes must have an edge direction.
    if not nm.all(mask[nodes]):
        raise ValueError('region %s has not complete edges!' % region.name)

    edge_dirs = nm.zeros((nodes.shape[0], dim), dtype=nm.float64)
    for path in paths:
        pcoors = coors[path]

        edirs = nm.diff(pcoors, axis=0)
        la.normalize_vectors(edirs, eps=1e-12)

        im = imap[nm.c_[path[:-1], path[1:]]]
        for ii, edir in enumerate(edirs):
            edge_dirs[im[ii]] += edir

    la.normalize_vectors(edge_dirs, eps=1e-12)

    if return_imap:
        return edge_dirs, imap

    else:
        return edge_dirs

def get_min_value(dofs):
    """
    Get a reasonable minimal value of DOFs suitable for extending over a
    whole domain.
    """
    if dofs.shape[1] > 1: # Vector.
        val = 0.0

    else: # Scalar.
        val = dofs.min()

    return val

def extend_cell_data(data, domain, rname, val=None, is_surface=False,
                     average_surface=True):
    """
    Extend cell data defined in a region to the whole domain.

    Parameters
    ----------
    data : array
        The data defined in the region.
    domain : FEDomain instance
        The FE domain.
    rname : str
        The region name.
    val : float, optional
        The value for filling cells not covered by the region. If not given,
        the smallest value in data is used.
    is_surface : bool
        If True, the data are defined on a surface region. In that case the
        values are averaged or summed into the cells containing the region
        surface faces (a cell can have several faces of the surface), see
        `average_surface`.
    average_surface : bool
        If True, the data defined on a surface region are averaged, otherwise
        the data are summed.

    Returns
    -------
    edata : array
        The data extended to all domain elements.
    """
    n_el = domain.shape.n_el
    if data.shape[0] == n_el: return data

    if val is None:
        if data.shape[2] > 1: # Vector.
            val = nm.amin(nm.abs(data))
        else: # Scalar.
            val = nm.amin(data)

    edata = nm.empty((n_el,) + data.shape[1:], dtype=data.dtype)
    edata.fill(val)

    region = domain.regions[rname]

    if not is_surface:
        edata[region.get_cells()] = data

    else:
        cells = region.get_cells(true_cells_only=False)
        ucells = nm.unique(cells)

        if len(cells) != len(region.facets):
            raise ValueError('region %s has an inner face!'
                             % region.name)

        if average_surface:
            avg = nm.bincount(cells, minlength=n_el)[ucells]

        else:
            avg = 1.0

        for ic in range(data.shape[2]):
            if nm.isrealobj(data):
                evals = nm.bincount(cells, weights=data[:, 0, ic, 0],
                                    minlength=n_el)[ucells]

            else:
                evals = (nm.bincount(cells, weights=data[:, 0, ic, 0].real,
                                     minlength=n_el)[ucells]
                         + 1j *
                         nm.bincount(cells, weights=data[:, 0, ic, 0].imag,
                                     minlength=n_el)[ucells])

            edata[ucells, 0, ic, 0] = evals / avg

    return edata

def refine_mesh(filename, level):
    """
    Uniformly refine `level`-times a mesh given by `filename`.

    The refined mesh is saved to a file with name constructed from base
    name of `filename` and `level`-times appended `'_r'` suffix.

    Parameters
    ----------
    filename : str
        The mesh file name.
    level : int
        The refinement level.
    """
    import os
    from sfepy.base.base import output
    from sfepy.discrete.fem import Mesh, FEDomain

    if level > 0:
        mesh = Mesh.from_file(filename)
        domain = FEDomain(mesh.name, mesh)
        for ii in range(level):
            output('refine %d...' % ii)
            domain = domain.refine()
            output('... %d nodes %d elements'
                   % (domain.shape.n_nod, domain.shape.n_el))

        suffix = os.path.splitext(filename)[1]
        filename = domain.name + suffix

        domain.mesh.write(filename, io='auto')

    return filename
