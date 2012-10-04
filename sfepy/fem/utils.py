import numpy as nm

import sfepy.linalg as la
from extmods.mappings import CSurfaceMapping

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
    """Nodal normals are computed by simple averaging of element normals of
    elements every node is contained in. """
    dim = field.shape[0]

    region.select_cells_of_surface()

    normals = nm.zeros( (nodes.shape[0], dim),
                        dtype = nm.float64 )
    mask = nm.zeros( (nodes.max()+1,), dtype = nm.int32 )
    imap = nm.empty_like( mask )
    imap.fill( nodes.shape[0] ) # out-of-range index for normals.
    imap[nodes] = nm.arange( nodes.shape[0], dtype = nm.int32 )
    
    for ig, fis in region.fis.iteritems():
        ap = field.aps[ig]
        n_fa = fis.shape[0]
        n_fp = ap.efaces.shape[1]
        face_type = 's%d' % n_fp

        faces = ap.efaces[fis[:,1]]
        ee = ap.econn[fis[:,0]]
        econn = nm.empty( faces.shape, dtype = nm.int32 )
        for ir, face in enumerate( faces ):
            econn[ir] = ee[ir,face]
        mask[econn] += 1
        # Unit normals -> weights = ones.
        ps = ap.interp.poly_spaces[face_type]
        weights = nm.ones((n_fp,), dtype=nm.float64)

        coors = ps.node_coors
        bf_sg = ps.eval_base(coors, diff=True)

        cmap = CSurfaceMapping(n_fa, n_fp, dim, n_fp)
        cmap.describe(field.get_coor(), econn, bf_sg, weights)

        e_normals = cmap.normal.squeeze()

        # normals[imap[econn]] += e_normals
        im = imap[econn]
        for ii, en in enumerate( e_normals ):
            normals[im[ii]] += en

    # All nodes must have a normal.
    if not nm.all( mask[nodes] > 0 ):
        raise ValueError( 'region %s has not complete faces!' % region.name )

    normals /= la.norm_l2_along_axis( normals )[:,nm.newaxis]

    if return_imap:
        return normals, imap

    else:
        return normals

def extend_cell_data(data, domain, rname, val=None, is_surface=False):
    """
    Extend cell data defined in a region to the whole domain.

    Parameters
    ----------
    data : array
        The data defined in the region.
    domain : Domain instance
        The FE domain.
    rname : str
        The region name.
    val : float, optional
        The value for filling cells not covered by the region. If not given,
        the smallest value in data is used.
    is_surface : bool
        If True, the data are defined on a surface region. In that case the
        values are averaged into the cells containing the region surface faces.

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

    edata = nm.empty((n_el,) + data.shape[1:], dtype = nm.float64)
    edata.fill(val)

    region = domain.regions[rname]
    eoffs = domain.get_cell_offsets()

    if not is_surface:
        offs = region.get_cell_offsets()
        for group in domain.iter_groups():
            ig = group.ig
            ii = eoffs[ig]
            if ig in region.igs:
                n_cell = region.shape[ig].n_cell
                ir = offs[ig]
                edata[ii+region.cells[ig]] = data[ir:ir+n_cell]

    else:
        for group in domain.iter_groups():
            ig = group.ig
            ii = eoffs[ig]
            if ig in region.igs:
                cells = region.cells[ig]
                ucells = nm.unique(cells)

                avg = nm.bincount(cells, minlength=group.shape.n_el)[ucells]
                for ic in xrange(data.shape[2]):
                    evals = nm.bincount(cells, weights=data[:, 0, ic, 0],
                                        minlength=group.shape.n_el)[ucells]

                    edata[ii+ucells, 0, ic, 0] = evals / avg

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
    from sfepy.fem import Mesh, Domain

    if level > 0:
        mesh = Mesh.from_file(filename)
        domain = Domain(mesh.name, mesh)
        for ii in range(level):
            output('refine %d...' % ii)
            domain = domain.refine()
            output('... %d nodes %d elements'
                   % (domain.shape.n_nod, domain.shape.n_el))

        suffix = os.path.splitext(filename)[1]
        filename = domain.name + suffix

        domain.mesh.write(filename, io='auto')

    return filename

