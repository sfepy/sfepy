from sfepy.base.base import *
import sfepy.base.la as la
from sfepy.fem.integrals import Integral
from extmods.geometry import SurfaceGeometry

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
        ap = field.aps.aps_per_group[ig]
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

        sg = SurfaceGeometry( n_fa, n_fp, dim, n_fp )
        sg.describe( field.aps.coors, econn, bf_sg, weights )

        e_normals = sg.variable( 0 ).squeeze()

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

def extend_cell_data( data, domain, rname, val = None ):
    """Extend cell data defined in a region rname to the whole domain using the
    value val, or the smallest value in data if val is None."""
    n_el = domain.shape.n_el
    if data.shape[0] == n_el: return data

    if val is None:
        if data.shape[2] > 1: # Vector.
            val = nm.amin( nm.abs( data ) )
        else: # Scalar.
            val = nm.amin( data )

    edata = nm.empty( (n_el,) + data.shape[1:], dtype = nm.float64 )
    edata.fill( val )

    region = domain.regions[rname]
    offs = region.get_cell_offsets()
    eoffs = domain.get_cell_offsets()
##     print offs
##     print eoffs
##     print domain.mat_ids_to_i_gs
##     pause()

    for group in domain.iter_groups():
        ig = group.ig
        ii = eoffs[ig]
        if ig in region.igs:
            n_cell = region.shape[ig].n_cell
            ir = offs[ig]
            edata[ii+region.cells[ig]] = data[ir:ir+n_cell]
    return edata
