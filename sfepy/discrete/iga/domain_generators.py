"""
IGA domain generators.
"""
from __future__ import absolute_import
import numpy as nm

from sfepy.base.base import assert_, output, Struct
import sfepy.discrete.iga as iga
from sfepy.discrete.iga.domain import NurbsPatch
from six.moves import range

def create_from_igakit(inurbs, verbose=False):
    """
    Create :class:`IGDomain <sfepy.discrete.iga.domain.IGDomain>` data from a
    given igakit NURBS object.

    Parameters
    ----------
    inurbs : igakit.nurbs.NURBS instance
        The igakit NURBS object.

    Returns
    -------
    nurbs : NurbsPatch instance
        The NURBS data. The igakit NURBS object is stored as `nurbs` attribute.
    bmesh : Struct instance
        The Bezier mesh data.
    regions : dict
        The patch surface regions.
    """
    # Compute Bezier extraction data.
    output('computing Bezier mesh...', verbose=verbose)
    cs = iga.compute_bezier_extraction(inurbs.knots, inurbs.degree)
    n_els = [len(ii) for ii in cs]
    conn, bconn = iga.create_connectivity(n_els, inurbs.knots, inurbs.degree)

    ccs = iga.combine_bezier_extraction(cs)

    dim = len(inurbs.shape)
    cps = inurbs.points[..., :dim].copy()
    cps = cps.reshape((-1, dim))
    bcps, bweights = iga.compute_bezier_control(cps, inurbs.weights.ravel(),
                                                ccs, conn, bconn)

    nurbs = NurbsPatch(inurbs.knots, inurbs.degree, cps,
                       inurbs.weights.ravel(), cs, conn)
    nurbs.nurbs = inurbs
    bmesh = Struct(name='bmesh', cps=bcps, weights=bweights, conn=bconn)
    output('...done', verbose=verbose)

    output('defining regions...', verbose=verbose)
    regions = iga.get_patch_box_regions(n_els, inurbs.degree)
    output('...done', verbose=verbose)

    return nurbs, bmesh, regions

def gen_patch_block_domain(dims, shape, centre, degrees, continuity=None,
                           cp_mode='greville', name='block', verbose=True):
    """
    Generate a single IGA patch block in 2D or 3D of given degrees and
    continuity using igakit.

    Parameters
    ----------
    dims : array of D floats
        Dimensions of the block.
    shape : array of D ints
        Numbers of unique knot values along each axis.
    centre : array of D floats
        Centre of the block.
    degrees : array of D floats
        NURBS degrees along each axis.
    continuity : array of D ints, optional
        NURBS continuity along each axis. If None, `degrees-1` is used.
    cp_mode : 'greville' or 'uniform'
        The control points mode. The default 'greville' results in a uniform
        Bezier mesh, while the 'uniform' mode results in a uniform grid of
        control points a finer Bezier mesh inside the block and a coarser
        Bezier mesh near the block boundary.
    name : string
        Domain name.
    verbose : bool
        If True, report progress of the domain generation.

    Returns
    -------
    nurbs : NurbsPatch instance
        The NURBS data. The igakit NURBS object is stored as `nurbs` attribute.
    bmesh : Struct instance
        The Bezier mesh data.
    regions : dict
        The patch surface regions.
    """
    import igakit.cad as cad

    dims = nm.asarray(dims, dtype=nm.float64)
    shape = nm.asarray(shape, dtype=nm.int32)
    centre = nm.asarray(centre, dtype=nm.float64)
    degrees = nm.asarray(degrees, dtype=nm.int32)

    assert_(cp_mode in ['greville', 'uniform'],
            "cp_mode has to be 'greville' or 'uniform'")

    if continuity is None:
        continuity = degrees - 1

    else:
        continuity = nm.asarray(continuity, dtype=nm.int32)

    dim = len(shape)

    output('generating NURBS...', verbose=verbose)

    dd = centre - 0.5 * dims
    block = cad.grid(shape - 1, degree=degrees, continuity=continuity)

    for ia in range(dim):
        block.scale(dims[ia], ia)

    for ia in range(dim):
        block.translate(dd[ia], ia)

    if cp_mode == 'uniform':
        # Force uniform control points. This prevents coarser resolution inside
        # the block.
        shape = nm.asarray(block.points.shape[:-1])
        n_nod = nm.prod(shape)
        x0 = centre - 0.5 * dims
        dd = dims / (shape - 1)

        ngrid = nm.mgrid[[slice(ii) for ii in shape]]
        ngrid.shape = (dim, n_nod)

        coors = x0 + ngrid.T * dd
        coors.shape = shape.tolist() + [dim]

        block.array[..., :dim] = coors

    output('...done', verbose=verbose)

    nurbs, bmesh, regions = create_from_igakit(block, verbose=verbose)

    return nurbs, bmesh, regions
