from __future__ import print_function
from __future__ import absolute_import
import numpy as nm
import sys
from six.moves import range
sys.path.append('.')

from sfepy.base.base import output, assert_
from sfepy.base.ioutils import ensure_path
from sfepy.linalg import cycle
from sfepy.discrete.fem.mesh import Mesh
from sfepy.mesh.mesh_tools import elems_q2t

def get_tensor_product_conn(shape):
    """
    Generate vertex connectivity for cells of a tensor-product mesh of the
    given shape.

    Parameters
    ----------
    shape : array of 2 or 3 ints
        Shape (counts of nodes in x, y, z) of the mesh.

    Returns
    -------
    conn : array
        The vertex connectivity array.
    desc : str
        The cell kind.
    """
    shape = nm.asarray(shape)
    dim = len(shape)
    assert_(1 <= dim <= 3)

    n_nod = nm.prod(shape)
    n_el = nm.prod(shape - 1)

    grid = nm.arange(n_nod, dtype=nm.int32)
    grid.shape = shape

    if dim == 1:
        conn = nm.zeros((n_el, 2), dtype=nm.int32)
        conn[:, 0] = grid[:-1]
        conn[:, 1] = grid[1:]
        desc = '1_2'

    elif dim == 2:
        conn = nm.zeros((n_el, 4), dtype=nm.int32)
        conn[:, 0] = grid[:-1, :-1].flat
        conn[:, 1] = grid[1:, :-1].flat
        conn[:, 2] = grid[1:, 1:].flat
        conn[:, 3] = grid[:-1, 1:].flat
        desc = '2_4'

    else:
        conn = nm.zeros((n_el, 8), dtype=nm.int32)
        conn[:, 0] = grid[:-1, :-1, :-1].flat
        conn[:, 1] = grid[1:, :-1, :-1].flat
        conn[:, 2] = grid[1:, 1:, :-1].flat
        conn[:, 3] = grid[:-1, 1:, :-1].flat
        conn[:, 4] = grid[:-1, :-1, 1:].flat
        conn[:, 5] = grid[1:, :-1, 1:].flat
        conn[:, 6] = grid[1:, 1:, 1:].flat
        conn[:, 7] = grid[:-1, 1:, 1:].flat
        desc = '3_8'

    return conn, desc

def gen_block_mesh(dims, shape, centre, mat_id=0, name='block',
                   coors=None, verbose=True):
    """
    Generate a 2D or 3D block mesh. The dimension is determined by the
    lenght of the shape argument.

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
    mesh : Mesh instance
    """
    dims = nm.asarray(dims, dtype=nm.float64)
    shape = nm.asarray(shape, dtype=nm.int32)
    centre = nm.asarray(centre, dtype=nm.float64)

    dim = shape.shape[0]

    centre = centre[:dim]
    dims = dims[:dim]

    n_nod = nm.prod(shape)
    output('generating %d vertices...' % n_nod, verbose=verbose)

    x0 = centre - 0.5 * dims
    dd = dims / (shape - 1)

    ngrid = nm.mgrid[[slice(ii) for ii in shape]]
    ngrid.shape = (dim, n_nod)

    if coors is None:
        coors = x0 + ngrid.T * dd
    output('...done', verbose=verbose)

    n_el = nm.prod(shape - 1)
    output('generating %d cells...' % n_el, verbose=verbose)

    mat_ids = nm.empty((n_el,), dtype=nm.int32)
    mat_ids.fill(mat_id)

    conn, desc = get_tensor_product_conn(shape)
    output('...done', verbose=verbose)

    mesh = Mesh.from_data(name, coors, None, [conn], [mat_ids], [desc])
    return mesh

def gen_cylinder_mesh(dims, shape, centre, axis='x', force_hollow=False,
                      is_open=False, open_angle=0.0, non_uniform=False,
                      name='cylinder', verbose=True):
    """
    Generate a cylindrical mesh along an axis. Its cross-section can be
    ellipsoidal.

    Parameters
    ----------
    dims : array of 5 floats
        Dimensions of the cylinder: inner surface semi-axes a1, b1, outer
        surface semi-axes a2, b2, length.
    shape : array of 3 ints
        Shape (counts of nodes in radial, circumferential and longitudinal
        directions) of the cylinder mesh.
    centre : array of 3 floats
        Centre of the cylinder.
    axis: one of 'x', 'y', 'z'
        The axis of the cylinder.
    force_hollow : boolean
        Force hollow mesh even if inner radii a1 = b1 = 0.
    is_open : boolean
        Generate an open cylinder segment.
    open_angle : float
        Opening angle in radians.
    non_uniform : boolean
        If True, space the mesh nodes in radial direction so that the element
        volumes are (approximately) the same, making thus the elements towards
        the outer surface thinner.
    name : string
        Mesh name.
    verbose : bool
        If True, show progress of the mesh generation.

    Returns
    -------
    mesh : Mesh instance
    """
    dims = nm.asarray(dims, dtype=nm.float64)
    shape = nm.asarray(shape, dtype=nm.int32)
    centre = nm.asarray(centre, dtype=nm.float64)

    a1, b1, a2, b2, length = dims
    nr, nfi, nl = shape
    origin = centre - nm.array([0.5 * length, 0.0, 0.0])

    dfi = 2.0 * (nm.pi - open_angle) / nfi
    if is_open:
        nnfi = nfi + 1
    else:
        nnfi = nfi

    is_hollow = force_hollow or not (max(abs(a1), abs(b1)) < 1e-15)

    if is_hollow:
        mr = 0
    else:
        mr = (nnfi - 1) * nl

    grid = nm.zeros((nr, nnfi, nl), dtype=nm.int32)

    n_nod = nr * nnfi * nl - mr
    coors = nm.zeros((n_nod, 3), dtype=nm.float64)

    angles = nm.linspace(open_angle, open_angle+(nfi)*dfi, nfi+1)
    xs = nm.linspace(0.0, length, nl)
    if non_uniform:
        ras = nm.zeros((nr,), dtype=nm.float64)
        rbs = nm.zeros_like(ras)
        advol = (a2**2 - a1**2) / (nr - 1)
        bdvol = (b2**2 - b1**2) / (nr - 1)
        ras[0], rbs[0] = a1, b1
        for ii in range(1, nr):
            ras[ii] = nm.sqrt(advol + ras[ii-1]**2)
            rbs[ii] = nm.sqrt(bdvol + rbs[ii-1]**2)
    else:
        ras = nm.linspace(a1, a2, nr)
        rbs = nm.linspace(b1, b2, nr)

    # This is 3D only...
    output('generating %d vertices...' % n_nod, verbose=verbose)
    ii = 0
    for ix in range(nr):
        a, b = ras[ix], rbs[ix]
        for iy, fi in enumerate(angles[:nnfi]):
            for iz, x in enumerate(xs):
                grid[ix,iy,iz] = ii
                coors[ii] = origin + [x, a * nm.cos(fi), b * nm.sin(fi)]
                ii += 1

                if not is_hollow and (ix == 0):
                    if iy > 0:
                        grid[ix,iy,iz] = grid[ix,0,iz]
                        ii -= 1
    assert_(ii == n_nod)
    output('...done', verbose=verbose)

    n_el = (nr - 1) * nfi * (nl - 1)
    conn = nm.zeros((n_el, 8), dtype=nm.int32)

    output('generating %d cells...' % n_el, verbose=verbose)
    ii = 0
    for (ix, iy, iz) in cycle([nr-1, nnfi, nl-1]):
        if iy < (nnfi - 1):
            conn[ii,:] = [grid[ix  ,iy  ,iz  ], grid[ix+1,iy  ,iz  ],
                          grid[ix+1,iy+1,iz  ], grid[ix  ,iy+1,iz  ],
                          grid[ix  ,iy  ,iz+1], grid[ix+1,iy  ,iz+1],
                          grid[ix+1,iy+1,iz+1], grid[ix  ,iy+1,iz+1]]
            ii += 1
        elif not is_open:
            conn[ii,:] = [grid[ix  ,iy  ,iz  ], grid[ix+1,iy  ,iz  ],
                          grid[ix+1,0,iz  ], grid[ix  ,0,iz  ],
                          grid[ix  ,iy  ,iz+1], grid[ix+1,iy  ,iz+1],
                          grid[ix+1,0,iz+1], grid[ix  ,0,iz+1]]
            ii += 1

    mat_id = nm.zeros((n_el,), dtype = nm.int32)
    desc = '3_8'

    assert_(n_nod == (conn.max() + 1))
    output('...done', verbose=verbose)

    if axis == 'z':
        coors = coors[:,[1,2,0]]
    elif axis == 'y':
        coors = coors[:,[2,0,1]]

    mesh = Mesh.from_data(name, coors, None, [conn], [mat_id], [desc])
    return mesh

def _spread_along_axis(axis, coors, tangents, grading_fun):
    """
    Spread the coordinates along the given axis using the grading function, and
    the tangents in the other two directions.
    """
    oo = list(set([0, 1, 2]).difference([axis]))
    c0, c1, c2 = coors[:, axis], coors[:, oo[0]], coors[:, oo[1]]

    out = nm.empty_like(coors)

    mi, ma = c0.min(), c0.max()
    nc0 = (c0 - mi) / (ma - mi)
    out[:, axis] = oc0 = grading_fun(nc0) * (ma - mi) + mi

    nc = oc0 - oc0.min()

    mi, ma = c1.min(), c1.max()
    n1 = 2 * (c1 - mi) / (ma - mi) - 1
    out[:, oo[0]] = c1 + n1 * nc * tangents[oo[0]]

    mi, ma = c2.min(), c2.max()
    n2 = 2 * (c2 - mi) / (ma - mi) - 1
    out[:, oo[1]] = c2 + n2 * nc * tangents[oo[1]]

    return out

def _get_extension_side(side, grading_fun, mat_id,
                        b_dims, b_shape, e_dims, e_shape, centre):
    """
    Get a mesh extending the given side of a block mesh.
    """
    # Pure extension dimensions.
    pe_dims = 0.5 * (e_dims - b_dims)
    coff = 0.5 * (b_dims + pe_dims)
    cc = centre + coff * nm.eye(3)[side]

    if side == 0: # x axis.
        dims = [pe_dims[0], b_dims[1], b_dims[2]]
        shape = [e_shape, b_shape[1], b_shape[2]]
        tangents = [0, pe_dims[1] / pe_dims[0], pe_dims[2] / pe_dims[0]]

    elif side == 1: # y axis.
        dims = [b_dims[0], pe_dims[1], b_dims[2]]
        shape = [b_shape[0], e_shape, b_shape[2]]
        tangents = [pe_dims[0] / pe_dims[1], 0, pe_dims[2] / pe_dims[1]]

    elif side == 2: # z axis.
        dims = [b_dims[0], b_dims[1], pe_dims[2]]
        shape = [b_shape[0], b_shape[1], e_shape]
        tangents = [pe_dims[0] / pe_dims[2], pe_dims[1] / pe_dims[2], 0]

    e_mesh = gen_block_mesh(dims, shape, cc, mat_id=mat_id, verbose=False)
    e_mesh.coors[:] = _spread_along_axis(side, e_mesh.coors, tangents,
                                         grading_fun)

    return e_mesh, shape

def gen_extended_block_mesh(b_dims, b_shape, e_dims, e_shape, centre,
                            grading_fun=None, name=None):
    """
    Generate a 3D mesh with a central block and (coarse) extending side meshes.

    The resulting mesh is again a block. Each of the components has a different
    material id.

    Parameters
    ----------
    b_dims : array of 3 floats
        The dimensions of the central block.
    b_shape : array of 3 ints
        The shape (counts of nodes in x, y, z) of the central block mesh.
    e_dims : array of 3 floats
        The dimensions of the complete block (central block + extensions).
    e_shape : int
        The count of nodes of extending blocks in the direction from the
        central block.
    centre : array of 3 floats
        The centre of the mesh.
    grading_fun : callable, optional
        A function of :math:`x \in [0, 1]` that can be used to shift nodes in
        the extension axis directions to allow smooth grading of element sizes
        from the centre. The default function is :math:`x**p` with :math:`p`
        determined so that the element sizes next to the central block have the
        size of the shortest edge of the central block.
    name : string, optional
        The mesh name.

    Returns
    -------
    mesh : Mesh instance
    """
    b_dims = nm.asarray(b_dims, dtype=nm.float64)
    b_shape = nm.asarray(b_shape, dtype=nm.int32)
    e_dims = nm.asarray(e_dims, dtype=nm.float64)
    centre = nm.asarray(centre, dtype=nm.float64)

    # Pure extension dimensions.
    pe_dims = 0.5 * (e_dims - b_dims)
    # Central block element sizes.
    dd = (b_dims / (b_shape - 1))
    # The "first x" going to grading_fun.
    nc = 1.0 / (e_shape - 1)
    # Grading power and function.
    power = nm.log(dd.min() / pe_dims.min()) / nm.log(nc)
    grading_fun = (lambda x: x**power) if grading_fun is None else grading_fun

    # Central block mesh.
    b_mesh = gen_block_mesh(b_dims, b_shape, centre, mat_id=0, verbose=False)

    # 'x' extension.
    e_mesh, xs = _get_extension_side(0, grading_fun, 10,
                                     b_dims, b_shape, e_dims, e_shape, centre)
    mesh = b_mesh + e_mesh

    # Mirror by 'x'.
    e_mesh.coors[:, 0] = (2 * centre[0]) - e_mesh.coors[:, 0]
    e_mesh.cmesh.cell_groups.fill(11)
    mesh = mesh + e_mesh

    # 'y' extension.
    e_mesh, ys = _get_extension_side(1, grading_fun, 20,
                                     b_dims, b_shape, e_dims, e_shape, centre)
    mesh = mesh + e_mesh

    # Mirror by 'y'.
    e_mesh.coors[:, 1] = (2 * centre[1]) - e_mesh.coors[:, 1]
    e_mesh.cmesh.cell_groups.fill(21)
    mesh = mesh + e_mesh

    # 'z' extension.
    e_mesh, zs = _get_extension_side(2, grading_fun, 30,
                                     b_dims, b_shape, e_dims, e_shape, centre)
    mesh = mesh + e_mesh

    # Mirror by 'z'.
    e_mesh.coors[:, 2] = (2 * centre[2]) - e_mesh.coors[:, 2]
    e_mesh.cmesh.cell_groups.fill(31)
    mesh = mesh + e_mesh

    if name is not None:
        mesh.name = name

    # Verify merging by checking the number of nodes.
    n_nod = (nm.prod(nm.maximum(b_shape - 2, 0)) + 2 * nm.prod(xs)
             + 2 * (max(ys[0] - 2, 0) * ys[1] * ys[2])
             + 2 * (max(zs[0] - 2, 0) * max(zs[1] - 2, 0) * zs[2]))
    if n_nod != mesh.n_nod:
        raise ValueError('Merge of meshes failed! (%d == %d)'
                         % (n_nod, mesh.n_nod))

    return mesh

def tiled_mesh1d(conn, coors, ngrps, idim, n_rep, bb, eps=1e-6, ndmap=False):
    from sfepy.discrete.fem.periodic import match_grid_plane

    s1 = nm.nonzero(coors[:,idim] < (bb[0] + eps))[0]
    s2 = nm.nonzero(coors[:,idim] > (bb[1] - eps))[0]

    if s1.shape != s2.shape:
        raise ValueError('incompatible shapes: %s == %s'\
              % (s1.shape, s2.shape))

    (nnod0, dim) = coors.shape
    nnod = nnod0 * n_rep - s1.shape[0] * (n_rep - 1)
    (nel0, nnel) = conn.shape
    nel = nel0 * n_rep

    dd = nm.zeros((dim,), dtype=nm.float64)
    dd[idim] = bb[1] - bb[0]

    m1, m2 = match_grid_plane(coors[s1], coors[s2], idim)

    oconn = nm.zeros((nel, nnel), dtype=nm.int32)
    ocoors = nm.zeros((nnod, dim), dtype=nm.float64)
    ongrps = nm.zeros((nnod,), dtype=nm.int32)

    if type(ndmap) is bool:
        ret_ndmap = ndmap

    else:
        ret_ndmap= True
        ndmap_out = nm.zeros((nnod,), dtype=nm.int32)

    el_off = 0
    nd_off = 0

    for ii in range(n_rep):
        if ii == 0:
            oconn[0:nel0,:] = conn
            ocoors[0:nnod0,:] = coors
            ongrps[0:nnod0] = ngrps.squeeze()
            nd_off += nnod0

            mapto = s2[m2]
            mask = nm.ones((nnod0,), dtype=nm.int32)
            mask[s1] = 0
            remap0 = nm.cumsum(mask) - 1
            nnod0r = nnod0 - s1.shape[0]
            cidx = nm.where(mask)
            if ret_ndmap:
                ndmap_out[0:nnod0] = nm.arange(nnod0)

        else:
            remap = remap0 + nd_off
            remap[s1[m1]] = mapto
            mapto = remap[s2[m2]]

            ocoors[nd_off:(nd_off + nnod0r),:] =\
              (coors[cidx,:] + ii * dd)
            ongrps[nd_off:(nd_off + nnod0r)] = ngrps[cidx].squeeze()
            oconn[el_off:(el_off + nel0),:] = remap[conn]
            if ret_ndmap:
                ndmap_out[nd_off:(nd_off + nnod0r)] = cidx[0]

            nd_off += nnod0r

        el_off += nel0

    if ret_ndmap:
        if ndmap is not None:
            max_nd_ref = nm.max(ndmap)
            idxs = nm.where(ndmap_out > max_nd_ref)
            ndmap_out[idxs] = ndmap[ndmap_out[idxs]]

        return oconn, ocoors, ongrps, ndmap_out

    else:
        return oconn, ocoors, ongrps

def gen_tiled_mesh(mesh, grid=None, scale=1.0, eps=1e-6, ret_ndmap=False):
    """
    Generate a new mesh by repeating a given periodic element
    along each axis.

    Parameters
    ----------
    mesh : Mesh instance
        The input periodic FE mesh.
    grid : array
        Number of repetition along each axis.
    scale : float, optional
        Scaling factor.
    eps : float, optional
        Tolerance for boundary detection.
    ret_ndmap : bool, optional
        If True, return global node map.

    Returns
    -------
    mesh_out : Mesh instance
        FE mesh.
    ndmap : array
        Maps: actual node id --> node id in the reference cell.
    """
    bbox = mesh.get_bounding_box()

    if grid is None:
        iscale = max(int(1.0 / scale), 1)
        grid = [iscale] * mesh.dim

    conn = mesh.get_conn(mesh.descs[0])
    mat_ids = mesh.cmesh.cell_groups

    coors = mesh.coors
    ngrps = mesh.cmesh.vertex_groups
    nrep = nm.prod(grid)
    ndmap = None

    output('repeating %s ...' % grid)
    nblk = 1
    for ii, gr in enumerate(grid):
        if ret_ndmap:
            (conn, coors,
             ngrps, ndmap0) = tiled_mesh1d(conn, coors, ngrps,
                                           ii, gr, bbox.transpose()[ii],
                                           eps=eps, ndmap=ndmap)
            ndmap = ndmap0

        else:
            conn, coors, ngrps = tiled_mesh1d(conn, coors, ngrps,
                                              ii, gr, bbox.transpose()[ii],
                                              eps=eps)
        nblk *= gr

    output('...done')

    mat_ids = nm.tile(mat_ids, (nrep,))
    mesh_out = Mesh.from_data('tiled mesh', coors * scale, ngrps,
                              [conn], [mat_ids], [mesh.descs[0]])

    if ret_ndmap:
        return mesh_out, ndmap
    else:
        return mesh_out

def gen_misc_mesh(mesh_dir, force_create, kind, args, suffix='.mesh',
                  verbose=False):
    """
    Create sphere or cube mesh according to `kind` in the given
    directory if it does not exist and return path to it.
    """
    import os
    from sfepy import data_dir

    defdir = os.path.join(data_dir, 'meshes')
    if mesh_dir is None:
        mesh_dir = defdir

    def retype(args, types, defaults):
        args=list(args)
        args.extend(defaults[len(args):len(defaults)])
        return tuple([type(value) for type, value in zip(types, args) ])

    if kind == 'sphere':
        default = [5, 41, args[0]]
        args = retype(args, [float, int, float], default)
        mesh_pattern = os.path.join(mesh_dir, 'sphere-%.2f-%.2f-%i')

    else:
        assert_(kind == 'cube')

        args = retype(args,
                      (int, float, int, float, int, float),
                      (args[0], args[1], args[0], args[1], args[0], args[1]))
        mesh_pattern = os.path.join(mesh_dir, 'cube-%i_%.2f-%i_%.2f-%i_%.2f')

    if verbose:
        output(args)

    filename = mesh_pattern % args
    if not force_create:
        if os.path.exists(filename): return filename
        if os.path.exists(filename + '.mesh') : return filename + '.mesh'
        if os.path.exists(filename + '.vtk'): return filename + '.vtk'

    if kind == 'cube':
        filename = filename + suffix
        ensure_path(filename)

        output('creating new cube mesh')
        output('(%i nodes in %.2f) x (%i nodes in %.2f) x (%i nodes in %.2f)'
               % args)
        output('to file %s...' % filename)

        mesh = gen_block_mesh(args[1::2], args[0::2],
                              (0.0, 0.0, 0.0), name=filename)
        mesh.write(filename, io='auto')
        output('...done')

    else:
        import subprocess, shutil, tempfile
        filename = filename + '.mesh'
        ensure_path(filename)

        output('creating new sphere mesh (%i nodes, r=%.2f) and gradation %d'
               % args)
        output('to file %s...' % filename)

        f = open(os.path.join(defdir, 'quantum', 'sphere.geo'))
        tmp_dir = tempfile.mkdtemp()
        tmpfile = os.path.join(tmp_dir, 'sphere.geo.temp')
        ff = open(tmpfile, "w")
        ff.write("""
R = %i.0;
n = %i.0;
dens = %f;
""" % args)
        ff.write(f.read())
        f.close()
        ff.close()
        subprocess.call(['gmsh', '-3', tmpfile, '-format', 'mesh',
                         '-o', filename])
        shutil.rmtree(tmp_dir)
        output('...done')

    return filename

def gen_mesh_from_string(mesh_name, mesh_dir):
    import re
    result = re.match('^\\s*([a-zA-Z]+)[:\\(]([^\\):]*)[:\\)](\\*)?\\s*$',
                      mesh_name)

    if result is None:
        return mesh_name

    else:
        args = re.split(',', result.group(2))
        kind = result.group(1)
        return gen_misc_mesh(mesh_dir, result.group(3)=='*', kind, args)

def gen_mesh_from_geom(geo, a=None, verbose=False, refine=False):
    """
    Runs mesh generator - tetgen for 3D or triangle for 2D meshes.

    Parameters
    ----------
    geo : geometry
        geometry description
    a : int, optional
        a maximum area/volume constraint
    verbose : bool, optional
        detailed information
    refine : bool, optional
        refines mesh

    Returns
    -------
    mesh : Mesh instance
        triangular or tetrahedral mesh
    """

    import os.path as op
    import pexpect
    import tempfile
    import shutil

    tmp_dir = tempfile.mkdtemp()
    polyfilename = op.join(tmp_dir, 'meshgen.poly')

    # write geometry to poly file
    geo.to_poly_file(polyfilename)
    meshgen_call = {2: ('triangle', ''), 3: ('tetgen', 'BFENk')}

    params = "-ACp"
    params += "q" if refine else ''
    params += "V" if verbose else "Q"
    params += meshgen_call[geo.dim][1]
    if a is not None:
        params += "a%f" % (a)
    params += " %s" % (polyfilename)

    cmd = "%s %s" % (meshgen_call[geo.dim][0], params)
    if verbose: print("Generating mesh using", cmd)

    p=pexpect.run(cmd, timeout=None)
    bname, ext = op.splitext(polyfilename)
    if geo.dim == 2:
        mesh = Mesh.from_file(bname + '.1.node')
    if geo.dim == 3:
        mesh = Mesh.from_file(bname + '.1.vtk')

    shutil.rmtree(tmp_dir)

    return mesh

def gen_mesh_from_voxels(voxels, dims, etype='q'):
    """
    Generate FE mesh from voxels (volumetric data).

    Parameters
    ----------
    voxels : array
        Voxel matrix, 1=material.
    dims : array
        Size of one voxel.
    etype : integer, optional
        'q' - quadrilateral or hexahedral elements
        't' - triangular or tetrahedral elements
    Returns
    -------
    mesh : Mesh instance
        Finite element mesh.
    """

    dims = nm.array(dims).squeeze()
    dim = len(dims)
    nddims = nm.array(voxels.shape) + 2

    nodemtx = nm.zeros(nddims, dtype=nm.int32)

    if dim == 2:
        #iy, ix = nm.where(voxels.transpose())
        iy, ix = nm.where(voxels)
        nel = ix.shape[0]

        if etype == 'q':
            nodemtx[ix,iy] += 1
            nodemtx[ix + 1,iy] += 1
            nodemtx[ix + 1,iy + 1] += 1
            nodemtx[ix,iy + 1] += 1

        elif etype == 't':
            nodemtx[ix,iy] += 2
            nodemtx[ix + 1,iy] += 1
            nodemtx[ix + 1,iy + 1] += 2
            nodemtx[ix,iy + 1] += 1
            nel *= 2

    elif dim == 3:
        #iy, ix, iz = nm.where(voxels.transpose(1, 0, 2))
        iy, ix, iz = nm.where(voxels)
        nel = ix.shape[0]

        if etype == 'q':
            nodemtx[ix,iy,iz] += 1
            nodemtx[ix + 1,iy,iz] += 1
            nodemtx[ix + 1,iy + 1,iz] += 1
            nodemtx[ix,iy + 1,iz] += 1
            nodemtx[ix,iy,iz + 1] += 1
            nodemtx[ix + 1,iy,iz + 1] += 1
            nodemtx[ix + 1,iy + 1,iz + 1] += 1
            nodemtx[ix,iy + 1,iz + 1] += 1

        elif etype == 't':
            nodemtx[ix,iy,iz] += 6
            nodemtx[ix + 1,iy,iz] += 2
            nodemtx[ix + 1,iy + 1,iz] += 2
            nodemtx[ix,iy + 1,iz] += 2
            nodemtx[ix,iy,iz + 1] += 2
            nodemtx[ix + 1,iy,iz + 1] += 2
            nodemtx[ix + 1,iy + 1,iz + 1] += 6
            nodemtx[ix,iy + 1,iz + 1] += 2
            nel *= 6

    else:
        msg = 'incorrect voxel dimension! (%d)' % dim
        raise ValueError(msg)

    ndidx = nm.where(nodemtx)
    coors = nm.array(ndidx).transpose() * dims
    nnod = coors.shape[0]

    nodeid = -nm.ones(nddims, dtype=nm.int32)
    nodeid[ndidx] = nm.arange(nnod)

    # generate elements
    if dim == 2:
        elems = nm.array([nodeid[ix,iy],
                          nodeid[ix + 1,iy],
                          nodeid[ix + 1,iy + 1],
                          nodeid[ix,iy + 1]]).transpose()

    elif dim == 3:
        elems = nm.array([nodeid[ix,iy,iz],
                          nodeid[ix + 1,iy,iz],
                          nodeid[ix + 1,iy + 1,iz],
                          nodeid[ix,iy + 1,iz],
                          nodeid[ix,iy,iz + 1],
                          nodeid[ix + 1,iy,iz + 1],
                          nodeid[ix + 1,iy + 1,iz + 1],
                          nodeid[ix,iy + 1,iz + 1]]).transpose()

    if etype == 't':
        elems = elems_q2t(elems)

    eid = etype + str(dim)
    eltab = {'q2': 4, 'q3': 8, 't2': 3, 't3': 4}

    mesh = Mesh.from_data('voxel_data',
                          coors, nm.ones((nnod,), dtype=nm.int32),
                          [nm.ascontiguousarray(elems)],
                          [nm.ones((nel,), dtype=nm.int32)],
                          ['%d_%d' % (dim, eltab[eid])])

    return mesh


def main():
    mesh = gen_block_mesh(nm.array((1.0, 2.0, 3.0)),
                          nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                          name='')
    mesh.write('0.mesh', io = 'auto')

    mesh = gen_cylinder_mesh(nm.array((1.0, 1.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=False, open_angle = 0.0,
                             name='')
    mesh.write('1.mesh', io = 'auto')
    mesh = gen_cylinder_mesh(nm.array((1.0, 1.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.0,
                             name='')
    mesh.write('2.mesh', io = 'auto')
    mesh = gen_cylinder_mesh(nm.array((1.0, 1.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5,
                             name='')
    mesh.write('3.mesh', io = 'auto')

    mesh = gen_cylinder_mesh(nm.array((0.0, 0.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=False, open_angle = 0.0,
                             name='')
    mesh.write('4.mesh', io = 'auto')

    mesh = gen_cylinder_mesh(nm.array((0.0, 0.0, 1.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5,
                             name='')
    mesh.write('5.mesh', io = 'auto')

    mesh = gen_cylinder_mesh(nm.array((0.0, 0.0, 1.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5, non_uniform=True,
                             name='')
    mesh.write('6.mesh', io = 'auto')

    mesh = gen_cylinder_mesh(nm.array((0.5, 0.5, 1.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5, non_uniform=True,
                             name='')
    mesh.write('7.mesh', io = 'auto')

if __name__ == '__main__':
    main()
