import numpy as nm

from sfepy.base.base import output, assert_
from sfepy.base.progressbar import MyBar
from sfepy.base.ioutils import ensure_path
from sfepy.linalg import cycle
from sfepy.fem.mesh import Mesh
from sfepy.mesh.mesh_tools import elems_q2t

def gen_block_mesh(dims, shape, centre, mat_id=0, name='block', verbose=True):
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

    x0 = centre - 0.5 * dims
    dd = dims / (shape - 1)

    grid = nm.zeros(shape, dtype = nm.int32)
    n_nod = nm.prod(shape)
    coors = nm.zeros((n_nod, dim), dtype = nm.float64)

    bar = MyBar("       nodes:", verbose=verbose)
    bar.init(n_nod)
    for ii, ic in enumerate(cycle(shape)):
        grid[tuple(ic)] = ii
        coors[ii] = x0 + ic * dd
        if not (ii % 100):
            bar.update(ii)
    bar.update(ii + 1)

    n_el = nm.prod(shape - 1)
    mat_ids = nm.empty((n_el,), dtype = nm.int32)
    mat_ids.fill(mat_id)

    if (dim == 2):
        conn = nm.zeros((n_el, 4), dtype = nm.int32)
        bar = MyBar("       elements:", verbose=verbose)
        bar.init(n_el)
        for ii, (ix, iy) in enumerate(cycle(shape - 1)):
            conn[ii,:] = [grid[ix  ,iy], grid[ix+1,iy  ],
                          grid[ix+1,iy+1], grid[ix  ,iy+1]]
            if not (ii % 100):
                bar.update(ii)
        bar.update(ii + 1)
        desc = '2_4'

    else:
        conn = nm.zeros((n_el, 8), dtype = nm.int32)
        bar = MyBar("       elements:", verbose=verbose)
        bar.init(n_el)
        for ii, (ix, iy, iz) in enumerate(cycle(shape - 1)):
            conn[ii,:] = [grid[ix  ,iy  ,iz  ], grid[ix+1,iy  ,iz  ],
                          grid[ix+1,iy+1,iz  ], grid[ix  ,iy+1,iz  ],
                          grid[ix  ,iy  ,iz+1], grid[ix+1,iy  ,iz+1],
                          grid[ix+1,iy+1,iz+1], grid[ix  ,iy+1,iz+1]]
            if not (ii % 100):
                bar.update(ii)
        bar.update(ii + 1)
        desc = '3_8'

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
    bar = MyBar("       nodes:", verbose=verbose)
    bar.init(n_nod)
    ii = 0
    for ix in range(nr):
        a, b = ras[ix], rbs[ix]
        for iy, fi in enumerate(angles[:nnfi]):
            for iz, x in enumerate(xs):
                grid[ix,iy,iz] = ii
                coors[ii] = origin + [x, a * nm.cos(fi), b * nm.sin(fi)]
                if not (ii % 100):
                    bar.update(ii)
                ii += 1

                if not is_hollow and (ix == 0):
                    if iy > 0:
                        grid[ix,iy,iz] = grid[ix,0,iz]
                        ii -= 1
    assert_(ii == n_nod)

    n_el = (nr - 1) * nnfi * (nl - 1)
    conn = nm.zeros((n_el, 8), dtype=nm.int32)

    bar = MyBar("       elements:", verbose=verbose)
    bar.init(n_el)
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

        if not (ii % 100):
            bar.update(ii)

    mat_id = nm.zeros((n_el,), dtype = nm.int32)
    desc = '3_8'

    assert_(n_nod == (conn.max() + 1))

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
    e_mesh.mat_ids[0].fill(11)
    mesh = mesh + e_mesh

    # 'y' extension.
    e_mesh, ys = _get_extension_side(1, grading_fun, 20,
                                     b_dims, b_shape, e_dims, e_shape, centre)
    mesh = mesh + e_mesh

    # Mirror by 'y'.
    e_mesh.coors[:, 1] = (2 * centre[1]) - e_mesh.coors[:, 1]
    e_mesh.mat_ids[0].fill(21)
    mesh = mesh + e_mesh

    # 'z' extension.
    e_mesh, zs = _get_extension_side(2, grading_fun, 30,
                                     b_dims, b_shape, e_dims, e_shape, centre)
    mesh = mesh + e_mesh

    # Mirror by 'z'.
    e_mesh.coors[:, 2] = (2 * centre[2]) - e_mesh.coors[:, 2]
    e_mesh.mat_ids[0].fill(31)
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

def tiled_mesh1d(conns, coors, ngrps, idim, n_rep, bb,
                 eps=1e-6, mybar=None, ndmap=False):
    from sfepy.fem.periodic import match_grid_plane

    s1 = nm.nonzero(coors[:,idim] < (bb[0] + eps))[0]
    s2 = nm.nonzero(coors[:,idim] > (bb[1] - eps))[0]

    if s1.shape != s2.shape:
        raise ValueError, 'incompatible shapes: %s == %s'\
              % (s1.shape, s2.shape)

    (nnod0, dim) = coors.shape
    nnod = nnod0 * n_rep - s1.shape[0] * (n_rep - 1)
    (nel0, nnel) = conns.shape
    nel = nel0 * n_rep

    dd = nm.zeros((dim,), dtype=nm.float64)
    dd[idim] = bb[1] - bb[0]

    m1, m2 = match_grid_plane(coors[s1], coors[s2], idim)

    oconns = nm.zeros((nel, nnel), dtype=nm.int32)
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
            oconns[0:nel0,:] = conns
            ocoors[0:nnod0,:] = coors
            ongrps[0:nnod0] = ngrps
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
            ongrps[nd_off:(nd_off + nnod0r)] = ngrps[cidx]
            oconns[el_off:(el_off + nel0),:] = remap[conns]
            if ret_ndmap:
                ndmap_out[nd_off:(nd_off + nnod0r)] = cidx[0]

            nd_off += nnod0r

        el_off += nel0


        if mybar is not None:
            mybar[0].update(mybar[1])

    if ret_ndmap:
        if ndmap is not None:
            max_nd_ref = nm.max(ndmap)
            idxs = nm.where(ndmap_out > max_nd_ref)
            ndmap_out[idxs] = ndmap[ndmap_out[idxs]]

        return oconns, ocoors, ongrps, ndmap_out

    else:
        return oconns, ocoors, ongrps

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

    conns = mesh.conns[0]
    for ii in mesh.conns[1:]:
        conns = nm.vstack((conns, ii))
    mat_ids = mesh.mat_ids[0]
    for ii in mesh.mat_ids[1:]:
        mat_ids = nm.hstack((mat_ids, ii))

    coors = mesh.coors
    ngrps = mesh.ngroups
    nrep = nm.prod(grid)
    ndmap = None

    bar = MyBar("       repeating:")
    bar.init(nrep)
    nblk = 1
    for ii, gr in enumerate(grid):
        if ret_ndmap:
            (conns, coors,
             ngrps, ndmap0) = tiled_mesh1d(conns, coors, ngrps,
                                           ii, gr, bbox.transpose()[ii],
                                           eps=eps, mybar=(bar, nblk),
                                           ndmap=ndmap)
            ndmap = ndmap0

        else:
            conns, coors, ngrps = tiled_mesh1d(conns, coors, ngrps,
                                               ii, gr, bbox.transpose()[ii],
                                               eps=eps, mybar=(bar, nblk))
        nblk *= gr

    bar.update(nblk)

    mat_ids = nm.tile(mat_ids, (nrep,))
    mesh_out = Mesh.from_data('tiled mesh', coors * scale, ngrps,
                              [conns], [mat_ids], [mesh.descs[0]])

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
        import subprocess
        filename = filename + '.mesh'
        ensure_path(filename)

        output('creating new sphere mesh (%i nodes, r=%.2f) and gradation %d'
               % args)
        output('to file %s...' % filename)

        f = open(os.path.join(defdir, 'quantum', 'sphere.geo'))
        tmpfile = os.path.join(data_dir, 'tmp', 'sphere.geo.temp')
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

def gen_mesh_from_goem(geo, a=None, quadratic=False, verbose=True,
                       refine=False, polyfilename='./meshgen.poly',
                       out='mesh', **kwargs):
    """
    Runs mesh generator - tetgen for 3D or triangle for 2D meshes.

    Parameters
    ----------
    geo : geometry
        geometry description
    a : int, optional
        a maximum area/volume constraint
    quadratic : bool, optional
        set True for quadratic elements
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

    # write geometry to poly file
    geo.to_poly_file(polyfilename)

    if not refine:
        params = "-Apq"
    else:
        params = "-Arq"
    if verbose:
        params = params + " -Q"
    if a != None and not refine:
        params = params + " -a%f" % (a)
    if refine:
        params = params + " -a"
    if quadratic:
        params = params + " -o2"
    params = params + " %s" % (polyfilename)

    meshgen_call = {2: 'triangle', 3: 'tetgen'}
    cmd = "%s %s" % (meshgen_call[geo.dim], params)
    if verbose: print "Generating mesh using", cmd
    if geo.dim == 2:
        p=pexpect.run(cmd, timeout=None)
        bname, ext = op.splitext(polyfilename)
        mesh = Mesh.from_file(bname + '.1.node')
        mesh.write(bname + '.' + out)
    if geo.dim == 3:
        p=pexpect.spawn(cmd, timeout=None)
        if not refine:
            p.expect("Opening %s." % (polyfilename))
        else:
            p.expect("Opening %s.node.\r\n" % (polyfilename))
            p.expect("Opening %s.ele.\r\n" % (polyfilename))
            p.expect("Opening %s.face.\r\n" % (polyfilename))
            p.expect("Opening %s.vol." % (polyfilename))
        assert p.before == ""
        p.expect(pexpect.EOF)
        if p.before != "\r\n":
            print p.before
            raise "Error when running mesh generator (see above for output): %s" % cmd

# http://www.cs.cmu.edu/~quake/triangle.html
#
# triangle [-prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh] input_file
#     -p  Triangulates a Planar Straight Line Graph (.poly file).
#     -r  Refines a previously generated mesh.
#     -q  Quality mesh generation.  A minimum angle may be specified.
#     -a  Applies a maximum triangle area constraint.
#     -u  Applies a user-defined triangle constraint.
#     -A  Applies attributes to identify triangles in certain regions.
#     -c  Encloses the convex hull with segments.
#     -D  Conforming Delaunay:  all triangles are truly Delaunay.
#     -j  Jettison unused vertices from output .node file.
#     -e  Generates an edge list.
#     -v  Generates a Voronoi diagram.
#     -n  Generates a list of triangle neighbors.
#     -g  Generates an .off file for Geomview.
#     -B  Suppresses output of boundary information.
#     -P  Suppresses output of .poly file.
#     -N  Suppresses output of .node file.
#     -E  Suppresses output of .ele file.
#     -I  Suppresses mesh iteration numbers.
#     -O  Ignores holes in .poly file.
#     -X  Suppresses use of exact arithmetic.
#     -z  Numbers all items starting from zero (rather than one).
#     -o2 Generates second-order subparametric elements.
#     -Y  Suppresses boundary segment splitting.
#     -S  Specifies maximum number of added Steiner points.
#     -i  Uses incremental method, rather than divide-and-conquer.
#     -F  Uses Fortune's sweepline algorithm, rather than d-and-c.
#     -l  Uses vertical cuts only, rather than alternating cuts.
#     -s  Force segments into mesh by splitting (instead of using CDT).
#     -C  Check consistency of final mesh.
#     -Q  Quiet:  No terminal output except errors.
#     -V  Verbose:  Detailed information on what I'm doing.
#     -h  Help:  Detailed instructions for Triangle.

# http://tetgen.berlios.de/
#
# tetgen [-prq_a_AiMYS_T_dzo_fenvgGOJBNEFICQVh] input_file
#     -p  Tetrahedralizes a piecewise linear complex (PLC).
#     -r  Reconstructs a previously generated mesh.
#     -q  Refines mesh (to improve mesh quality).
#     -a  Applies a maximum tetrahedron volume constraint.
#     -A  Assigns attributes to tetrahedra in different regions.
#     -i  Inserts a list of additional points into mesh.
#     -M  No merge of coplanar facets.
#     -Y  No splitting of input boundaries (facets and segments).
#     -S  Specifies maximum number of added points.
#     -T  Sets a tolerance for coplanar test (default 1e-8).
#     -d  Detects self-intersections of facets of the PLC.
#     -z  Numbers all output items starting from zero.
#     -o2 Generates second-order subparametric elements.
#     -f  Outputs all faces to .face file.
#     -e  Outputs all edges to .edge file.
#     -n  Outputs tetrahedra neighbors to .neigh file.
#     -v  Outputs Voronoi diagram to files.
#     -g  Outputs mesh to .mesh file for viewing by Medit.
#     -G  Outputs mesh to .msh file for viewing by Gid.
#     -O  Outputs mesh to .off file for viewing by Geomview.
#     -K  Outputs mesh to .vtk file for viewing by Paraview.
#     -J  No jettison of unused vertices from output .node file.
#     -B  Suppresses output of boundary information.
#     -N  Suppresses output of .node file.
#     -E  Suppresses output of .ele file.
#     -F  Suppresses output of .face file.
#     -I  Suppresses mesh iteration numbers.
#     -C  Checks the consistency of the final mesh.
#     -Q  Quiet:  No terminal output except errors.
#     -V  Verbose:  Detailed information, more terminal output.
#     -h  Help:  A brief instruction for using TetGen.

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

    dims = dims.squeeze()
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
                          {0: nm.ascontiguousarray(elems)},
                          {0: nm.ones((nel,), dtype=nm.int32)},
                          {0: '%d_%d' % (dim, eltab[eid])})

    return mesh

def gen_mesh_from_poly(filename, verbose=True):
    """
    Import mesh generated by tetgen or triangle.

    Parameters
    ----------
    filename : string
        file name

    Returns
    -------
    mesh : Mesh instance
        triangular or tetrahedral mesh
    """

    def getnodes(fnods,up):
        f=file(fnods)
        l=[int(x) for x in f.readline().split()]
        npoints,dim,nattrib,nbound=l
        if verbose: up.init(npoints)
        nodes=[]
        for line in f:
            if line[0]=="#": continue
            l=[float(x) for x in line.split()]
            l = l[:(dim + 1)]
            l[0]=int(l[0])
            nodes.append(tuple(l))
            assert l[0]==len(nodes)
        assert npoints==len(nodes)
        return nodes

    def getele(fele,up):
        f=file(fele)
        l=[int(x) for x in f.readline().split()]
        nele,nnod,nattrib=l
        #we have either linear or quadratic tetrahedra:
        if nnod in [4,10]:
            elem = 'tetra'
            linear = (nnod == 4)
        if nnod in [3, 7]:
            elem = 'tri'
            linear = (nnod == 3)

        # if nattrib!=1:
        #     raise "tetgen didn't assign an entity number to each element (option -A)"
        els=[]
        regions={}
        for line in f:
            if line[0]=="#": continue
            l=[int(x) for x in line.split()]
            if elem == 'tri':
                if linear:
                    assert (len(l) - 1 - nattrib) == 3
                    els.append((l[0],l[1],l[2],l[3]))
                    regionnum=l[5]
                else:
                    assert len(l)-2 == 10
                    els.append((l[0],54,l[1],l[2],l[3],l[4],
                                l[5],l[6],l[7],l[8],l[9],l[10]))
                    regionnum=l[11]
            if elem == 'tetra':
                if linear:
                    assert len(l)-2 == 4
                    els.append((l[0],54,l[1],l[2],l[3],l[4]))
                    regionnum=l[5]
                else:
                    assert len(l)-2 == 10
                    els.append((l[0],54,l[1],l[2],l[3],l[4],
                                l[5],l[6],l[7],l[8],l[9],l[10]))
                    regionnum=l[11]
            if regionnum==0:
                print "see %s, element # %d"%(fele,l[0])
                raise "there are elements not belonging to any physical entity"
            if regions.has_key(regionnum):
                regions[regionnum].append(l[0])
            else:
                regions[regionnum]=[l[0]]
            assert l[0]==len(els)
            if verbose: up.update(l[0])
        return els,regions,linear

    def getBCfaces(ffaces,up):
        f=file(ffaces)
        l=[int(x) for x in f.readline().split()]
        nfaces,nattrib=l
        if nattrib!=1:
            raise "tetgen didn't assign an entity number to each face \
(option -A)"
        if verbose: up.init(nfaces)
        faces={}
        for line in f:
            if line[0]=="#": continue
            l=[int(x) for x in line.split()]
            assert len(l)==5
            regionnum=l[4]
            if regionnum==0: continue
            if faces.has_key(regionnum):
                faces[regionnum].append((l[1],l[2],l[3]))
            else:
                faces[regionnum]=[(l[1],l[2],l[3])]
            if verbose: up.update(l[0])
        return faces

    def calculatexyz(nodes, els):
        """Calculate the missing xyz values in place"""
        def avg(i,j,n4,nodes):
            a=nodes[n4[i-1]-1]
            b=nodes[n4[j-1]-1]
            return (a[1]+b[1])/2, (a[2]+b[2])/2, (a[3]+b[3])/2
        def getxyz(i,n4,nodes):
            if i+5==5: return avg(1,2,n4,nodes)
            if i+5==6: return avg(2,3,n4,nodes)
            if i+5==7: return avg(1,3,n4,nodes)
            if i+5==8: return avg(1,4,n4,nodes)
            if i+5==9: return avg(2,4,n4,nodes)
            if i+5==10: return avg(3,4,n4,nodes)
            raise "wrong topology"
        for e in els:
            n4=e[2:2+4]
            n6=e[2+4:2+4+10]
            for i,n in enumerate(n6):
                x,y,z=getxyz(i,n4,nodes)
                nodes[n-1]=(n,x,y,z)

    if verbose: print "Reading geometry from poly file..."
    m=Mesh()
    m.nodes=getnodes(filename+".node")
    m.elements,m.regions, lin=getele(filename+".ele")
    if not lin:
        #tetgen doesn't compute xyz coordinates of the aditional 6 nodes
        #(only of the 4 corner nodes) in tetrahedra.
        calculatexyz(m.nodes,m.elements)
    m.faces=getBCfaces(filename+".face")

    return m

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
