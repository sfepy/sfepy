import numpy as nm

from sfepy.base.base import output, assert_
from sfepy.base.progressbar import MyBar
from sfepy.linalg import cycle
from sfepy.fem.mesh import Mesh
from sfepy.fem.mesh import find_map, merge_mesh, make_mesh, fix_double_nodes
from sfepy.fem.mesh import get_min_edge_size, get_min_vertex_distance

def gen_block_mesh(dims, shape, centre, name='block'):
    """Generate a 2D or 3D block mesh. The dimension is determined by the
    lenght of the shape argument. 

    Parameters
    ----------

    dims : array of 2 or 3 floats
        Dimensions of the block.
    shape : array of 2 or 3 ints
        Shape (counts of nodes in x, y, z) of the block mesh.
    centre : array of 2 or 3 floats
        Centre of the block.

    name : string
        Mesh name.

    Returns
    -------

    mesh : Mesh instance
    """
    dims = nm.asarray(dims)
    shape = nm.asarray(shape)
    centre = nm.asarray(centre)

    dim = shape.shape[0]

    centre = centre[:dim]
    dims = dims[:dim]

    x0 = centre - 0.5 * dims
    dd = dims / (shape - 1)

    grid = nm.zeros( shape, dtype = nm.int32 )
    n_nod = nm.prod( shape )
    coors = nm.zeros( (n_nod, dim), dtype = nm.float64 )
    
    bar = MyBar( "       nodes:" )
    bar.init( n_nod )
    for ii, ic in enumerate( cycle( shape ) ):
        grid[tuple(ic)] = ii
        coors[ii] = x0 + ic * dd
        if not (ii % 100):
            bar.update( ii )
    bar.update(ii + 1)

    n_el = nm.prod( shape - 1 )
    mat_id = nm.zeros( (n_el,), dtype = nm.int32 )

    if (dim == 2):
        conn = nm.zeros( (n_el, 4), dtype = nm.int32 )
        bar = MyBar( "       elements:" )
        bar.init( n_el )
        for ii, (ix, iy) in enumerate( cycle( shape - 1 ) ):
            conn[ii,:] = [grid[ix  ,iy], grid[ix+1,iy  ],
                          grid[ix+1,iy+1], grid[ix  ,iy+1]]
            if not (ii % 100):
                bar.update( ii )
        bar.update(ii + 1)
        desc = '2_4'

    else:
        conn = nm.zeros( (n_el, 8), dtype = nm.int32 )
        bar = MyBar( "       elements:" )
        bar.init( n_el )
        for ii, (ix, iy, iz) in enumerate( cycle( shape - 1 ) ):
            conn[ii,:] = [grid[ix  ,iy  ,iz  ], grid[ix+1,iy  ,iz  ],
                          grid[ix+1,iy+1,iz  ], grid[ix  ,iy+1,iz  ],
                          grid[ix  ,iy  ,iz+1], grid[ix+1,iy  ,iz+1],
                          grid[ix+1,iy+1,iz+1], grid[ix  ,iy+1,iz+1]]
            if not (ii % 100):
                bar.update( ii )
        bar.update(ii + 1)
        desc = '3_8'

    mesh = Mesh.from_data(name, coors, None, [conn], [mat_id], [desc])
    return mesh

def gen_cylinder_mesh(dims, shape, centre, axis='x', force_hollow=False,
                      is_open=False, open_angle=0.0, non_uniform=False,
                      name='cylinder'):
    """Generate a cylindrical mesh along an axis. Its cross-section can be
    ellipsoidal.

    Parameters
    ----------

    axis: one of 'x', 'y', 'z'
        The axis of the cylinder.
        
    dims : array of 5 floats
        Dimensions of the cylinder: inner surface semi-axes a1, b1, outer
        surface semi-axes a2, b2, length.
    shape : array of 3 ints
        Shape (counts of nodes in radial, circumferential and longitudinal
        directions) of the cylinder mesh.
    centre : array of 3 floats
        Centre of the cylinder.

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

    Returns
    -------

    mesh : Mesh instance
    """
    dims = nm.asarray(dims)
    shape = nm.asarray(shape)
    centre = nm.asarray(centre)

    a1, b1, a2, b2, length = dims
    nr, nfi, nl = shape
    origin = centre - nm.array([0.5 * length, 0.0, 0.0])

    da = (a2 - a1) / (nr - 1)
    db = (b2 - b1) / (nr - 1)
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

       
##     print dfi * 180.0 / nm.pi
##     print angles * 180.0 / nm.pi
##     print xs
##     print ras
##     print rbs

    # This is 3D only...
    bar = MyBar( "       nodes:" )
    bar.init( n_nod )
    ii = 0
    for ix in range(nr): 
        a, b = ras[ix], rbs[ix]
        for iy, fi in enumerate(angles[:nnfi]):
#            print iy, fi * 180.0 / nm.pi
            for iz, x in enumerate(xs):
##                 print ix, iy, iz, ii
                grid[ix,iy,iz] = ii
                coors[ii] = origin + [x, a * nm.cos(fi), b * nm.sin(fi)]
                if not (ii % 100):
                    bar.update( ii )
                ii += 1

                if not is_hollow and (ix == 0):
                    if iy > 0:
                        grid[ix,iy,iz] = grid[ix,0,iz]
                        ii -= 1
    print
    assert_(ii == n_nod)

    n_el = (nr - 1) * nnfi * (nl - 1)
    conn = nm.zeros((n_el, 8), dtype=nm.int32)

    bar = MyBar( "       elements:" )
    bar.init(n_el)
    ii = 0
    for (ix, iy, iz) in cycle([nr-1, nnfi, nl-1]):
#        print ii, ix, iy, iz
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
            bar.update( ii )
    print
    mat_id = nm.zeros( (n_el,), dtype = nm.int32 )
    desc = '3_8'

##     print n_nod, n_el, conn.max()
    assert_(n_nod == (conn.max() + 1))

    if axis == 'z':
        coors = coors[:,[1,2,0]]
    elif axis == 'y':
        coors = coors[:,[2,0,1]]

    mesh = Mesh.from_data(name, coors, None, [conn], [mat_id], [desc])
    return mesh

def compose_periodic_mesh(mesh_in, scale, repeat, eps, check_mvd=False):
    """
    Create a new mesh from a periodic mesh by scaling it and repeating
    along each axis,

    Parameters
    ----------
    mesh_in : Mesh instance
        The input periodic mesh.
    scale : float
        The scale parameter.
    repeat : array
        The number of repetitions of `mesh_in` along each axis.
    eps : float
        Tolerance for finding matching mesh vertices.
    check_mvd : bool
        If True, verify periodicity by checking approximate minimum
        vertex distance and minimum edge length.

    Returns
    -------
    mesh_out : Mesh instance
        The composed periodic mesh.
    """
    dim = mesh_in.dim
    bbox = mesh_in.get_bounding_box()
    mscale = bbox[1] - bbox[0]
    output('bbox:\n', bbox)

    centre0 = 0.5 * (bbox[1] + bbox[0])
    output('centre:\n', centre0)

    scale = float(scale)
    if repeat is None:
        repeat = [scale] * dim

    repeat = nm.asarray(repeat)

    # Normalize original coordinates.
    coor0 = (mesh_in.coors - centre0) / (mscale)

    aux = fix_double_nodes(coor0, mesh_in.ngroups, mesh_in.conns, eps)
    coor0, ngroups0, mesh_in.conns = aux
    
    if check_mvd:
        mes0 = get_min_edge_size(coor0, mesh_in.conns)
        mvd0 = get_min_vertex_distance(coor0, mes0)
        if mes0 > (mvd0 + eps):
            output('          original min. "edge" length: %.5e' % mes0)
            output('original approx. min. vertex distance: %.5e' % mvd0)
            output('-> still double nodes in input mesh!')
            output('try increasing eps')
            raise ValueError('cannot fix double nodes!')

    output('composing periodic mesh...')
    for indx in cycle(repeat):
        aindx = nm.array(indx, dtype=nm.float64)
        centre = 0.5 * (2.0 * aindx - repeat + 1.0)
        output(indx, centre)

        if aindx.sum() == 0:
            coor = coor0 + centre
            ngroups = ngroups0
            conns = mesh_in.conns
        else:
            coor1 = coor0 + centre
            ngroups1 = ngroups0
            conns1 = mesh_in.conns

            cmap = find_map(coor, coor1, eps=eps)
            if not cmap.size:
                raise ValueError('non-periodic mesh!')

            else:
                output(cmap.size / 2)

            coor, ngroups, conns = merge_mesh(coor, ngroups, conns,
                                              coor1, ngroups1, conns1,
                                              cmap, eps=eps)

    if check_mvd:
        mes = get_min_edge_size( coor, conns )
        mvd = get_min_vertex_distance( coor, mes0 )

        output('          original min. "edge" length: %.5e' % mes0)
        output('             final min. "edge" length: %.5e' % mes)
        output('original approx. min. vertex distance: %.5e' % mvd0)
        output('   final approx. min. vertex distance: %.5e' % mvd)
        if mvd < 0.99999 * mvd0:
            if mvd0 < (mes0 - eps):
                output('-> probably non-periodic input mesh!')
                output('   - adjacent sides were not connected!')
                output('   try increasing eps')
            else:
                output('-> input mesh might be periodic')
                output('   try increasing eps')
        else:
            output('-> input mesh looks periodic')
    else:
        output('non-periodic input mesh detection skipped!')

    # Renormalize.
    coor = (coor * mscale) / scale
    mesh_out = make_mesh(coor, ngroups, conns, mesh_in)
    output('...done')

    return mesh_out

def main():
    mesh = gen_block_mesh(nm.array((1.0, 2.0, 3.0)),
                          nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                          name='')
    mesh.write('0.mesh', io = 'auto' )

    mesh = gen_cylinder_mesh(nm.array((1.0, 1.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=False, open_angle = 0.0,
                             name='')
    mesh.write('1.mesh', io = 'auto' )
    mesh = gen_cylinder_mesh(nm.array((1.0, 1.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.0,
                             name='')
    mesh.write('2.mesh', io = 'auto' )
    mesh = gen_cylinder_mesh(nm.array((1.0, 1.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5,
                             name='')
    mesh.write('3.mesh', io = 'auto' )

    mesh = gen_cylinder_mesh(nm.array((0.0, 0.0, 2.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=False, open_angle = 0.0,
                             name='')
    mesh.write('4.mesh', io = 'auto' )

    mesh = gen_cylinder_mesh(nm.array((0.0, 0.0, 1.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5,
                             name='')
    mesh.write('5.mesh', io = 'auto' )

    mesh = gen_cylinder_mesh(nm.array((0.0, 0.0, 1.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5, non_uniform=True,
                             name='')
    mesh.write('6.mesh', io = 'auto' )

    mesh = gen_cylinder_mesh(nm.array((0.5, 0.5, 1.0, 2.0, 3)),
                             nm.array((10,10,10)), nm.array((1.0, 2.0, 3.0)),
                             is_open=True, open_angle = 0.5, non_uniform=True,
                             name='')
    mesh.write('7.mesh', io = 'auto' )

if __name__ == '__main__':
    main()
