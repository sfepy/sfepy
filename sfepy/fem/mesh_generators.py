from sfepy.base.base import *
from sfepy.base.progressbar import MyBar
from sfepy.base.la import cycle
from sfepy.fem import Mesh

def gen_block_mesh(dims, shape, centre, name='block'):
    """Generate a block mesh.

    Parameters
    ----------

    dims : array of 3 floats
        Dimensions of the block.
    shape : array of 3 ints
        Shape (counts of nodes in x, y, z) of the block mesh.
    centre : array of 3 floats
        Centre of the block.

    name : string
        Mesh name.

    Returns
    -------

    mesh : Mesh instance
    """
    dim = shape.shape[0]

    x0 = centre - 0.5 * dims
    dd = dims / (shape - 1)

    grid = nm.zeros( shape, dtype = nm.int32 )
    n_nod = nm.prod( shape )
    coors = nm.zeros( (n_nod, dim), dtype = nm.float64 )

    # This is 3D only...
    bar = MyBar( "       nodes:" )
    bar.init( n_nod )
    for ii, ic in enumerate( cycle( shape ) ):
        ix, iy, iz = ic
        grid[ix,iy,iz] = ii
        coors[ii] = x0 + ic * dd
        if not (ii % 100):
            bar.update( ii )
    print
    n_el = nm.prod( shape - 1 )
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
    print
    mat_id = nm.zeros( (n_el,), dtype = nm.int32 )
    desc = '3_8'

    mesh = Mesh.from_data(name, coors, None, [conn], [mat_id], [desc])
    return mesh

def gen_cylinder_mesh(dims, shape, centre, force_hollow=False,
                      is_open=False, open_angle=0.0, non_uniform=False,
                      name='cylinder'):
    """Generate a cylindrical mesh along the x axis. Its cross-section can be
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

    mesh = Mesh.from_data(name, coors, None, [conn], [mat_id], [desc])
    return mesh

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
