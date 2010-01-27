def gen_concentric( filename, a, el_size_out, el_size_in, r0, r1, n_circ ):
    fd = open( filename, 'w' )
    fd.write( 'a = %f;\n' % a )
    fd.write( 'out = %f;\n' % el_size_out )
    fd.write( 'in = %f;\n' % el_size_in )
    fd.write( 'Point(1) = {-a,-a,0, out};\n' )
    fd.write( 'Point(2) = {-a,a,0, out};\n' )
    fd.write( 'Point(3) = {a, a,0, out};\n' )
    fd.write( 'Point(4) = {a,-a,0, out};\n' )
    fd.write( 'Point(5) = {0,0,0, in};\n' )
    fd.write( 'Line(1) = {1,4};\n' )
    fd.write( 'Line(2) = {4,3};\n' )
    fd.write( 'Line(3) = {3,2};\n' )
    fd.write( 'Line(4) = {2,1};\n' )
    fd.write( 'Line Loop(1) = {2,3,4,1};\n' )

    r = r0
    if n_circ > 1:
        dr = (r1 - r0) / (n_circ - 1)
    else:
        dr = r1 - r0

    ii = 6
    il = 2
    in_groups = []
    rs = []
    for ir in xrange( n_circ ):
        fd.write( 'r = %f;\n' % r )
        fd.write( 'Point(%d) = {r,0,0, in};\n' % ii )
        fd.write( 'Point(%d) = {0,r,0, in};\n' % (ii+1) )
        fd.write( 'Point(%d) = {-r, 0,0, in};\n' % (ii+2) )
        fd.write( 'Point(%d) = {0,-r,0, in};\n' % (ii+3) )
        fd.write( 'Circle(%d) = {%d,5,%d};\n' % (ii,ii,ii+1) )
        fd.write( 'Circle(%d) = {%d,5,%d};\n' % (ii+1,ii+1,ii+2) )
        fd.write( 'Circle(%d) = {%d,5,%d};\n' % (ii+2,ii+2,ii+3) )
        fd.write( 'Circle(%d) = {%d,5,%d};\n' % (ii+3,ii+3,ii) )
        fd.write( 'Line Loop(%d) = {%d,%d,%d,%d};\n' % (il,
                                                        ii, ii+1, ii+2, ii+3) )
        in_groups.append( il )
        rs.append( r )
        ii += 4
        il += 1
        r += dr

    fd.write( 'Plane Surface(1) = {1,%d};\n' % (il - 1) )
    fd.write( 'Plane Surface(2) = {2};\n' )
    for ir in range( 1, n_circ ):
        fd.write( 'Plane Surface(%d) = {%d,%d};\n' % (ir+2,ir+2,ir+1) )

    fd.write( 'Physical Surface(1) = {1};\n' )
    aux = ['%d' % ii for ii in in_groups]
    fd.write( 'Physical Surface(2) = {%s};\n' % ','.join( aux ) )

    fd.close()

    return in_groups, rs
