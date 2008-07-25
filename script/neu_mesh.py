#!/usr/bin/python
# 13.01.2005, c
# 17.01.2005
# 03.02.2005
# 31.01.2006
# 29.05.2007
import sys
import fileinput

if (len( sys.argv ) == 3):
    filename_in = sys.argv[1];
    filename_out = sys.argv[2];
else:
    print 'Two args required!'
    raise ValueError

if (filename_out == '-'):
    file_out = sys.stdout
else:
    file_out = open( filename_out, "w" ); 

mode = 0
nod = []
els = {3 : {'tetra4' : [], 'pyram5' : [], 'wedge6' : [], 'brick8' : []},
       2 : {'tria3' : [], 'quad4' : []}}

group_ids = []
group_n_els = []
groups = []

input = fileinput.input( filename_in )
for line in input:
    row = line.split()
    if len( row ) == 0: continue

    if (row[0] == 'NUMNP'):
        row = input.readline().split()
        n_nod, n_el, dim = row[0], row[1], int( row[4] )
        print n_nod, n_el, dim
        el = els[dim]

    elif (row[0] == 'NODAL'):
        row = input.readline().split()
        while not( row[0] == 'ENDOFSECTION' ):
            nod.append( row[1:] )
            row = input.readline().split()

    elif (row[0] == 'ELEMENTS/CELLS'):
        if dim == 3:
            row = input.readline().split()
            while not( row[0] == 'ENDOFSECTION' ):
    #            print row
                if (int( row[2] ) == 4):
                    el['tetra4'].append( row )
                if (int( row[2] ) == 5):
                    el['pyram5'].append( row )
                if (int( row[2] ) == 6):
                    el['wedge6'].append( row )
                elif (int( row[2] ) == 8):
                    rr = row[1:]
                    if (len( rr ) < 10):
                        rr.extend( input.readline().split() )
                    el['brick8'].append( rr )
                row = input.readline().split()
        else:
            row = input.readline().split()
            while not( row[0] == 'ENDOFSECTION' ):
    #            print row
                if (int( row[2] ) == 3):
                    el['tria3'].append( row )
                if (int( row[2] ) == 4):
                    el['quad4'].append( row )
                row = input.readline().split()

    elif (row[0] == 'GROUP:'):
        group_ids.append( row[1] )
        g_n_el = int( row[3] )
        group_n_els.append( g_n_el )
        name = input.readline().strip()
        print name, group_ids[-1]
        
        els = []
        row = input.readline().split()
        row = input.readline().split()
        while not( row[0] == 'ENDOFSECTION' ):
            els.extend( row )
            row = input.readline().split()
        if g_n_el != len( els ):
            print 'wrong number of group elements! (%d == %d)'\
                  % (n_el, len( els ))
            raise ValueError
        groups.append( els )
        
if int( n_el ) != sum( group_n_els ):
    print 'wrong total number of group elements! (%d == %d)'\
          % (int( n_el ), len( group_n_els ))

mat_ids = [None] * int( n_el )
for ii, els in enumerate( groups ):
    for iel in els:
        mat_ids[int( iel ) - 1] = group_ids[ii]

out = 'els: '
total = 0
for key, lst in el.iteritems():
    num = len( lst )
    out += ' %s: %d,' % (key, num)
    total += num
out += ' total: %d' % total
print out

file_out.write( """MeshVersionFormatted 1
Dimension %d
""" % dim )
    
file_out.write( "Vertices\n%d\n" % len( nod ) )
for nn in nod:
    file_out.write( " ".join( nn ) + " 0\n" )

if dim == 3:
    if (len( el['tetra4'] ) > 0):
        file_out.write( "Tetrahedra\n%d\n" % len( el['tetra4'] ) )
        for ee in el['tetra4']:
            ii = int( ee[0] ) - 1
            file_out.write( " ".join( ee[3:] ) + " " + mat_ids[ii] + "\n" )

    n_hex = len( el['pyram5'] ) + len( el['wedge6'] ) + len( el['brick8'] )
    if (n_hex > 0):
        file_out.write( "Hexahedra\n%d\n" % n_hex )
        for ee in el['brick8']:
            ii = int( ee[0] ) - 1
            file_out.write( " ".join( (ee[3], ee[4], ee[6], ee[5],
                                      ee[7], ee[8], ee[10], ee[9]) )
                           + " " + mat_ids[ii] + "\n" )
        for ee in el['wedge6']:
            ii = int( ee[0] ) - 1
            file_out.write( " ".join( (ee[3], ee[5], ee[8], ee[6],
                                      ee[4], ee[4], ee[7], ee[7]) )
                           + " " + mat_ids[ii] + "\n" )
        for ee in el['pyram5']:
            ii = int( ee[0] ) - 1
            file_out.write( " ".join( (ee[5], ee[3], ee[4], ee[6],
                                      ee[7], ee[7], ee[7], ee[7]) )
                           + " " + mat_ids[ii] + "\n" )
else:
    if (len( el['tria3'] ) > 0):
        file_out.write( "Triangles\n%d\n" % len( el['tria3'] ) )
        for ee in el['tria3']:
            ii = int( ee[0] ) - 1
            file_out.write( " ".join( ee[3:] ) + " " + mat_ids[ii] + "\n" )

    if (len( el['quad4'] ) > 0):
        file_out.write( "Quadrilaterals\n%d\n" % len( el['quad4'] ) )
        for ee in el['quad4']:
            ii = int( ee[0] ) - 1
            file_out.write( " ".join( ee[3:] ) + " " + mat_ids[ii] + "\n" )
    
if (filename_out != '-'):
    file_out.close()
