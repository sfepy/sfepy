#!/usr/bin/python
# 13.01.2005, c
# 17.01.2005
# 03.02.2005
# 31.01.2006
# 29.05.2007
import sys
import fileinput

if (len( sys.argv ) == 3):
    fileNameIn = sys.argv[1];
    fileNameOut = sys.argv[2];
else:
    print 'Two args required!'
    raise ValueError

if (fileNameOut == '-'):
    fileOut = sys.stdout
else:
    fileOut = open( fileNameOut, "w" ); 

mode = 0
nod = []
els = {3 : {'tetra4' : [], 'pyram5' : [], 'wedge6' : [], 'brick8' : []},
       2 : {'tria3' : [], 'quad4' : []}}

groupIds = []
groupNEls = []
groups = []

input = fileinput.input( fileNameIn )
for line in input:
    row = line.split()
    if len( row ) == 0: continue

    if (row[0] == 'NUMNP'):
        row = input.readline().split()
        nNod, nEl, dim = row[0], row[1], int( row[4] )
        print nNod, nEl, dim
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
        groupIds.append( row[1] )
        gNEl = int( row[3] )
        groupNEls.append( gNEl )
        name = input.readline().strip()
        print name, groupIds[-1]
        
        els = []
        row = input.readline().split()
        row = input.readline().split()
        while not( row[0] == 'ENDOFSECTION' ):
            els.extend( row )
            row = input.readline().split()
        if gNEl != len( els ):
            print 'wrong number of group elements! (%d == %d)'\
                  % (nEl, len( els ))
            raise ValueError
        groups.append( els )
        
if int( nEl ) != sum( groupNEls ):
    print 'wrong total number of group elements! (%d == %d)'\
          % (int( nEl ), len( groupNEls ))

matIds = [None] * int( nEl )
for ii, els in enumerate( groups ):
    for iel in els:
        matIds[int( iel ) - 1] = groupIds[ii]

out = 'els: '
total = 0
for key, lst in el.iteritems():
    num = len( lst )
    out += ' %s: %d,' % (key, num)
    total += num
out += ' total: %d' % total
print out

fileOut.write( """MeshVersionFormatted 1
Dimension %d
""" % dim )
    
fileOut.write( "Vertices\n%d\n" % len( nod ) )
for nn in nod:
    fileOut.write( " ".join( nn ) + " 0\n" )

if dim == 3:
    if (len( el['tetra4'] ) > 0):
        fileOut.write( "Tetrahedra\n%d\n" % len( el['tetra4'] ) )
        for ee in el['tetra4']:
            ii = int( ee[0] ) - 1
            fileOut.write( " ".join( ee[3:] ) + " " + matIds[ii] + "\n" )

    nHex = len( el['pyram5'] ) + len( el['wedge6'] ) + len( el['brick8'] )
    if (nHex > 0):
        fileOut.write( "Hexahedra\n%d\n" % nHex )
        for ee in el['brick8']:
            ii = int( ee[0] ) - 1
            fileOut.write( " ".join( (ee[3], ee[4], ee[6], ee[5],
                                      ee[7], ee[8], ee[10], ee[9]) )
                           + " " + matIds[ii] + "\n" )
        for ee in el['wedge6']:
            ii = int( ee[0] ) - 1
            fileOut.write( " ".join( (ee[3], ee[5], ee[8], ee[6],
                                      ee[4], ee[4], ee[7], ee[7]) )
                           + " " + matIds[ii] + "\n" )
        for ee in el['pyram5']:
            ii = int( ee[0] ) - 1
            fileOut.write( " ".join( (ee[5], ee[3], ee[4], ee[6],
                                      ee[7], ee[7], ee[7], ee[7]) )
                           + " " + matIds[ii] + "\n" )
else:
    if (len( el['tria3'] ) > 0):
        fileOut.write( "Triangles\n%d\n" % len( el['tria3'] ) )
        for ee in el['tria3']:
            ii = int( ee[0] ) - 1
            fileOut.write( " ".join( ee[3:] ) + " " + matIds[ii] + "\n" )

    if (len( el['quad4'] ) > 0):
        fileOut.write( "Quadrilaterals\n%d\n" % len( el['quad4'] ) )
        for ee in el['quad4']:
            ii = int( ee[0] ) - 1
            fileOut.write( " ".join( ee[3:] ) + " " + matIds[ii] + "\n" )
    
if (fileNameOut != '-'):
    fileOut.close()
